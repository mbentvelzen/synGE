# Internal helpers shared by the synthesizers and evaluators.
# Functions starting with "." are internal (not exported).

# Package imports
# Functions from other packages that are called without the 'pkg::' prefix.
# devtools::document() writes these to the NAMESPACE file, so no library()
# calls are needed.

#' @import ggplot2
#' @importFrom dplyr %>% across arrange desc if_else mutate row_number select slice where
#' @importFrom tibble rownames_to_column
#' @importFrom purrr is_empty
#' @importFrom glmnet cv.glmnet assess.glmnet
#' @importFrom caret trainControl
#' @importFrom fasano.franceschini.test fasano.franceschini.test
#' @importFrom stats aggregate ks.test median model.matrix na.omit predict
#'   quantile reorder sd setNames var
#' @importFrom graphics hist
NULL

# Data preparation

# Convert all character columns to factors
.chr_to_factor <- function(df) {
  df[] <- lapply(df, function(x) if (is.character(x)) factor(x) else x)
  df
}

#' Detect variable types
#'
#' Classifies each column as "binary", "categorical" or "numeric".
#'
#' @param df A data.frame.
#' @param max_unique_for_cat Numeric columns with at most this many unique
#'   values are treated as categorical.
#' @return A data.frame with columns `variable` and `type`.
#' @export
detect_types <- function(df, max_unique_for_cat = 2) {
  types <- vapply(df, function(v) {
    n_uniq <- length(unique(stats::na.omit(v)))
    if (is.logical(v) || n_uniq <= 2) return("binary")
    if (is.numeric(v)) {
      if (n_uniq <= max_unique_for_cat) return("categorical")
      return("numeric")
    }
    "categorical"
  }, character(1))
  data.frame(variable = names(types), type = unname(types),
             stringsAsFactors = FALSE)
}

# Divergences

#' Kullback-Leibler divergence (base 2)
#' @param P,Q Probability vectors of equal length.
#' @param epsilon Small constant added to avoid log(0).
#' @export
KL_div <- function(P, Q, epsilon = 1e-12) {
  P <- P + epsilon
  Q <- Q + epsilon
  sum(P * log2(P / Q))
}

#' Jensen-Shannon divergence (base 2)
#' @inheritParams KL_div
#' @export
JS_div <- function(P, Q, epsilon = 1e-12) {
  P <- P + epsilon
  Q <- Q + epsilon
  M <- 0.5 * (P + Q)
  0.5 * KL_div(P, M, epsilon = 0) + 0.5 * KL_div(Q, M, epsilon = 0)
}

# Quick quality score (used for tuning CTGAN / TVAE)

# Total variation distance between two probability vectors
.compute_tvd <- function(p, q) sum(abs(p - q)) / 2

# TVD between the category frequencies of two vectors
.evaluate_tvd_cat <- function(real_vec, synth_vec) {
  lv <- union(unique(as.character(real_vec)), unique(as.character(synth_vec)))
  p  <- prop.table(table(factor(as.character(real_vec),  levels = lv)))
  q  <- prop.table(table(factor(as.character(synth_vec), levels = lv)))
  .compute_tvd(as.numeric(p), as.numeric(q))
}

# Kolmogorov-Smirnov statistic between two numeric vectors
.evaluate_ks <- function(real_vec, synth_vec) {
  as.numeric(suppressWarnings(stats::ks.test(real_vec, synth_vec)$statistic))
}

# Mean of KS (numeric columns) and TVD (discrete columns); lower is better
.evaluate_quality <- function(real_df, synth_df, numeric_cols, discrete_cols) {
  ks  <- vapply(numeric_cols,  function(col) .evaluate_ks(real_df[[col]], synth_df[[col]]), numeric(1))
  tvd <- vapply(discrete_cols, function(col) .evaluate_tvd_cat(real_df[[col]], synth_df[[col]]), numeric(1))
  mean(c(if (length(ks)) mean(ks), if (length(tvd)) mean(tvd)), na.rm = TRUE)
}

# Python (CTGAN / TVAE) helpers

# Import the Python 'ctgan' module with a clear error if it is missing
.import_ctgan <- function() {
  if (!reticulate::py_module_available("ctgan")) {
    stop("Python module 'ctgan' not found. Install it with ",
         "reticulate::py_install('ctgan') or select an environment that has it.",
         call. = FALSE)
  }
  reticulate::import("ctgan")
}

# Set the seed in R, NumPy, Python 'random' and torch
.set_seed_all <- function(seed) {
  if (is.null(seed)) return(invisible(NULL))
  set.seed(seed)
  reticulate::import("numpy")$random$seed(as.integer(seed))
  reticulate::import("random")$seed(as.integer(seed))
  if (reticulate::py_module_available("torch")) {
    reticulate::import("torch")$manual_seed(as.integer(seed))
  }
  invisible(NULL)
}

# Make sure a sampled table is an R data.frame
.as_r_df <- function(x) {
  if (inherits(x, "python.builtin.object")) x <- reticulate::py_to_r(x)
  as.data.frame(x, stringsAsFactors = FALSE)
}

# Tuning helpers (shared by synth_ctgan and synth_vae)

# Loop over the rows of a tuning grid, fit a model for each row and keep the best.
# fit_fun(p) must return list(synth = <data.frame>, model = <python model>).
.tune_loop <- function(tune_grid, fit_fun, real_df, numeric_cols,
                       discrete_cols, tune_eval = c("MWS", "quality")) {
  tune_eval <- match.arg(tune_eval)
  tune_results <- tune_grid
  tune_results$eval_metric <- NA_real_
  tune_results$error <- NA_character_
  best <- list(metric = Inf, params = NULL, syn = NULL, model = NULL)
  
  for (i in seq_len(nrow(tune_grid))) {
    p <- tune_grid[i, , drop = FALSE]
    res <- tryCatch({
      fit <- fit_fun(p)
      m <- if (tune_eval == "MWS") {
        eval_MWS(real_df, fit$synth, metric = "gower")$MWS
      } else {
        .evaluate_quality(real_df, fit$synth, numeric_cols, discrete_cols)
      }
      list(fit = fit, metric = m, error = NA_character_)
    }, error = function(e) list(fit = NULL, metric = NA_real_, error = conditionMessage(e)))
    
    tune_results$eval_metric[i] <- res$metric
    tune_results$error[i] <- res$error
    
    if (is.finite(res$metric) && res$metric < best$metric) {
      best <- list(metric = res$metric, params = as.list(p),
                   syn = res$fit$synth, model = res$fit)
    }
  }
  
  if (is.null(best$syn)) {
    warning("All tuning runs failed. See the 'error' column of the tuning results.")
  }
  best$tune <- tune_results
  best
}

# Build the return value after tuning
.tune_output <- function(best, return_model, return_tune) {
  if (!return_model && !return_tune) return(best$syn)
  out <- list(syn = best$syn)
  if (return_model) out$model <- best$model
  out$param <- best$params
  if (return_tune) out$tune <- best$tune
  out
}

# Distance helpers (used by DCR and NNDR)

# Prepare (scale and encode) the real and synthetic data for distance calculations
.prepare_distance_data <- function(real_df, syn_df, scale_numeric = TRUE, encoding = c("onehot", "mixed"), max_unique_for_cat = 2)
{
  encoding <- match.arg(encoding)
  
  warning("Onehot encoding results in calculating a Euclidean distance using a k-nearest neighbor searching algorithm. \n Mixed encoding uses the Gower distance/dissimilarity instead.")
  
  if (!is.data.frame(real_df) || !is.data.frame(syn_df)) {
    stop("Both real_df and syn_df must be data.frames.")
  }
  
  if (!setequal(names(real_df), names(syn_df))) {
    missing_in_syn  <- setdiff(names(real_df), names(syn_df))
    missing_in_real <- setdiff(names(syn_df), names(real_df))
    stop(paste(c("real_df and syn_df must have the same variable names. ",
                 c(if (length(missing_in_syn))  paste0("Missing in syn_df: ",  paste(missing_in_syn,  collapse=", ")),
                   if (length(missing_in_real)) paste0("Missing in real_df: ", paste(missing_in_real, collapse=", ")))),
               sep = "\n"))
  }
  
  syn_df <- syn_df[, names(real_df), drop = FALSE]
  
  if (anyNA(real_df) || anyNA(syn_df)) {
    warning("Rows with missing values are removed before distance calculations.")
    real_df <- stats::na.omit(real_df)
    syn_df  <- stats::na.omit(syn_df)
  }
  
  if (nrow(real_df) < 2) stop("real_df must contain at least 2 rows.")
  if (nrow(syn_df)  < 1) stop("syn_df must contain at least 1 row.")
  
  
  # Column types
  types_real_df <- detect_types(real_df, max_unique_for_cat)  
  types_syn_df <- detect_types(syn_df, max_unique_for_cat)  
  num_cols <- types_real_df$variable[types_real_df$type == "numeric"]
  cat_cols <- setdiff(names(real_df), num_cols)
  
  if(any(!types_real_df[match(types_real_df$variable, types_real_df$variable),]$type == types_syn_df[match(types_syn_df$variable, types_real_df$variable),]$type))  {
    warning(paste(c(paste0("real_df and syn_df have different type for variables: ", 
                           toString(types_real_df$variable[!types_real_df[match(types_real_df$variable, types_real_df$variable),]$type == types_syn_df[match(types_syn_df$variable, types_real_df$variable),]$type])))))
  }
  
  # combine to ensure consistent dummy encoding
  df_all <- rbind(real_df, syn_df)
  n_real <- nrow(real_df)
  
  for (col in cat_cols) df_all[[col]] <- as.factor(df_all[[col]])
  
  # Standardize numeric columns using statistics from real data only
  if (scale_numeric && length(num_cols) > 0) {
    mu  <- vapply(real_df[, num_cols, drop = FALSE], mean,      numeric(1))
    sdv <- vapply(real_df[, num_cols, drop = FALSE], stats::sd, numeric(1))
    sdv[sdv == 0 | is.na(sdv)] <- 1
    df_all[, num_cols] <- sweep(df_all[, num_cols, drop = FALSE], 2, mu,  "-")
    df_all[, num_cols] <- sweep(df_all[, num_cols, drop = FALSE], 2, sdv, "/")
  }
  
  if (encoding == "onehot") {
    X <- stats::model.matrix(~ . - 1, data = df_all)
    list(X_real = X[seq_len(n_real), , drop = FALSE],
         X_syn  = X[-seq_len(n_real), , drop = FALSE],
         encoding = "onehot")
  } else {
    list(X_real = df_all[seq_len(n_real), , drop = FALSE],
         X_syn  = df_all[-seq_len(n_real), , drop = FALSE],
         encoding = "mixed")
  }
}

# k nearest neighbors of each row of 'query' in 'ref'
.knn_query <- function(query, ref, k, encoding) {
  if (encoding == "onehot") {
    nn <- FNN::get.knnx(data = ref, query = query, k = k)
    list(dist = nn$nn.dist, index = nn$nn.index)
  } else {
    nn <- gower::gower_topn(x = query, y = ref, n = k)
    list(dist = t(nn$distance), index = t(nn$index))
  }
}

# k nearest neighbors of each row within 'data' (excluding the row itself)
.knn_self <- function(data, k, encoding) {
  if (encoding == "onehot") {
    nn <- FNN::get.knn(data = data, k = k)   # excludes self automatically
    list(dist = nn$nn.dist, index = nn$nn.index)
  } else {
    nn <- gower::gower_topn(x = data, y = data, n = k + 1)
    n  <- nrow(data)
    dist_mat  <- matrix(NA_real_,    nrow = n, ncol = k)
    index_mat <- matrix(NA_integer_, nrow = n, ncol = k)
    for (i in seq_len(n)) {
      keep <- which(nn$index[, i] != i)[seq_len(k)]
      dist_mat[i, ]  <- nn$distance[keep, i]
      index_mat[i, ] <- nn$index[keep, i]
    }
    list(dist = dist_mat, index = index_mat)
  }
}

# Summary table (min, quantiles, mean, max) for a named list of numeric vectors
.summary_stats <- function(values_list, digits = 3) {
  do.call(rbind, lapply(names(values_list), function(nm) {
    x <- values_list[[nm]]
    data.frame(
      measure = nm,
      n       = length(x),
      min     = round(min(x), digits),
      p05     = round(stats::quantile(x, 0.05), digits),
      q25     = round(stats::quantile(x, 0.25), digits),
      median  = round(stats::median(x),         digits),
      mean    = round(mean(x),                  digits),
      q75     = round(stats::quantile(x, 0.75), digits),
      max     = round(max(x), digits),
      row.names = NULL
    )
  }))
}

utils::globalVariables(c(
  "JSD",
  "Resample",
  "Variable",
  "abs_diff",
  "discrepancy",
  "fake",
  "label_to_show",
  "moment",
  "mtry",
  "neglog10_p",
  "obs",
  "p_label",
  "p_value",
  "pair",
  "pred",
  "real",
  "statistic",
  "syn",
  "type",
  "value",
  "var1",
  "var2",
  "variable"
))