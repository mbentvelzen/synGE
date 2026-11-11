# Multivariate evaluators: multivariate Kolmogorov-Smirnov, maximum mean
# discrepancy, PCA cosine similarity, isolation forest and Wasserstein distance.

# Multivariate Kolmogorov-Smirnov test

#' Multivariate Kolmogorov-Smirnov (Fasano-Franceschini) test
#'
#' Numeric variables are standardised with the mean and sd of the real data
#' and compared with [fasano.franceschini.test::fasano.franceschini.test()].
#' Binary and categorical variables are omitted.
#'
#' @param real_df,syn_df data.frames with the real and the synthetic data
#'   (same column names).
#' @param nPermute,seed,verbose,threads Passed to
#'   [fasano.franceschini.test::fasano.franceschini.test()].
#' @param max_unique_for_cat Numeric columns with at most this many unique
#'   values are treated as categorical (see [detect_types()]).
#' @return The test result (statistic and permutation p-value), with the
#'   names of the omitted variables in `omitted_var`.
#' @examples
#' \dontrun{
#' res_MKS <- eval_MKS(real_df, syn_df, threads = 12)
#' }
#' @export
eval_MKS <- function(real_df, syn_df, nPermute = 2000, seed = NULL, verbose = TRUE, threads = "auto", max_unique_for_cat = 5){
  if (!is.data.frame(real_df) || !is.data.frame(syn_df)) {
    stop("Both real_df and syn_df must be data.frames.")
  }
  
  if (!setequal(names(real_df), names(syn_df))) {
    missing_in_syn  <- setdiff(names(real_df), names(syn_df))
    missing_in_real <- setdiff(names(syn_df), names(real_df))
    stop(paste(c("real_df and syn_df must have the same variable names. ",
                 c(if (length(missing_in_syn))  paste0("Missing in syn_df: ",  paste(missing_in_syn,  collapse=", ")),
                   if (length(missing_in_real)) paste0("Missing in real_df: ", paste(missing_in_real, collapse=", ")))), 
               sep="\n"))
  }    
  
  types_real_df <- detect_types(real_df, max_unique_for_cat)  
  types_syn_df <- detect_types(syn_df, max_unique_for_cat)  
  
  if (sum(setNames(types_real_df$type, types_real_df$variable) != setNames(types_syn_df$type,  types_syn_df$variable))>0) {
    warning(paste(c("Error: real_df and syn_df have different type for variable: ", 
                    toString(names(setNames(types_real_df$type, types_real_df$variable)[setNames(types_real_df$type, types_real_df$variable) != setNames(types_syn_df$type,  types_syn_df$variable)])))))  
  }
  
  omitted_var <- types_real_df[!types_real_df$type == "numeric","variable" ]
  
  if (length(omitted_var)>0) {
    warning("One or multiple variables are categorical and/or binary and will be omitted from the analysis: ", 
            paste0(omitted_var, collapse = ", "))  
  }
  
  # Select only numeric variables and scale
  num_vars <- types_real_df$variable[types_real_df$type == "numeric"]
  
  real_df <- real_df[, num_vars, drop = FALSE]
  syn_df  <- syn_df[, num_vars, drop = FALSE]
  
  mu <- colMeans(real_df, na.rm = TRUE)
  sdv <- apply(real_df, 2, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  
  real_df <- sweep(real_df, 2, mu, "-")
  real_df <- sweep(real_df, 2, sdv, "/")
  syn_df  <- sweep(syn_df, 2, mu, "-")
  syn_df  <- sweep(syn_df, 2, sdv, "/")
  
  result <- fasano.franceschini.test(real_df, syn_df, seed = seed, verbose = verbose, nPermute = nPermute, threads = threads)
  
  result$omitted_var <- omitted_var
  print(result)
  # Output contains the permutation test p-value, the value of the test statistic, and the names of the omitted variables.
}

# Maximum Mean Discrepancy

#' Maximum mean discrepancy (MMD)
#'
#' All variables are one-hot encoded and scaled with the mean and sd of the
#' real data, then compared with [DataSimilarity::MMD()]. MMD = 0 and
#' p-value = 1 mean that the distributions are equal.
#'
#' @param real_df,syn_df data.frames with the real and the synthetic data
#'   (same column names).
#' @param kernel Kernel function, passed to [DataSimilarity::MMD()] as
#'   `kernel.fun`.
#' @param seed,n.perm,n.times Passed to [DataSimilarity::MMD()].
#' @return The output of [DataSimilarity::MMD()].
#' @examples
#' \dontrun{
#' res_MMD <- eval_MMD(real_df, syn_df)
#' }
#' @export
eval_MMD <- function(real_df, syn_df, kernel = 'rbfdot', seed = NULL, n.perm = 100, n.times = 1000) {
  if (!setequal(names(real_df), names(syn_df))) {
    missing_in_syn  <- setdiff(names(real_df), names(syn_df))
    missing_in_real <- setdiff(names(syn_df), names(real_df))
    stop(paste(c("real_df and syn_df must have the same variable names. ",
                 c(if (length(missing_in_syn))  paste0("Missing in syn_df: ",  paste(missing_in_syn,  collapse=", ")),
                   if (length(missing_in_real)) paste0("Missing in real_df: ", paste(missing_in_real, collapse=", ")))), 
               sep="\n"))
  }
  
  # One-hot encode all variables
  both <- rbind(
    mutate(real_df, .src = "real"),
    mutate(syn_df, .src = "syn")
  )
  X <- model.matrix(.src ~ . - 1, data = both)
  X_real <- X[both$.src == "real", , drop = FALSE]
  X_syn  <- X[both$.src == "syn",  , drop = FALSE]

    # filter 0 variance variables before scaling data.
  # scale using real-data moments only
  sds <- apply(X_real, 2, sd)
  mus <- colMeans(X_real)
  
  keep_scale <- which(sds > 0 & !is.na(sds))
  keep_const <- which(sds == 0 | is.na(sds))
  
  X_real_sc <- X_real
  X_syn_sc  <- X_syn

  if (length(keep_scale) > 0) {
    X_real_sc[, keep_scale] <- sweep(X_real[, keep_scale, drop = FALSE], 2, mus[keep_scale], "-")
    X_real_sc[, keep_scale] <- sweep(X_real_sc[, keep_scale, drop = FALSE], 2, sds[keep_scale], "/")
    
    X_syn_sc[, keep_scale] <- sweep(X_syn[, keep_scale, drop = FALSE], 2, mus[keep_scale], "-")
    X_syn_sc[, keep_scale] <- sweep(X_syn_sc[, keep_scale, drop = FALSE], 2, sds[keep_scale], "/")
  }
  
  real_df_m <- X_real_sc
  syn_df_m <- X_syn_sc
  not_scale <- c(colnames(real_df_m)[apply(real_df_m, 2, sd) ==0], colnames(syn_df_m)[apply(syn_df_m, 2, sd) ==0])
  
  if (sum(apply(real_df_m, 2, sd) ==0) + sum(apply(syn_df_m, 2, sd) ==0)>0) {warning(paste("The following variables had 0 variance:", toString(not_scale), sep = " ")) }
  
  DataSimilarity::MMD(real_df_m, syn_df_m, kernel.fun = kernel, n.perm = n.perm, n.times = n.times, seed = seed)
  # MMD = 0 and p-value = 1: distributions are equal
}

# PCA cosine similarity

#' Cosine similarity between the principal components
#'
#' All variables are one-hot encoded and scaled with the mean and sd of the
#' real data. The principal components of the real and the synthetic data
#' are then compared with the absolute cosine similarity (1 = same direction,
#' 0 = completely different).
#'
#' @param real_df,syn_df data.frames with the real and the synthetic data
#'   (same column names).
#' @param var_threshold Use the smallest number of components that explains
#'   at least this fraction of the variance in the real data.
#' @param n_comp Optional fixed number of components (overrules
#'   `var_threshold`).
#' @return A list with `k`, `var_threshold`, `var_expl_real`, `cum_var_real`,
#'   `cosine_per_comp` and `weighted_cosine` (the main summary statistic).
#' @examples
#' \dontrun{
#' res_PCA <- eval_PCA(real_df, syn_df)
#' }
#' @export
eval_PCA <- function(real_df, syn_df, var_threshold = 0.80, n_comp = NULL){
  both <- rbind(
    mutate(real_df, .src = "real"),
    mutate(syn_df, .src = "syn")
  )
  
  X <- model.matrix(.src ~ . - 1, data = both)
  X_real <- X[both$.src == "real", , drop = FALSE]
  X_syn  <- X[both$.src == "syn",  , drop = FALSE]
  
  #Remove zero variance var:
  zero_var <- c(names(which(apply(X , 2, var) == 0)), names(which(apply(X_real , 2, var) == 0)), names(which(apply(X_syn , 2, var) == 0)))
  if (!is_empty(zero_var)) {
    warning(paste0("Zero variance variables are removed: ", zero_var ))
    X <- X[,!colnames(X) %in% zero_var]
    X_real <- X[both$.src == "real", , drop = FALSE]
    X_syn  <- X[both$.src == "syn",  , drop = FALSE]
  }
  
  # center/scale using real data only
  mu_real   <- colMeans(X_real, na.rm = TRUE)
  sd_real   <- apply(X_real, 2, sd, na.rm = TRUE)
  X_real_sc <- sweep(sweep(X_real, 2, mu_real, "-"), 2, sd_real, "/")
  X_syn_sc  <- sweep(sweep(X_syn, 2, mu_real, "-"), 2, sd_real, "/")
  
  # Determine number of components
  pca0 <- svd(X_real_sc)
  var_expl_real <- pca0$d^2 / sum(pca0$d^2)
  cum_var       <- cumsum(var_expl_real)
  k             <- which(cum_var >= var_threshold)[1]
  if (is.na(k))          k <- length(var_expl_real) 
  if (!is.null(n_comp))  k <- min(n_comp, ncol(X_real_sc), ncol(X_syn_sc))# min(k, n_comp)
  
  pca1 <- svd(X_syn_sc,  nu = k, nv = k)
  
  # Absolute cosine similarity per component
  cos_per_comp <- sapply(seq_len(k), function(ii) {
    v0 <- pca0$v[, ii]; v1 <- pca1$v[, ii]
    abs(sum(v0 * v1) / sqrt(sum(v0^2) * sum(v1^2)))
  })
  
  # Variance-weighted average of the cosine similarity, weights = real eigenvalues (d^2)
  lambda_real     <- pca0$d[seq_len(k)]^2
  weighted_cosine <- sum(lambda_real * cos_per_comp) / sum(lambda_real)
  
  list(
    k                = k,
    var_threshold    = var_threshold,
    var_expl_real    = var_expl_real[seq_len(k)],
    cum_var_real     = cum_var[seq_len(k)],
    cosine_per_comp  = cos_per_comp,  
    weighted_cosine  = weighted_cosine  # main summary statistic
  )
  
} # 1 component is equal; -1 component is equal but flipped sign; ~0 component is completely different

# Anomaly detection (isolation forest)

#' Compare anomaly scores with an isolation forest
#'
#' Fits an isolation forest ([isotree::isolation.forest()]) on the real data
#' and compares the outlier scores and distances of the real and the
#' synthetic data.
#'
#' @param real_df,syn_df data.frames with the real and the synthetic data
#'   (same column names).
#' @param detect_type If `TRUE`, binary and categorical columns are detected
#'   with [detect_types()] and passed on as `categ_cols`.
#' @param max_unique_for_cat Numeric columns with at most this many unique
#'   values are treated as categorical (see [detect_types()]).
#' @param ndim,ntry,scoring_metric,ntrees,sample_size,max_depth,missing_action,categ_split_type,recode_categ,penalize_range,seed,...
#'   Passed to [isotree::isolation.forest()].
#' @return A list with the fitted `model`, summaries, Kolmogorov-Smirnov
#'   tests and outlier rates (threshold: 95th percentile of the real scores),
#'   and the scores and distances of the real and the synthetic data.
#' @examples
#' \dontrun{
#' res_IF <- eval_IF(real_df, syn_df)
#' eval_IF_plot(res_IF)
#' }
#' @export
eval_IF <- function(real_df, syn_df, detect_type = FALSE, max_unique_for_cat = 2,
                    ndim = 2, ntry = 10, scoring_metric = "depth", ntrees = 100,
                    sample_size = 500, max_depth = 10,
                    missing_action = "fail", categ_split_type = "subset",
                    recode_categ = TRUE, penalize_range = TRUE,
                    seed = 12345, ... ) {
  
  if (!is.data.frame(real_df) || !is.data.frame(syn_df)) {
    stop("Both real_df and syn_df must be data.frames.")
  }
  
  if (!setequal(names(real_df), names(syn_df))) {
    missing_in_syn  <- setdiff(names(real_df), names(syn_df))
    missing_in_real <- setdiff(names(syn_df), names(real_df))
    stop(paste(c("real_df and syn_df must have the same variable names. ",
                 c(if (length(missing_in_syn))  paste0("Missing in syn_df: ",  paste(missing_in_syn,  collapse=", ")),
                   if (length(missing_in_real)) paste0("Missing in real_df: ", paste(missing_in_real, collapse=", ")))), 
               sep="\n"))
  }
  
  vars <- names(real_df)
  syn_df <- syn_df[vars]
  
  if (detect_type) {
    types_real_df <- detect_types(real_df, max_unique_for_cat)
    types_syn_df <- detect_types(syn_df, max_unique_for_cat)
    
    if (sum(setNames(types_real_df$type, types_real_df$variable) != setNames(types_syn_df$type,  types_syn_df$variable))>0) {
      warning(paste(c("Error: real_df and syn_df have different type for variable: ", 
                      toString(names(setNames(types_real_df$type, types_real_df$variable)[setNames(types_real_df$type, types_real_df$variable) != setNames(types_syn_df$type,  types_syn_df$variable)])))))  
    }
    
    categ_cols <- types_real_df$variable[types_real_df$type %in% c("binary", "categorical")]
  } else {
    categ_cols <- NULL
  }
  
  set.seed(seed)
  
  model <- isotree::isolation.forest(ntrees = ntrees,
                                     real_df,
                                     seed = seed,
                                     ndim = ndim,
                                     ntry = ntry,
                                     sample_size = sample_size, 
                                     max_depth = max_depth,
                                     scoring_metric = scoring_metric,
                                     missing_action = missing_action,
                                     categ_split_type = categ_split_type,
                                     recode_categ = recode_categ,
                                     penalize_range = penalize_range,
                                     categ_cols = categ_cols, ...
  )
  
  scores_real <- predict(model, real_df, type = "score")
  scores_syn  <- predict(model, syn_df,  type = "score")
  dist_real   <- predict(model, real_df, type = "dist")
  dist_syn    <- predict(model, syn_df,  type = "dist")
  
  # Outlier threshold (95th percentile of real data)
  threshold <- quantile(scores_real, 0.95, na.rm = TRUE)
  real_outlier_rate <- mean(scores_real > threshold, na.rm = TRUE)
  syn_outlier_rate  <- mean(scores_syn  > threshold, na.rm = TRUE)
  
  # KS test to compare score distributions
  ks_scores <- ks.test(scores_real, scores_syn)
  ks_dist   <- ks.test(dist_real, dist_syn)
  
  # Summary statistics  
  summ <- function(x) {
    c(mean   = mean(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      p95    = quantile(x, 0.95, na.rm = TRUE),
      max    = max(x, na.rm = TRUE)
    )
  }
  
  score_summary <- rbind(
    Real = summ(scores_real),
    Synthetic = summ(scores_syn)
  )
  
  dist_summary <- rbind(
    Real = summ(dist_real),
    Synthetic = summ(dist_syn)
  )
  
  # Output
  list(
    model           = model,
    score_summary   = score_summary,
    dist_summary    = dist_summary,
    ks_score_test   = ks_scores,
    ks_dist_test    = ks_dist,
    threshold_95    = threshold,
    outlier_rate_real = real_outlier_rate,
    outlier_rate_syn  = syn_outlier_rate,
    scores_real     = scores_real,
    scores_syn      = scores_syn,
    dist_real       = dist_real,
    dist_syn        = dist_syn,
    categ_cols_used = categ_cols,
    columns_used    = vars,
    note            = "Isolation Forest anomaly comparison completed."
  )
}

#' Plot the isolation forest results
#'
#' @param res_IF Output of [eval_IF()].
#' @param plot `"all"`, `"density"`, `"boxplot"`, `"violin"` or `"text"`.
#' @param metric `"score"` (outlier scores) or `"dist"` (distances).
#' @return A list with the requested plots and/or the text summary.
#' @export
eval_IF_plot <- function(res_IF, 
                         plot = c("all", "density", "boxplot", "violin", "text"),
                         metric = c("score", "dist")) {
  plot <- match.arg(plot)
  metric <- match.arg(metric)
  
  # Extract selected values
  if (metric == "score") {
    real_vals <- res_IF$scores_real
    syn_vals  <- res_IF$scores_syn
    plot_title <- "Isolation Forest Outlier Scores"
  } else {
    real_vals <- res_IF$dist_real
    syn_vals  <- res_IF$dist_syn
    plot_title <- "Isolation Forest Pairwise Distances"
  }
  
  df <- data.frame(
    value = c(real_vals, syn_vals),
    type  = factor(rep(c("Real", "Synthetic"),
                       c(length(real_vals), length(syn_vals))))
  )
  
  density_plot <- ggplot2::ggplot(df, aes(x = value, fill = type)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0(plot_title, " - Density Plot"),
                  x = metric, y = "Density")
  
  box_plot <- ggplot2::ggplot(df, aes(x = type, y = value, fill = type)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0(plot_title, " - Boxplot"),
                  x = "", y = metric)
  
  violin_plot <- ggplot2::ggplot(df, aes(x = type, y = value, fill = type)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.6) +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0(plot_title, " - Violin Plot"),
                  x = "", y = metric)
  
  text_out <- NULL
  
  # Outlier-rate
  or_real <- round(res_IF$outlier_rate_real * 100, 2)
  or_syn  <- round(res_IF$outlier_rate_syn * 100, 2)
  
  # KS test interpretation
  ks_stat <- if (metric == "score") res_IF$ks_score_test$statistic else res_IF$ks_dist_test$statistic
  ks_pval <- if (metric == "score") res_IF$ks_score_test$p.value else res_IF$ks_dist_test$p.value
  ks_pval_r <- signif(ks_pval, 3)
  
  text_out <- paste0(
    "KS test statistic of Synthetic vs Real (metric = ", metric, "):\n",
    "- The KS test statistic is ", round(ks_stat, 3),
    " with p-value = ", ks_pval_r, "\n",
    if (ks_pval < 0.05) {
      "  \u2192 The distributions differ significantly.\n"
    } else {
      "  \u2192 No significant difference in distributions detected.\n"
    },
    "Outlier rate (threshold is the 95th percentile of the scores in the real data).\n",
    "- Real-data outlier rate: ", paste0(or_real, "%"), "\n",
    "- Synthetic-data outlier rate: ", paste0(or_syn, "%"), "\n",
    if (or_syn > or_real * 2) {
      "  \u2192 Synthetic data produces substantially more anomalies.\n\n"
    } else if (or_syn < or_real * 0.5) {
      "  \u2192 Synthetic data appears over-smoothed (too few anomalies).\n\n"
    } else {
      "  \u2192 Synthetic and real anomaly rates are comparable.\n"
    },
    sep = ""
  )
  
  output <- list()
  
  if (plot %in% c("density", "all")) output$density <- density_plot
  if (plot %in% c("boxplot", "all")) output$boxplot <- box_plot
  if (plot %in% c("violin", "all")) output$violin <- violin_plot
  if (plot %in% c("text", "all"))   output$text <- text_out
  
  # If only 1 item requested, print directly
  if (plot != "all" && length(output) == 1) {
    print(output[[1]])
  }
  
  return(output)
}

# Wasserstein distance

#' Wasserstein distance between real and synthetic data
#'
#' Computes the optimal transport (Wasserstein) distance between the rows of
#' the real and synthetic data, using a Gower or Euclidean cost.
#' Also used by [synth_ctgan()] and [synth_vae()] to pick the best tuning run.
#'
#' @param real_df,syn_df data.frames with the same column names.
#' @param metric `"gower"` or `"euclidean"`.
#' @param p Order of the distance; the cost is ||x - y||^p.
#' @param scale_numeric Scale numeric columns (Euclidean only).
#' @param weights Optional named numeric vector of column weights (Gower only).
#' @param max_n Optional maximum number of rows per dataset (subsampled).
#' @param seed Seed used for subsampling.
#' @param method Solver passed to [transport::transport()].
#' @param verbose Print the size of the cost matrix.
#' @return A one row data.frame with the distance in column `MWS`.
#' @export
eval_MWS <- function(real_df,
                     syn_df,
                     metric = c("gower", "euclidean"),
                     p = 1,                      # cost is ||x-y||^p; the returned value is (sum mass*cost)^(1/p)
                     scale_numeric = TRUE,
                     weights = NULL,             # optional column weights (named numeric vector)
                     max_n = Inf,                # optional subsample cap per dataset (e.g., 2000)
                     seed = 1,
                     method = c("shortsimplex", "primaldual", "revsimplex"),
                     verbose = FALSE) {
  
  metric <- match.arg(metric)
  method <- match.arg(method)
  
  if (!is.data.frame(real_df) || !is.data.frame(syn_df)) {
    stop("Both real_df and syn_df must be data.frames.")
  }
  
  if (!setequal(names(real_df), names(syn_df))) {
    missing_in_syn  <- setdiff(names(real_df), names(syn_df))
    missing_in_real <- setdiff(names(syn_df), names(real_df))
    stop(paste(
      "real_df and syn_df must have the same variable names.",
      if (length(missing_in_syn))  paste0("Missing in syn_df: ",  paste(missing_in_syn,  collapse=", ")),
      if (length(missing_in_real)) paste0("Missing in real_df: ", paste(missing_in_real, collapse=", ")),
      sep = "\n"
    ))
  }
  
  # keep same column order
  syn_df <- syn_df[, names(real_df), drop = FALSE]
  
  # optional subsampling to control n^2
  set.seed(seed)
  if (is.finite(max_n)) {
    if (nrow(real_df) > max_n) real_df <- real_df[sample.int(nrow(real_df), max_n), , drop = FALSE]
    if (nrow(syn_df)  > max_n) syn_df  <- syn_df[sample.int(nrow(syn_df),  max_n), , drop = FALSE]
  }
  
  nR <- nrow(real_df)
  nS <- nrow(syn_df)
  if (nR < 1 || nS < 1) stop("Both data.frames must have at least 1 row.")
  
  # helpers
  
  # make factor levels consistent across datasets for factor/character columns
  harmonize_levels <- function(r, s) {
    for (nm in names(r)) {
      if (is.character(r[[nm]]) || is.factor(r[[nm]]) ||
          is.character(s[[nm]]) || is.factor(s[[nm]])) {
        
        r_chr <- as.character(r[[nm]])
        s_chr <- as.character(s[[nm]])
        levs <- sort(unique(c(r_chr[!is.na(r_chr)], s_chr[!is.na(s_chr)])))
        
        r[[nm]] <- factor(r_chr, levels = levs)
        s[[nm]] <- factor(s_chr, levels = levs)
      }
    }
    list(r = r, s = s)
  }
  
  # build numeric feature matrices for euclidean option (one-hot + scaling)
  build_euclidean_features <- function(r, s, scale_numeric = TRUE) {
    hs <- harmonize_levels(r, s)
    r <- hs$r; s <- hs$s
    
    # identify numeric columns
    is_num <- vapply(r, is.numeric, logical(1))
    num_names <- names(r)[is_num]
    
    # optional scaling for numeric columns (based on real_df stats; fallback to combined)
    if (scale_numeric && length(num_names) > 0) {
      for (nm in num_names) {
        mu <- mean(r[[nm]], na.rm = TRUE)
        sdv <- stats::sd(r[[nm]], na.rm = TRUE)
        if (!is.finite(sdv) || sdv == 0) {
          # fallback: use combined sd
          comb <- c(r[[nm]], s[[nm]])
          sdv2 <- stats::sd(comb, na.rm = TRUE)
          if (is.finite(sdv2) && sdv2 > 0) sdv <- sdv2 else sdv <- 1
        }
        r[[nm]] <- (r[[nm]] - mu) / sdv
        s[[nm]] <- (s[[nm]] - mu) / sdv
      }
    }
    
    # replace NA:
    # - numeric: 0 after scaling (roughly mean-imputed)
    # - factors: add explicit "NA" level
    for (nm in names(r)) {
      if (is.numeric(r[[nm]])) {
        r[[nm]][is.na(r[[nm]])] <- 0
        s[[nm]][is.na(s[[nm]])] <- 0
      } else if (is.factor(r[[nm]])) {
        # add NA as level, then set missing
        levs <- levels(r[[nm]])
        if (!("NA" %in% levs)) {
          levels(r[[nm]]) <- c(levs, "NA")
          levels(s[[nm]]) <- c(levs, "NA")
        }
        r[[nm]][is.na(r[[nm]])] <- "NA"
        s[[nm]][is.na(s[[nm]])] <- "NA"
      }
    }
    
    # one-hot encode factors/characters via model.matrix
    Xr <- stats::model.matrix(~ . - 1, data = r)
    Xs <- stats::model.matrix(~ . - 1, data = s)
    
    list(Xr = Xr, Xs = Xs)
  }
  
  # build cross cost matrix with Gower (mixed types)
  build_gower_cost <- function(r, s, weights = NULL) {
    if (!requireNamespace("cluster", quietly = TRUE)) {
      stop("Package 'cluster' is required for metric='gower'. Install it via install.packages('cluster').")
    }
    
    hs <- harmonize_levels(r, s)
    r <- hs$r; s <- hs$s
    
    # cluster::daisy handles NA; for factors, ensure they are factors (done above)
    combined <- rbind(r, s)
    
    # weights: can be NULL or numeric vector aligned to columns
    w <- NULL
    if (!is.null(weights)) {
      if (is.null(names(weights))) {
        if (length(weights) != ncol(combined)) stop("If weights is unnamed, it must have length = ncol(data).")
        w <- as.numeric(weights)
      } else {
        if (!all(names(combined) %in% names(weights))) {
          miss <- setdiff(names(combined), names(weights))
          stop(paste0("weights missing for columns: ", paste(miss, collapse = ", ")))
        }
        w <- as.numeric(weights[names(combined)])
      }
    } else { w <- rep.int(1, ncol(combined))}
    
    d <- cluster::daisy(combined, metric = "gower", weights = w)
    D <- as.matrix(d)
    
    # cross block: real rows vs syn rows
    C <- D[seq_len(nR), nR + seq_len(nS), drop = FALSE]
    C
  }
  
  # build cross cost matrix with Euclidean (after encoding)
  build_euclidean_cost <- function(r, s, scale_numeric = TRUE, p = 1) {
    if (!requireNamespace("proxy", quietly = TRUE)) {
      stop("Package 'proxy' is required for metric='euclidean'. Install it via install.packages('proxy').")
    }
    feats <- build_euclidean_features(r, s, scale_numeric = scale_numeric)
    Xr <- feats$Xr
    Xs <- feats$Xs
    
    # proxy::dist can compute cross distances:
    # dist(x, y) returns an nrow(x) x nrow(y) "dist" object; coerce to matrix
    C <- as.matrix(proxy::dist(Xr, Xs, method = "Euclidean"))
    
    if (p != 1) C <- C^p
    C
  }
  
  # solve OT and compute Wasserstein value
  solve_ot <- function(C, p = 1, method = "shortsimplex") {
    if (!requireNamespace("transport", quietly = TRUE)) {
      stop("Package 'transport' is required. Install it via install.packages('transport').")
    }
    a <- rep(1 / nrow(C), nrow(C))   # uniform weights for real
    b <- rep(1 / ncol(C), ncol(C))   # uniform weights for syn
    
    plan <- transport::transport(a, b, costm = C, method = method)
    
    # plan has columns: from, to, mass
    # cost for each moved mass:
    moved_cost <- C[cbind(plan$from, plan$to)]
    val_p <- sum(plan$mass * moved_cost)
    
    # return W_p (i.e., (E[cost])^(1/p)) when using ||.||^p
    if (p == 1) val_p else val_p^(1 / p)
  }
  
  #  main
  
  if (metric == "gower") {
    C <- build_gower_cost(real_df, syn_df, weights = weights)
    # Gower distance is already a proper distance; interpret p as exponent if you want
    if (p != 1) C <- C^p
  } else {
    C <- build_euclidean_cost(real_df, syn_df, scale_numeric = scale_numeric, p = p)
  }
  
  if (verbose) {
    message(sprintf("Cost matrix: %d x %d (metric=%s, p=%s)", nrow(C), ncol(C), metric, p))
  }
  
  mws <- solve_ot(C, p = p, method = method)
  
  out <- data.frame(
    metric = metric,
    p = p,
    n_real = nR,
    n_syn = nS,
    MWS = mws,
    note = if (metric == "gower")
      "Mixed-type OT using Gower cross-distance cost matrix."
    else
      "OT using Euclidean distance on one-hot + (optional) scaled numeric features.",
    stringsAsFactors = FALSE
  )
  
  out
}
