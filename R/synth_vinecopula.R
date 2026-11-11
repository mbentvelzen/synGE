#' Generate synthetic data with a vine copula
#' 
#' Fits kernel density marginals ([kde1d::kde1d()]) and a vine copula
#' ([rvinecopulib::vinecop()]), then simulates new rows. Factor and character
#' columns are treated as discrete variables and mapped back to their original labels.
#'
#' @param real_df A data.frame with the real data.
#' @param seed Optional random seed.
#' @param n Number of synthetic rows.
#' @param oversample_factor Simulate `oversample_factor * n` rows, drop rows
#'   with invalid category codes, then keep the first `n`.
#' @param return_model If `TRUE`, return a list with the synthetic data and
#'   the fitted model.
#' @param ... Further arguments passed to [rvinecopulib::vinecop()].
#' @return A data.frame, or (if `return_model = TRUE`) a list with `syn_df` and `model`.
#' @examples
#' \dontrun{
#' syn_df <- synth_vc(iris, seed = 234)
#' }
#' @export


synth_vc <- function(real_df,
                     seed = NULL,
                     n = nrow(real_df),
                     oversample_factor = 1.5,
                     return_model = FALSE,
                     ...) {
  stopifnot(is.data.frame(real_df))
  if (!is.null(seed)) set.seed(seed)
  
  # 1) characters to factors, remember factor levels
  df_fac <- .chr_to_factor(real_df)
  factor_cols <- names(df_fac)[vapply(df_fac, is.factor, logical(1))]
  factor_levels <- lapply(df_fac[factor_cols], levels)
  
  # 2) factors to integer codes
  df_num <- df_fac
  df_num[factor_cols] <- lapply(df_num[factor_cols], as.integer)
  df_num <- as.data.frame(df_num)
  
  # 3) variable types: continuous ("c") or discrete ("d")
  var_types <- ifelse(names(df_num) %in% factor_cols, "d", "c")
  names(var_types) <- names(df_num)
 
  # 4) fit and simulate
  fit_vc <- .vc_estimate(df_num, var_types = var_types, ...)
  sim_df <- .vc_simulate(fit_vc, ceiling(oversample_factor * n))
  
  # 5) keep rows with valid codes and map back to labels
  if (length(factor_cols) > 0) {
    keep <- rep(TRUE, nrow(sim_df))
    for (col in factor_cols) {
      keep <- keep & sim_df[[col]] %in% seq_along(factor_levels[[col]])
    }
    sim_df <- sim_df[keep, , drop = FALSE]
    for (col in factor_cols) {
      lvls <- factor_levels[[col]]
      sim_df[[col]] <- factor(sim_df[[col]], levels = seq_along(lvls), labels = lvls)
    }
  }
  
  # 6) trim to n rows
  if (nrow(sim_df) < n) {
    warning("Fewer valid rows than requested after relabeling; returning available rows.")
  } else {
    sim_df <- sim_df[seq_len(n), , drop = FALSE]
  }
  sim_df <- as.data.frame(sim_df)
  rownames(sim_df) <- NULL
  
  if (return_model) {
    return(list(syn_df = sim_df,
                model = list(fit = fit_vc,
                             factor_cols = factor_cols,
                             factor_levels = factor_levels)))
  }
  sim_df
}

#' @rdname synth_vc
#' @description `synth_vc_model()` is the same as `synth_vc()` but returns
#'   the fitted model by default.
#' @export
synth_vc_model <- function(real_df, seed = NULL, n = nrow(real_df),
                           oversample_factor = 1.5, return_model = TRUE, ...) {
  synth_vc(real_df, seed = seed, n = n, oversample_factor = oversample_factor,
           return_model = return_model, ...)
}

# Internal: kernel density marginal 

.vc_spline_marginal <- function(covariate, xmin = NaN) {
  covariate <- covariate[!is.na(covariate)]
  param <- kde1d::kde1d(covariate, xmin = xmin)
  param$x <- NULL
  list(pdf     = function(u) kde1d::qkde1d(u, param),
       pit     = function(x) kde1d::pkde1d(x, param),
       rdist   = function(n) kde1d::rkde1d(n, param),
       density = function(x) kde1d::dkde1d(x, param),
       dist    = param)
}

# Internal: estimate the vine copula 
# polynomial = TRUE fits per-subject quadratic random effects (needs lme4).

.vc_estimate <- function(dat, variables_of_interest = NULL, polynomial = FALSE,
                         ID_name = NULL, time_name = NULL, keep_data = FALSE,
                         var_types = NULL, ...) {
  if (is.null(variables_of_interest)) {
    variables_of_interest <- setdiff(colnames(dat), c(ID_name, time_name))
  }
  
  if (polynomial) {
    if (is.null(ID_name) || is.null(time_name)) {
      stop("ID_name and time_name should be specified when polynomial = TRUE")
    }
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("Package 'lme4' is required when polynomial = TRUE.")
    }
    time_range <- range(dat[, time_name])
    dat_out <- as.data.frame(matrix(NA, nrow = length(unique(dat[, ID_name])),
                                    ncol = 3 * length(variables_of_interest)))
    names(dat_out) <- paste0("b", 0:2, "_", rep(variables_of_interest, each = 3))
    for (variable in variables_of_interest) {
      formula_v <- stats::as.formula(paste(
        variable, "~ poly(", time_name, ", 2, raw = TRUE) +",
        "(1 + poly(", time_name, ", 2, raw = TRUE)|", ID_name, ")"))
      re_lm <- lme4::lmer(formula_v, data = dat, REML = TRUE)
      dat_out[grep(paste0("_", variable), names(dat_out))] <- stats::coef(re_lm)[[ID_name]]
    }
  } else {
    dat_out <- dat[, variables_of_interest, drop = FALSE]
    time_range <- NULL
  }
  
  # marginals and probability integral transform
  marginals <- list()
  dat_unif <- dat_out
  for (variable in names(dat_unif)) {
    marginals[[variable]] <- .vc_spline_marginal(dat_out[, variable])
    ind_na <- is.na(dat_unif[, variable])
    dat_unif[!ind_na, variable] <- marginals[[variable]]$pit(dat_unif[!ind_na, variable])
  }
  
  # discrete variables need the left limit F(x - 1) as extra columns
  if (!is.null(var_types) && any(var_types == "d")) {
    discrete_vars <- names(dat_out)[var_types == "d"]
    discrete_data <- dat_out[, discrete_vars, drop = FALSE]
    for (variable in discrete_vars) {
      ind_na <- is.na(dat_out[, variable])
      discrete_data[!ind_na, variable] <- marginals[[variable]]$pit(dat_out[!ind_na, variable] - 1)
    }
    names(discrete_data) <- paste0(names(discrete_data), "_d")
    dat_unif <- cbind.data.frame(dat_unif, discrete_data)
  }
  
  vine_coefs <- if (is.null(var_types)) {
    rvinecopulib::vinecop(dat_unif, ...)
  } else {
    rvinecopulib::vinecop(dat_unif, var_types = unname(var_types), ...)
  }
  
  out <- list(vine_copula = vine_coefs, marginals = marginals,
              polynomial = polynomial,
              names = c(ID_name = ID_name, time_name = time_name),
              time_range = time_range,
              variables_of_interest = variables_of_interest)
  if (keep_data) {
    out <- append(out, list(original_data = dat_out, uniform_data = dat_unif), after = 3)
  }
  class(out) <- "estVineCopula"
  out
}

# Internal: simulate from the fitted vine copula 
# value_only = FALSE returns both the longitudinal values and the parameters.

.vc_simulate <- function(vine_output, n, value_only = TRUE) {
  vt <- vine_output$vine_copula$var_types
  is_disc <- vt == "d"
  
  if (any(is_disc)) {
    vine_dist <- rvinecopulib::vinecop_dist(vine_output$vine_copula$pair_copulas,
                                            vine_output$vine_copula$structure,
                                            var_types = vt)
    dat_sim <- as.data.frame(rvinecopulib::rvinecop(n, vine_dist))
  } else {
    dat_sim <- as.data.frame(rvinecopulib::rvinecop(n, vine_output$vine_copula))
  }
  names(dat_sim) <- names(vine_output$marginals)[seq_len(ncol(dat_sim))]
  
  for (variable in names(dat_sim)) {
    ind_na <- is.na(dat_sim[, variable])
    dat_sim[!ind_na, variable] <- vine_output$marginals[[variable]]$pdf(dat_sim[!ind_na, variable])
  }
  if (any(is_disc)) {
    dat_sim[, is_disc] <- round(dat_sim[, is_disc])
  }
  
  if (vine_output$polynomial) {
    id_name   <- vine_output$names["ID_name"]
    time_name <- vine_output$names["time_name"]
    gest_times <- seq(vine_output$time_range[1], vine_output$time_range[2], length.out = 100)
    time_data <- as.data.frame(tidyr::expand_grid(rownames(dat_sim), gest_times))
    names(time_data) <- vine_output$names
    df_sim <- suppressMessages(
      dplyr::right_join(tibble::rownames_to_column(dat_sim, id_name), time_data)
    )
    for (variable in vine_output$variables_of_interest) {
      col_ind <- grep(paste0("_", variable), names(df_sim))
      df_sim[, variable] <- df_sim[, col_ind[1]] +
        df_sim[, col_ind[2]] * df_sim[, time_name] +
        df_sim[, col_ind[3]] * df_sim[, time_name]^2
    }
    output <- df_sim[, c(as.character(vine_output$names), vine_output$variables_of_interest)]
    if (value_only) return(output)
    return(list(values = output, parameters = dat_sim))
  }
  
  dat_sim
}