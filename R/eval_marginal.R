# Marginal (per variable) evalutation: moments, Jensen-Shannon divergence,
# Kolmogorov-Smirnov, Cramer-von Mises, Anderson-Darling and pMSE.

# Moments: mean, var, sd, skewness, kurtosis

#' Compare the moments of real and synthetic data
#'
#' Compares the mean, variance, standard deviation, skewness and (excess)
#' kurtosis of every numeric variable. Non-numeric variables are omitted.
#'
#' @param real_df,syn_df data.frames with the real and the synthetic data
#'   (make sure both have the same column names).
#' @param moment `"all"`, or one or more of `"mean"`, `"var"`, `"sd"`,
#'   `"skewness"` and `"kurtosis"`.
#' @param na_rm Remove missing values before computing the moments.
#' @return A list with `by_moment` (one row per variable and moment) and
#'   `meta`. The column `discrepancy` is the absolute difference (mean,
#'   skewness, kurtosis) or the absolute log2 ratio (var, sd); 0 is best.
#' @examples
#' \dontrun{
#' res_moments <- eval_moments(real_df, syn_df, moment = "all", na_rm = TRUE)
#' eval_moments_heatmap(res_moments, top_var = 10)
#' eval_moments_calplot(res_moments)
#' }
#' @export
eval_moments <- function(real_df, syn_df, moment = "all", na_rm = TRUE) {
  ## - basic checks - ##
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
  
  names_no_num_real <- names(real_df)[!vapply(real_df, is.numeric, logical(1))]
  names_no_num_syn  <- names(syn_df)[!vapply(syn_df,  is.numeric, logical(1))]
  if (length(c(names_no_num_real, names_no_num_syn) ) > 0) {
    warning(paste0("One or multiple variables are categorical and/or binary and will be omitted from the analysis: ", paste(unique(c(names_no_num_real, names_no_num_syn)), collapse = ", ")))
    real_df <- real_df[vapply(real_df, is.numeric, logical(1))]
    syn_df  <- syn_df[vapply(syn_df, is.numeric, logical(1))]
  }
  
  # re-order synthetic data to have the same order as the real data
  syn_df <- syn_df[, names(real_df), drop = FALSE]
  rows_real <- nrow(real_df)
  rows_syn  <- nrow(syn_df)  
  
  moments <- c("mean","var","sd","skewness","kurtosis")
  if (!all(moment %in% c(moments, "all"))) {
    stop("moment must be 'all' or one/multiple of: 'mean','var','sd','skewness','kurtosis'.")
  } 
  if (identical(moment, "all")) { moment <- moments }
  
  ## compute variables ##
  result_list <- lapply(names(real_df), function(var) {
    xr <- real_df[[var]]
    xs <- syn_df[[var]]
    
    mr <- mean(xr, na.rm = na_rm)
    ms <- mean(xs, na.rm = na_rm)
    diff_mean <- ms - mr
    
    vr <- stats::var(xr, na.rm = na_rm)
    sdr <- sqrt(vr)
    vs <- stats::var(xs, na.rm = na_rm)
    sds <- sqrt(vs)
    diff_sd <- sds - sdr
    
    skr <- e1071::skewness(xr, type = 1, na.rm = na_rm) #g1 = m3/m2^(3/2)
    sks <- e1071::skewness(xs, type = 1, na.rm = na_rm)
    diff_sk <- sks - skr
    
    # Excess Kurtosis: normal distribution has a kurtosis of 0.
    kur_r <- e1071::kurtosis(xr, type = 1, na.rm = na_rm) # g2 = m4/m2^2 - 3
    kur_s <- e1071::kurtosis(xs, type = 1, na.rm = na_rm)
    diff_kur <- kur_s - kur_r
    
    ## smd: standardised mean difference; 0 is best ##
    sp <- sqrt(((rows_real - 1) * vr + (rows_syn - 1) * vs) / (rows_real + rows_syn - 2))
    smd <- if (is.finite(sp) && sp > 0) (ms - mr) / sp else NA_real_
    
    ##  lvr and lsr: log2 variance and sd ratio ##
    lvr <- if (is.finite(vr) && vr > 0 && is.finite(vs) && vs > 0) log2(vs/vr) else NA_real_
    lsr <- if (is.finite(sdr) && sdr > 0 && is.finite(sds) && sds > 0) log2(sds/sdr) else NA_real_
    
    out <- list(
      mean      = data.frame(variable=var, moment="mean", real=mr, syn=ms,
                             diff=diff_mean, ratio=NA_real_, smd=smd, discrepancy = abs(diff_mean)),
      var       = data.frame(variable=var, moment="var", real=vr, syn=vs,
                             diff=vs-vr, ratio=lvr,smd=NA_real_, discrepancy = abs(lvr)),
      sd        = data.frame(variable=var, moment="sd",real=sdr, syn=sds,
                             diff=diff_sd, ratio=lsr, smd=NA_real_, discrepancy = abs(lsr)),
      skewness  = data.frame(variable=var, moment="skewness", real=skr, syn=sks,
                             diff=diff_sk, ratio=NA_real_, smd=NA_real_, discrepancy = abs(diff_sk)),
      kurtosis  = data.frame(variable=var, moment="kurtosis", real=kur_r, syn=kur_s,
                             diff=diff_kur, ratio=NA_real_, smd=NA_real_, discrepancy = abs(diff_kur))
    )
    
    # keep only indicated metrics
    keep_metrics <- intersect(names(out), moment)
    do.call(rbind, out[keep_metrics])
  }   )
  
  results_by_moment <- do.call(rbind, result_list)
  rownames(results_by_moment) <- NULL
  

  list(
    by_moment   = results_by_moment,
    meta        = list(n_real = rows_real, n_syn = rows_syn, moments_included = unique(results_by_moment$moment))
  )
}

#' Heatmap of the moment discrepancies
#'
#' @param res_moments Output of [eval_moments()].
#' @param top_var Number of variables to show (those with the largest
#'   standardised mean difference).
#' @return A ggplot object.
#' @export
eval_moments_heatmap <- function(res_moments, top_var = 10) {
  if (top_var > length(unique(res_moments$by_moment$variable))  || is.na(top_var) ) {
    warning("Number of 'top_var' is ", top_var, " and the number of variables is n = ", length(unique(res_moments$by_moment$variable)), ".\n",
            "top_var is set to n = ", length(unique(res_moments$by_moment$variable))) 
    top_var <- length(unique(res_moments$by_moment$variable))
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("ggnewscale", quietly = TRUE)) {
    stop('Please install.packages(c("ggplot2","ggnewscale"))')
  }
  
  df <- res_moments$by_moment
  # order variables by worst (largest) standardized difference in mean
  ord <- aggregate(abs(smd) ~ variable, df[df$moment == "mean",], max)
  ord <- ord[order(ord$`abs(smd)`, decreasing = TRUE), "variable"]
  
  df$variable <- factor(df$variable, levels = ord)
  df <- subset(df, variable %in% ord[1:top_var])
  
  df$moment <- factor(df$moment, levels = c("mean","sd","var","skewness","kurtosis"))
  
  lab_map <- c(mean = "\u0394 mean", sd = "ratio sd", var = "ratio var", skewness = "\u0394 skewness", kurtosis = "\u0394 kurtosis"  )
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = moment, y = variable)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(x = "Moment", y = "Variable",
                  title = "Real vs Synthetic: Moment Discrepancies (lower is better)")
  
  for (m in levels(df$moment)) {
    p <- p +
      ggplot2::geom_tile(
        data = df[df$moment == m,],
        ggplot2::aes(fill = discrepancy)
      ) +
      ggplot2::scale_fill_viridis_c(name = lab_map[m], option = "D") +
      ggnewscale::new_scale_fill()
  }
  p
}

#' Calibration plot of the moments (real against synthetic)
#'
#' @param res_moments Output of [eval_moments()].
#' @return A ggplot object.
#' @export
eval_moments_calplot <- function(res_moments) {
  ggplot2::ggplot(
    res_moments$by_moment,
    ggplot2::aes(x = .data[["real"]], y = .data[["syn"]])
  ) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, linetype = 2
    ) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::facet_wrap(~ moment, scales = "free") +
    ggplot2::labs(
      x = "Real",
      y = "Synthetic",
      title = "Calibration of moments"
    )
}
# Jensen-Shannon divergence

#' Jensen-Shannon divergence per variable
#'
#' Numeric variables are binned with a histogram; binary and categorical
#' variables use their category frequencies. The variable types come from
#' [detect_types()].
#'
#' @param real_df,syn_df data.frames with the real and the synthetic data
#'   (same column names).
#' @param max_unique_for_cat Numeric columns with at most this many unique
#'   values are treated as categorical (see [detect_types()]).
#' @param num_bins Number of histogram break points for numeric variables.
#'   If `NA`, `max(10, floor(n / 4))` is used, with `n` the smallest sample
#'   size.
#' @param epsilon Small constant added to avoid log(0).
#' @return A data.frame with one row per variable and the columns
#'   `variable`, `type`, `JSD`, `KLD` and `note`.
#' @examples
#' \dontrun{
#' res_JSD <- eval_JSD(real_df, syn_df)
#' eval_JSD_plot(res_JSD)
#' }
#' @export
eval_JSD <- function(real_df, syn_df, max_unique_for_cat = 5, num_bins = NA, epsilon = 1e-12){ 
  if (!is.data.frame(real_df) || !is.data.frame(syn_df)) {
    stop("Both real_df and syn_df must be data.frames.")
  }
  
  if (!setequal(names(real_df), names(syn_df))) {
    missing_in_syn  <- setdiff(names(real_df), names(syn_df))
    missing_in_real <- setdiff(names(syn_df), names(real_df))
    warning(paste(c("real_df and syn_df must have the same variable names. ",
                    c(if (length(missing_in_syn))  paste0("Missing in syn_df: ",  paste(missing_in_syn,  collapse=", ")),
                      if (length(missing_in_real)) paste0("Missing in real_df: ", paste(missing_in_real, collapse=", ")))), 
                  sep="\n"))
  }
  
  vars <- names(real_df)
  
  types_real_df <- detect_types(real_df, max_unique_for_cat)  
  types_syn_df <- detect_types(syn_df, max_unique_for_cat)  
  
  results <- lapply(vars, function(var) {    ## CHECK: LAPPLY OR SAPPLY?
    var_real <- real_df[[var]] 
    var_syn <- syn_df[[var]]
    combined <- c(var_real, var_syn)
    
    if (!(types_real_df[types_real_df$variable == var,]$type == types_syn_df[types_syn_df$variable == var,]$type)) {
      warning(paste(c(paste0("Error: real_df and syn_df have different type for variable: ", paste(var))),
                    c(paste0("real_df type is ", paste(types_real_df[types_real_df$variable == var,]$type))),
                    c(paste0("syn_df type is ", paste(types_syn_df[types_syn_df$variable == var,]$type))), sep = "\n"))
    }
    
    type <- types_real_df[types_real_df$variable == var,]$type
    
    if (type == "numeric") {  
      if (is.na(num_bins)) { 
        num_bins <- max(10, floor(min(length(var_real), length(var_syn)) / 4)) 
      } else { 
        num_bins <- as.integer(num_bins) }
      
      
      breaks <- seq(min(combined), max(combined), length.out = num_bins)
      hist_r <- hist(var_real, breaks = breaks, plot = F)
      hist_s <- hist(var_syn, breaks = breaks, plot = F)
      
      P <- hist_r$density / sum(hist_r$density)
      Q <- hist_s$density / sum(hist_s$density)
      
      note <- "numeric with histogram for estimate distribution"
      
    } else if (type == "binary") {
      if (is.logical(c(var_real, var_syn))) {
        var_real <- as.numeric(var_real)
        var_syn  <- as.numeric(var_syn)
      } else if (is.numeric(var_real) & is.numeric(var_syn)) {
        value_0  <- sort(unique(na.omit(c(var_real, var_syn))))[1]
        var_real <- ifelse(var_real == value_0, 0, 1)
        var_syn  <- ifelse(var_syn == value_0, 0, 1)
      } else {
        real_chr <- as.character(var_real)
        syn_chr  <- as.character(var_syn)
        categories <- sort(unique(c(real_chr[!is.na(real_chr)], syn_chr[!is.na(syn_chr)])))
        var_real <- ifelse(var_real == categories[1], 0, 1)
        var_syn  <- ifelse(var_syn == categories[1], 0, 1)
      }
      
      var_real_0 <- sum(var_real == 0, na.rm = TRUE); var_real_1 <- sum(var_real == 1, na.rm = TRUE)
      var_syn_0 <- sum(var_syn == 0, na.rm = TRUE); var_syn_1 <- sum(var_syn == 1, na.rm = TRUE)
      
      P <- c(var_real_0, var_real_1) / (var_real_0 + var_real_1)
      Q <- c(var_syn_0, var_syn_1) / (var_syn_0 + var_syn_1)
      
      note <- "Binary distribution"
      
    } else if (type == "categorical") {
      var_real <- as.character(var_real)
      var_syn <- as.character(var_syn)
      
      categories <- unique(c(na.omit(var_real), na.omit(var_syn)))
      
      var_real_c <- as.numeric(table(factor(var_real, levels = categories), useNA = "no"))
      var_syn_c <- as.numeric(table(factor(var_syn, levels = categories), useNA = "no"))
      
      P <- var_real_c / sum(var_real_c)
      Q <- var_syn_c / sum(var_syn_c)
      
      note <- "Categorical variable"
      
    } else { note <- paste0("unknown_type: ", type)} 
    
    
    JSD_value <- JS_div(P, Q, epsilon = epsilon)
    
    data.frame(variable = var, type = type, JSD = JSD_value, KLD = KL_div(P, Q), note = note, stringsAsFactors = FALSE) 
  }  
  )
  
  res_df <- do.call(rbind, results)
}

#' Bar plot of the Jensen-Shannon divergence per variable
#'
#' @param res_JSD Output of [eval_JSD()].
#' @param top_var Number of variables (highest JSD first) to show. `NA`
#'   shows all variables.
#' @return A ggplot object.
#' @export
eval_JSD_plot <- function(res_JSD, top_var = NA) {
  if (top_var > nrow(res_JSD) || is.na(top_var) ) {
    top_var <- nrow(res_JSD)
    warning("Number of 'top_var' is ", top_var, " and the number of variables is n = ", nrow(res_JSD), ".\n",
            "top_var is set to n = ", nrow(res_JSD)) 
  }
  res_JSD <- res_JSD |> arrange(desc(JSD))
  res_JSD <- res_JSD[1:top_var,]
  
  ggplot(res_JSD, aes(x = reorder(variable, JSD), y = JSD, fill = type)) +
    geom_bar(stat = "identity") +
    coord_flip() +  
    labs(title = "Jensen-Shannon Divergence per Variable",
         x = "Variable",
         y = "JSD (bits)") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
}

# two sample) Kolmogorov-Smirnov test

# Internal: two sample Kolmogorov-Smirnov test statistic
KS_test_stat <- function(real_vec, syn_vec, power = 1){
  n_real <- length(real_vec)
  n_syn <- length(syn_vec)
  n  <- n_real + n_syn
  
  joint_vector <- c(real_vec, syn_vec)
  # increments for each emperical cummulative distribution function, separately for real and synthetic data.
  r <- c(rep(1/n_real, n_real), rep(0, n_syn))
  s <- c(rep(0, n_real), rep(1/n_syn, n_syn))
  
  vector_order <- order(joint_vector, decreasing = F) # Order from smallest to largest
  d   <- joint_vector[vector_order]
  r   <- r[vector_order]
  s   <- s[vector_order]
  
  KS_statistic <- 0
  Rcur   <- 0 
  Scur   <- 0 
  height <- 0 
  
  for (i in 1:(n - 1)) {
    Rcur   <- Rcur + r[i]
    Scur   <- Scur + s[i]
    if (d[i] != d[i + 1]) {
      height = abs(Scur - Rcur)
    }
    if (height > KS_statistic) {
      KS_statistic <- height
    }  
  }
  KS_statistic^power
}

#' Kolmogorov-Smirnov test per variable
#'
#' Computes the two sample Kolmogorov-Smirnov statistic for every numeric
#' variable. The p-value is obtained by randomly reshuffling the real and
#' synthetic values `nboots` times. Binary and categorical variables are
#' skipped.
#'
#' @param real_df,syn_df data.frames with the real and the synthetic data
#'   (same column names).
#' @param nboots Number of random reshuffles used for the p-value.
#' @param max_unique_for_cat Numeric columns with at most this many unique
#'   values are treated as categorical (see [detect_types()]).
#' @return A data.frame with one row per variable and the columns
#'   `variable`, `type`, `statistic`, `p_value`, `n_real`, `n_syn`, `nboots`
#'   and `note`.
#' @examples
#' \dontrun{
#' res_KS <- eval_KS(real_df, syn_df)
#' eval_KS_plot(res_KS)
#' }
#' @export
eval_KS <- function(real_df, syn_df, nboots = 2000, max_unique_for_cat = 5) {
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
  types_real_df <- detect_types(real_df, max_unique_for_cat)  
  types_syn_df <- detect_types(syn_df, max_unique_for_cat)  
  
  if ("categorical" %in% types_real_df$type || "binary" %in% types_real_df$type) {
    warning(paste(c(paste0("One or multiple variables are categorical and/or binary and will be omitted from the analysis: ", 
                           paste(toString(types_real_df[types_syn_df$type == "categorical" | types_syn_df$type == "binary", "variable"]))))))
  }
  
  results_list <- lapply(vars, function(var) {
    var_real <- real_df[[var]] 
    var_syn <- syn_df[[var]]
    
    type <- types_real_df[types_real_df$variable == var,]$type
    
    if (!(types_real_df[types_real_df$variable == var,]$type == types_syn_df[types_syn_df$variable == var,]$type)) {
      warning(paste(c(paste0("Error: real_df and syn_df have different type for variable: ", paste(var))),
                    c(paste0("real_df type is ", paste(types_real_df[types_real_df$variable == var,]$type))),
                    c(paste0("syn_df type is ", paste(types_syn_df[types_syn_df$variable == var,]$type))), sep = "\n"))
    }
    
    if (type != "numeric") {
      # Non-continuous: mark as missing
      return(data.frame(
        variable  = var,
        type      = type,
        statistic = NA_real_,
        p_value   = NA_real_,
        n_real    = sum(!is.na(var_real)),
        n_syn     = sum(!is.na(var_syn)),
        nboots    = nboots,
        note      = "Non-continuous (binary/categorical) - KS not computed",
        stringsAsFactors = FALSE))
    }
    
    var_real <- sort(as.numeric(var_real), decreasing = FALSE, na.last = NA) #NAs are removed
    var_syn <- sort(as.numeric(var_syn), decreasing = FALSE, na.last = NA)  #NAs are removed
    n_real <- length(var_real)
    n_syn <- length(var_syn)
    
    # Need at least 2 values in each sample
    if (n_real < 2 || n_syn < 2) {
      return(data.frame(
        variable  = var,
        type      = type,
        statistic = NA_real_,
        p_value   = NA_real_,
        n_real    = n_real,
        n_syn     = n_syn,
        nboots    = nboots,
        note      = "Not enough non-missing values for KS",
        stringsAsFactors = FALSE))
    }
    
    # Observed statistic (on sorted vectors)
    KS <- KS_test_stat(var_real, var_syn)
    
    # Bootstrap for p-value
    comb <- c(var_real, var_syn)
    count_larger <- 0L
    
    for (b in seq_len(nboots)) {
      idx1      <- sample.int(n = length(comb), size = n_real, replace = FALSE) # Indices to reshuffle real and synthetic data.
      a_b       <- comb[idx1]
      b_b       <- comb[-idx1]
      boot_KS <- KS_test_stat(a_b, b_b)
      if (boot_KS >= KS) count_larger <- count_larger + 1L
    }
    
    p_val <- (count_larger + 1) / (nboots + 1)
    if (p_val == 0) p_val <- 1 / (2 * nboots)  
    
    data.frame(
      variable  = var,
      type      = type,
      statistic = KS,
      p_value   = p_val,
      n_real    = n_real,
      n_syn     = n_syn,
      nboots    = nboots,
      note      = NA_character_,
      stringsAsFactors = FALSE)
  })
  
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL
  results
  
}

#' Plot the Kolmogorov-Smirnov results
#'
#' @param res_KS Output of [eval_KS()].
#' @param top_var Number of variables (lowest p-values first) that get a
#'   label (`method = "map"`) or are shown (`method = "bar"`). `NA` uses all
#'   variables.
#' @param method `"map"` (statistic against -log10 p-value) or `"bar"`.
#' @return A ggplot object.
#' @export
eval_KS_plot <- function(res_KS, top_var = 20, method = "map")  {
  if (top_var > nrow(res_KS) || is.na(top_var) ) {
    top_var <- nrow(res_KS)
    warning("Number of 'top_var' is ", top_var, " and the number of variables is n = ", nrow(res_KS), ".\n",
            "top_var is set to n = ", nrow(res_KS)) 
  }
  
  res_plot <- subset(res_KS, !is.na(statistic))
  res_plot$neglog10_p <- -log10(res_plot$p_value)
  res_plot <- res_plot %>% arrange(p_value) %>% mutate(label_to_show = if_else(row_number() <= ifelse(is.na(top_var), nrow(res_plot), top_var), variable, ""))
  
  if (method == "map") {
    ggplot(res_plot, aes(x = statistic, y = neglog10_p, color = p_value, label = variable)) +
      geom_point(size = 3) +
      ggrepel::geom_text_repel(aes(label = label_to_show), size = 3, max.overlaps = ifelse(is.na(top_var), 20, Inf)) +
      scale_color_gradientn(
        colours = c("darkblue", "purple3", "red"),
        values  = scales::rescale(c(1, 0.1, 0)), 
        limits = c(0,1),
        name    = "p-value") +
      labs(
        x = "Kolmogorov-Smirnov test",
        y = expression(-log[10](p)),
        title = "Real vs synthetic similarity (per variable)",
        subtitle = "Lower KS and higher p-values indicate better similarity" ) +
      theme_minimal()
  } else if (method == "bar") {
    res_top <- res_plot %>% arrange(p_value) %>% slice(1:ifelse(is.na(top_var), nrow(res_plot), top_var))
    res_top <- res_top %>% mutate(p_label = if_else(p_value <= 0.05, sprintf("p=%.3f", p_value), ""))
    
    ggplot(res_top, aes(x = reorder(variable, statistic), y = statistic, fill = p_value)) +
      geom_col() +
      geom_text(aes(label = p_label), hjust = -0.1, size = 3) +
      coord_flip(ylim = c(0, max(res_top$statistic) * 1.15)) +
      scale_fill_gradientn(
        colours = c("darkblue", "purple3", "red"),
        values  = scales::rescale(c(1, 0.1, 0)), 
        limits = c(0,1),
        name    = "p-value"
      ) +
      labs(
        x = "Variable",
        y = "Kolmogorov-Smirnov test",
        title = "Per-variable KS statistics: real vs synthetic",
        subtitle = "Bar height = KS; color = p-value"
      ) +
      theme_minimal()
  } else { stop("Method does not exist. Choose 'map' or 'bar'") }
  
}

# Cramer-von Mises (two sample) test

# Internal: two sample Cramer-von Mises test statistic
CvM_test_stat <- function(real_vec, syn_vec, power = 2) {
  n_real <- length(real_vec)
  n_syn <- length(syn_vec)
  n  <- n_real + n_syn
  
  joint_vector <- c(real_vec, syn_vec)
  
  # increments for each emperical cummulative distribution function, separately for real and synthetic data.
  r <- c(rep(1/n_real, n_real), rep(0, n_syn))
  s <- c(rep(0, n_real), rep(1/n_syn, n_syn))
  
  vector_order <- order(joint_vector, decreasing = F) # Order from smallest to largest
  d   <- joint_vector[vector_order]
  r   <- r[vector_order]
  s   <- s[vector_order]
  
  CvM_statistic <- 0
  Rcur   <- 0 # current value of empirical CDF of real sample
  Scur   <- 0 # current value of empirical CDF of synthetic sample
  height <- 0 # absolute difference (|Rcur - Scur|0)
  duplicates   <- 1 # count of duplicates of a value
  
  for (i in 1:(n - 1)) {
    Rcur   <- Rcur + r[i]
    Scur   <- Scur + s[i]
    height <- abs(Rcur - Scur)
    
    if (d[i] != d[i + 1]) {
      CvM_statistic  <- CvM_statistic + (height^power) * duplicates
      duplicates <- 1
    } else {  duplicates <- duplicates + 1  }
  }
  CvM_statistic
}

#' Cramer-von Mises test per variable
#'
#' Same as [eval_KS()], but with the two sample Cramer-von Mises statistic.
#'
#' @inheritParams eval_KS
#' @inherit eval_KS return
#' @examples
#' \dontrun{
#' res_CvM <- eval_CvM(real_df, syn_df)
#' eval_CvM_plot(res_CvM)
#' }
#' @export
eval_CvM <- function(real_df, syn_df, nboots = 2000, max_unique_for_cat = 5) {
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
  types_real_df <- detect_types(real_df, max_unique_for_cat)  
  types_syn_df <- detect_types(syn_df, max_unique_for_cat)  
  
  if ("categorical" %in% types_real_df$type || "binary" %in% types_real_df$type) {
    warning(paste(c(paste0("One or multiple variables are categorical and/or binary and will be omitted from the analysis: ", 
                           paste(toString(types_real_df[types_syn_df$type == "categorical" | types_syn_df$type == "binary", "variable"]))))))
  }
  
  results_list <- lapply(vars, function(var) {
    var_real <- real_df[[var]] 
    var_syn <- syn_df[[var]]
    
    type <- types_real_df[types_real_df$variable == var,]$type
    
    if (!(types_real_df[types_real_df$variable == var,]$type == types_syn_df[types_syn_df$variable == var,]$type)) {
      warning(paste(c(paste0("Error: real_df and syn_df have different type for variable: ", paste(var))),
                    c(paste0("real_df type is ", paste(types_real_df[types_real_df$variable == var,]$type))),
                    c(paste0("syn_df type is ", paste(types_syn_df[types_syn_df$variable == var,]$type))), sep = "\n"))
    }
    
    if (type != "numeric") {
      # Non-continuous: mark as missing
      return(data.frame(
        variable  = var,
        type      = type,
        statistic = NA_real_,
        p_value   = NA_real_,
        n_real    = sum(!is.na(var_real)), # number of non-missing
        n_syn     = sum(!is.na(var_syn)),
        nboots    = nboots,
        note      = "Non-continuous (binary/categorical) - CvM not computed",
        stringsAsFactors = FALSE))
    }
    
    # Continuous numeric: compute CvM
    var_real <- sort(as.numeric(var_real), decreasing = FALSE, na.last = NA) #NAs are removed
    var_syn <- sort(as.numeric(var_syn), decreasing = FALSE, na.last = NA)  #NAs are removed
    n_real <- length(var_real)
    n_syn <- length(var_syn)
    
    # Need at least 2 values in each sample
    if (n_real < 2 || n_syn < 2) {
      return(data.frame(
        variable  = var,
        type      = type,
        statistic = NA_real_,
        p_value   = NA_real_,
        n_real    = n_real,
        n_syn     = n_syn,
        nboots    = nboots,
        note      = "Not enough non-missing values for CvM",
        stringsAsFactors = FALSE))
    }
    
    # Observed statistic (on sorted vectors)
    CvM <- CvM_test_stat(var_real, var_syn, power = 2)
    
    # Bootstrap for p-value
    comb <- c(var_real, var_syn)
    count_larger <- 0L
    
    for (b in seq_len(nboots)) {
      idx1      <- sample.int(n = length(comb), size = n_real, replace = FALSE) # Indices to reshuffle real and synthetic data.
      a_b       <- comb[idx1]
      b_b       <- comb[-idx1]
      boot_CvM <- CvM_test_stat(a_b, b_b, power = 2)
      if (boot_CvM >= CvM) count_larger <- count_larger + 1L
    }
    
    p_val <- (count_larger + 1) / (nboots + 1)
    if (p_val == 0) p_val <- 1 / (2 * nboots)  
    
    data.frame(
      variable  = var,
      type      = type,
      statistic = CvM,
      p_value   = p_val,
      n_real    = n_real,
      n_syn     = n_syn,
      nboots    = nboots,
      note      = NA_character_,
      stringsAsFactors = FALSE)
  })
  
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL
  results
}

#' Plot the Cramer-von Mises results
#'
#' @param res_CvM Output of [eval_CvM()].
#' @inheritParams eval_KS_plot
#' @return A ggplot object.
#' @export
eval_CvM_plot <- function(res_CvM, top_var = NA, method = "map")  {
  if (top_var > nrow(res_CvM) || is.na(top_var) ) {
    top_var <- nrow(res_CvM)
    warning("Number of 'top_var' is ", top_var, " and the number of variables is n = ", nrow(res_CvM), ".\n",
            "top_var is set to n = ", nrow(res_CvM)) 
  }
  
  res_plot <- subset(res_CvM, !is.na(statistic))
  res_plot$neglog10_p <- -log10(res_plot$p_value)
  res_plot <- res_plot %>% arrange(p_value) %>% mutate(label_to_show = if_else(row_number() <= ifelse(is.na(top_var), nrow(res_plot), top_var), variable, ""))
  
  if (method == "map") {
    ggplot(res_plot, aes(x = statistic, y = neglog10_p, color = p_value, label = variable)) +
      geom_point(size = 3) +
      ggrepel::geom_text_repel(aes(label = label_to_show), size = 3, max.overlaps = ifelse(is.na(top_var), 20, Inf)) +
      scale_color_gradientn(
        colours = c("darkblue", "purple3", "red"),
        values  = scales::rescale(c(1, 0.1, 0)), 
        name    = "p-value") +
      labs(
        x = "Cramer-von Mises statistic",
        y = expression(-log[10](p)),
        title = "Real vs synthetic similarity (per variable)",
        subtitle = "Lower CvM and higher p-values indicate better similarity" ) +
      theme_minimal()
  } else if (method == "bar") {
    res_top <- res_plot %>% arrange(p_value) %>% slice(1:ifelse(is.na(top_var), nrow(res_plot), top_var))
    res_top <- res_top %>% mutate(p_label = if_else(p_value <= 0.05, sprintf("p=%.3f", p_value), ""))
    
    ggplot(res_top, aes(x = reorder(variable, statistic), y = statistic, fill = p_value)) +
      geom_col() +
      geom_text(aes(label = p_label), hjust = -0.1, size = 3) +
      coord_flip(ylim = c(0, max(res_top$statistic) * 1.15)) +
      scale_fill_gradientn(
        colours = c("darkblue", "purple3", "red"),
        values  = scales::rescale(c(1, 0.1, 0)), 
        limits = c(0,1),
        name    = "p-value"
      ) +
      labs(
        x = "Variable",
        y = "Cramer-von Mises statistic",
        title = "Per-variable CvM statistics: real vs synthetic",
        subtitle = "Bar height = CvM; color = p-value"
      ) +
      theme_minimal()
  } else { stop("Method does not exist. Choose 'map' or 'bar'") }
}

# (two sample) Anderson-Darling test

# Internal: two sample Anderson-Darling test statistic
AD_test_stat <- function(real_vec, syn_vec, power = 2){
  n_real <- length(real_vec)
  n_syn <- length(syn_vec)
  n  <- n_real + n_syn
  
  joint_vector <- c(real_vec, syn_vec)
  # increments for each emperical cummulative distribution function, separately for real and synthetic data.
  r <- c(rep(1/n_real, n_real), rep(0, n_syn))
  s <- c(rep(0, n_real), rep(1/n_syn, n_syn))
  
  vector_order <- order(joint_vector, decreasing = F) # Order from smallest to largest
  d   <- joint_vector[vector_order]
  r   <- r[vector_order]
  s   <- s[vector_order]
  
  AD_statistic <- 0
  Rcur   <- 0 # current value of empirical CDF of real sample
  Scur   <- 0 # current value of empirical CDF of synthetic sample
  Jcur   <- 0 # ECDF of joint sample
  height <- 0 # absolute difference (|Rcur - Scur|0)
  duplicates <- 1
  
  for (i in 1:(n - 1)) {
    Rcur   <- Rcur + r[i]
    Scur   <- Scur + s[i]
    Jcur   <- Jcur + 1/n
    sd <- (2*Jcur*(1-Jcur)/n)^0.5
    height <- abs(Scur - Rcur)
    if (d[i] != d[i + 1]) {
      AD_statistic <- AD_statistic + ((height / sd) ^ power) * duplicates
      duplicates <- 1
    } else if (d[i] == d[i + 1]) {
      duplicates = duplicates + 1
    }  
  }
  AD_statistic
}

#' Anderson-Darling test per variable
#'
#' Same as [eval_KS()], but with the two sample Anderson-Darling statistic.
#'
#' @inheritParams eval_KS
#' @inherit eval_KS return
#' @examples
#' \dontrun{
#' res_AD <- eval_AD(real_df, syn_df)
#' eval_AD_plot(res_AD)
#' }
#' @export
eval_AD <- function(real_df, syn_df, nboots = 2000, max_unique_for_cat = 5) {
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
  types_real_df <- detect_types(real_df, max_unique_for_cat)  
  types_syn_df <- detect_types(syn_df, max_unique_for_cat)  
  
  if ("categorical" %in% types_real_df$type || "binary" %in% types_real_df$type) {
    warning(paste(c(paste0("One or multiple variables are categorical and/or binary and will be omitted from the analysis: ", 
                           paste(toString(types_real_df[types_syn_df$type == "categorical" | types_syn_df$type == "binary", "variable"]))))))
  }
  
  results_list <- lapply(vars, function(var) {
    var_real <- real_df[[var]] 
    var_syn <- syn_df[[var]]
    
    type <- types_real_df[types_real_df$variable == var,]$type
    
    if (!(types_real_df[types_real_df$variable == var,]$type == types_syn_df[types_syn_df$variable == var,]$type)) {
      warning(paste(c(paste0("Error: real_df and syn_df have different type for variable: ", paste(var))),
                    c(paste0("real_df type is ", paste(types_real_df[types_real_df$variable == var,]$type))),
                    c(paste0("syn_df type is ", paste(types_syn_df[types_syn_df$variable == var,]$type))), sep = "\n"))
    }
    
    if (type != "numeric") {
      # Non-continuous: mark as missing
      return(data.frame(
        variable  = var,
        type      = type,
        statistic = NA_real_,
        p_value   = NA_real_,
        n_real    = sum(!is.na(var_real)),
        n_syn     = sum(!is.na(var_syn)),
        nboots    = nboots,
        note      = "Non-continuous (binary/categorical) - AD not computed",
        stringsAsFactors = FALSE))
    }
    
    var_real <- sort(as.numeric(var_real), decreasing = FALSE, na.last = NA) #NAs are removed
    var_syn <- sort(as.numeric(var_syn), decreasing = FALSE, na.last = NA)  #NAs are removed
    n_real <- length(var_real)
    n_syn <- length(var_syn)
    
    # Need at least 2 values in each sample
    if (n_real < 2 || n_syn < 2) {
      return(data.frame(
        variable  = var,
        type      = type,
        statistic = NA_real_,
        p_value   = NA_real_,
        n_real    = n_real,
        n_syn     = n_syn,
        nboots    = nboots,
        note      = "Not enough non-missing values for AD",
        stringsAsFactors = FALSE))
    }
    
    # Observed statistic (on sorted vectors)
    AD <- AD_test_stat(var_real, var_syn)
    
    # Bootstrap for p-value
    comb <- c(var_real, var_syn)
    count_larger <- 0L
    
    for (b in seq_len(nboots)) {
      idx1      <- sample.int(n = length(comb), size = n_real, replace = FALSE) # Indices to reshuffle real and synthetic data.
      a_b       <- comb[idx1]
      b_b       <- comb[-idx1]
      boot_AD <- AD_test_stat(a_b, b_b)
      if (boot_AD >= AD) count_larger <- count_larger + 1L
    }
    
    p_val <- (count_larger + 1) / (nboots + 1)
    if (p_val == 0) p_val <- 1 / (2 * nboots)  # small correction
    
    data.frame(
      variable  = var,
      type      = type,
      statistic = AD,
      p_value   = p_val,
      n_real    = n_real,
      n_syn     = n_syn,
      nboots    = nboots,
      note      = NA_character_,
      stringsAsFactors = FALSE)
  })
  
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL
  results
  
}

#' Plot the Anderson-Darling results
#'
#' @param res_AD Output of [eval_AD()].
#' @inheritParams eval_KS_plot
#' @return A ggplot object.
#' @export
eval_AD_plot <- function(res_AD, top_var = 20, method = "map")  {
  if (top_var > nrow(res_AD) || is.na(top_var)) {
    top_var <- nrow(res_AD)
    warning("Number of 'top_var' is ", top_var, " and the number of variables is n = ", nrow(res_AD), ".\n",
            "top_var is set to n = ", nrow(res_AD)) 
  }
  
  res_plot <- subset(res_AD, !is.na(statistic))
  res_plot$neglog10_p <- -log10(res_plot$p_value)
  res_plot <- res_plot %>% arrange(p_value) %>% mutate(label_to_show = if_else(row_number() <= ifelse(is.na(top_var), nrow(res_plot), top_var), variable, ""))
  
  if (method == "map") {
    ggplot(res_plot, aes(x = statistic, y = neglog10_p, color = p_value, label = variable)) +
      geom_point(size = 3) +
      ggrepel::geom_text_repel(aes(label = label_to_show), size = 3, max.overlaps = ifelse(is.na(top_var), 20, Inf)) +
      scale_color_gradientn(
        colours = c("darkblue", "purple3", "red"),
        values  = scales::rescale(c(1, 0.1, 0)), 
        name    = "p-value") +
      labs(
        x = "Anderson-Darling",
        y = expression(-log[10](p)),
        title = "Real vs synthetic similarity (per variable)",
        subtitle = "Lower AD and higher p-values indicate better similarity" ) +
      theme_minimal()
  } else if (method == "bar") {
    res_top <- res_plot %>% arrange(p_value) %>% slice(1:ifelse(is.na(top_var), nrow(res_plot), top_var))
    res_top <- res_top %>% mutate(p_label = if_else(p_value <= 0.05, sprintf("p=%.3f", p_value), ""))
    
    ggplot(res_top, aes(x = reorder(variable, statistic), y = statistic, fill = p_value)) +
      geom_col() +
      geom_text(aes(label = p_label), hjust = -0.1, size = 3) +
      coord_flip(ylim = c(0, max(res_top$statistic) * 1.15)) +
      scale_fill_gradientn(
        colours = c("darkblue", "purple3", "red"),
        values  = scales::rescale(c(1, 0.1, 0)), 
        name    = "p-value"
      ) +
      labs(
        x = "Variable",
        y = "Anderson-Darling",
        title = "Per-variable AD statistics: real vs synthetic",
        subtitle = "Bar height = AD; color = p-value"
      ) +
      theme_minimal()
  } else { stop("Method does not exist. Choose 'map' or 'bar'") }
  
}

# pMSE

#' Propensity score mean squared error (pMSE) per variable
#'
#' Wrapper around [synthpop::compare()].
#'
#' @param real_df,syn_df data.frames with the real and the synthetic data
#'   (same column names).
#' @param stat,plot,utility.stats,... Passed to [synthpop::compare()].
#' @param max_unique_for_cat Numeric columns with at most this many unique
#'   values are treated as categorical (see [detect_types()]).
#' @return A list with `comparison` (the output of [synthpop::compare()])
#'   and `result_table` (the requested utility statistics per variable).
#' @examples
#' \dontrun{
#' res_pMSE <- eval_pMse(real_df, syn_df)
#' }
#' @export
eval_pMse <- function(real_df, syn_df, stat = "percents", plot = TRUE, utility.stats = c("pMSE", "S_pMSE"), max_unique_for_cat = 5, ...) {
  
  if (!is.data.frame(real_df) || !is.data.frame(syn_df)) {
    stop("Both real_df and syn_df must be data.frames.")
  }
  real_df <- data.frame(real_df)
  syn_df  <- data.frame(syn_df)
  if (!setequal(names(real_df), names(syn_df))) {
    missing_in_syn  <- setdiff(names(real_df), names(syn_df))
    missing_in_real <- setdiff(names(syn_df), names(real_df))
    stop(paste(c("real_df and syn_df must have the same variable names. ",
                 c(if (length(missing_in_syn))  paste0("Missing in syn_df: ",  paste(missing_in_syn,  collapse=", ")),
                   if (length(missing_in_real)) paste0("Missing in real_df: ", paste(missing_in_real, collapse=", ")))), 
               sep="\n"))
  }
  
  vars <- names(real_df)
  types_real_df <- detect_types(real_df, max_unique_for_cat)  
  types_syn_df <- detect_types(syn_df, max_unique_for_cat)  
  
  if(any(!types_real_df[match(types_real_df$variable, types_real_df$variable),]$type == types_syn_df[match(types_syn_df$variable, types_real_df$variable),]$type))  {
    warning(paste(c(paste0("real_df and syn_df have different type for variables: ", 
                           toString(types_real_df$variable[!types_real_df[match(types_real_df$variable, types_real_df$variable),]$type == types_syn_df[match(types_syn_df$variable, types_real_df$variable),]$type])))))
  }
  
  comp   <- synthpop::compare(syn_df, real_df, stat = stat, plot = plot, utility.stats = utility.stats, ...)
  result <- as.data.frame(comp$tab.utility) %>%  rownames_to_column("Variable") %>% select(Variable, paste0(utility.stats))
  return(list(comparison = comp, result_table = result))
}
