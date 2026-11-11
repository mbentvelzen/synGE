# Bivariate evaluation: difference in pairwise correlations.

#' Compare the pairwise correlations of real and synthetic data
#'
#' Computes the correlation matrix of the real and of the synthetic data and
#' the (absolute) difference between the two.
#'
#' @param real_df,syn_df data.frames with the real and the synthetic data
#'   (same column names).
#' @param method,use Passed to [stats::cor()].
#' @param max_unique_for_cat Numeric columns with at most this many unique
#'   values are treated as categorical (see [detect_types()]).
#' @param make_numeric If `TRUE`, categorical and binary variables are
#'   recoded to dummy columns and included. If `FALSE`, non-numeric
#'   variables are omitted.
#' @return A list with `result` (one row per pair of variables, sorted by
#'   `abs_diff`, with a `"summary"` attribute) and `cor_dif` (matrix of
#'   absolute correlation differences).
#' @examples
#' \dontrun{
#' res_Cor <- eval_Cor(real_df, syn_df)
#' eval_Cor_plot(res_Cor)
#' }
#' @export
eval_Cor <- function(real_df, syn_df, method = "pearson", use = "pairwise.complete.obs", max_unique_for_cat = 5, make_numeric = FALSE) {
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
  
  if(any(!types_real_df[match(types_real_df$variable, types_real_df$variable),]$type == types_syn_df[match(types_syn_df$variable, types_real_df$variable),]$type))  {
    warning(paste(c(paste0("real_df and syn_df have different type for variables: ", 
                           toString(types_real_df$variable[!types_real_df[match(types_real_df$variable, types_real_df$variable),]$type == types_syn_df[match(types_syn_df$variable, types_real_df$variable),]$type])))))
  }
  
  if ("categorical" %in% types_real_df$type || "binary" %in% types_real_df$type) {
    
    if (make_numeric) {
      warning(paste(c(paste0("One or multiple variables are categorical and/or binary. These variables are recoded via model.matrix(~., data) to dummy columns to incorporate them in the correlation-difference analysis: ", 
                             paste(toString(types_real_df[types_syn_df$type == "categorical" | types_syn_df$type == "binary", "variable"]))))))
      # Create dummy variables for categorical or binary variables
      for (v in vars) {
        if (is.factor(real_df[[v]]) || is.character(real_df[[v]]) ||
            is.factor(syn_df[[v]])  || is.character(syn_df[[v]])) {
          
          real_chr <- as.character(real_df[[v]])
          syn_chr  <- as.character(syn_df[[v]])
          categories <- sort(unique(c(real_chr[!is.na(real_chr)], syn_chr[!is.na(syn_chr)])))
          real_df[[v]] <- factor(real_chr, levels = categories)
          syn_df[[v]]  <- factor(syn_chr,  levels = categories)
        }
      }
      real_df <- model.matrix(~ 0 + ., real_df)
      syn_df  <- model.matrix(~ 0 + ., syn_df)
    } else {
      warning(paste(c(paste0("One or multiple variables are categorical and/or binary and will be omitted from the analysis: ", 
                             paste(toString(types_real_df[types_syn_df$type == "categorical" | types_syn_df$type == "binary", "variable"]))))))
      real_df <- real_df[vapply(real_df, is.numeric, logical(1))]
      syn_df  <- syn_df[vapply(syn_df, is.numeric, logical(1))]
    }
  }
  
  
  cor_real <- stats::cor(real_df, method = method, use = use)
  cor_syn  <- stats::cor(syn_df,  method = method, use = use)
  cor_dif  <- abs(cor_real - cor_syn)

  index <- which(upper.tri(cor_dif), arr.ind = TRUE)
  result <- data.frame(
    var1     = rownames(cor_dif)[index[, 1]],
    var2     = colnames(cor_dif)[index[, 2]],
    cor_real = cor_real[index],
    cor_syn  = cor_syn[index],
    abs_diff = cor_dif[index],
    method   = method,
    stringsAsFactors = FALSE
  )
  
  # Reorder by highest difference
  result <- result[order(result$abs_diff, decreasing = TRUE), , drop = FALSE]
  
  attr(result, "summary") <- list(
    n_obs_real = nrow(real_df),
    n_obs_syn  = nrow(syn_df),
    n_features = ncol(real_df),
    mean_abs_diff = mean(result$abs_diff, na.rm = TRUE),
    median_abs_diff = stats::median(result$abs_diff, na.rm = TRUE),
    max_abs_diff = max(result$abs_diff, na.rm = TRUE)
  )
  list(result = result, cor_dif = cor_dif)
}

#' Plot the correlation differences
#'
#' @param res_Cor Output of [eval_Cor()].
#' @param top_var Number of variable pairs (largest difference first) to
#'   show.
#' @param plot `"all"`, `"barplot"` or `"heatmap"`.
#' @param max_label_length Longer variable names are truncated.
#' @return A ggplot object, or a list with `barplot` and `heatmap` when both
#'   are requested.
#' @export
eval_Cor_plot <- function(res_Cor, top_var = 20, plot = c("all"), max_label_length = 35) {
  if (top_var > nrow(res_Cor$cor_dif) || is.na(top_var) ) {
    top_var <- nrow(res_Cor$cor_dif)
    warning("Number of 'top_var' is ", top_var, " and the number of variables is n = ", nrow(res_Cor$cor_dif), ".\n",
            "top_var is set to n = ", nrow(res_Cor$cor_dif)) 
  }
  
  plot_list <- list()
  
  truncate_labels <- function(x, max_len = max_label_length) {
    ifelse(
      nchar(x) > max_len,
      paste0(substr(x, 1, max_len), "..."),
      x
    )
  }
  
  if ("all" %in% plot || "barplot" %in% plot) {
    top_res <- res_Cor$result
    res_Cor$result <- res_Cor$result |> arrange(desc(abs_diff))
    top_res <- res_Cor$result[1:top_var, ] 
    
    top_res$pair <- paste(truncate_labels(top_res$var1), truncate_labels(top_res$var2), sep = " \u00D7 ")
    
    plot_list$barplot <- ggplot(top_res, aes(x = reorder(pair, abs_diff), y = abs_diff)) +
      geom_col() +
      coord_flip() +
      labs(
        x = "Feature pair",
        y = "Absolute correlation difference",
        title = "Real vs synthetic correlation differences",
        subtitle = "Top pairs by |cor_real - cor_syn|"
      ) +
      theme_minimal() 
  }
  
  if ("all" %in% plot || "heatmap" %in% plot) {
    cor_matrix <- res_Cor$cor_dif
    diag(cor_matrix) <- 0
    
    res_Cor$result <- res_Cor$result |> arrange(desc(abs_diff))
    top_res <- res_Cor$result[1:top_var, ] 
    keep_vars <- unique(c(top_res$var1, top_res$var2))
    cor_matrix <- cor_matrix[keep_vars, keep_vars, drop = FALSE]
    
    rownames(cor_matrix) <- truncate_labels(rownames(cor_matrix), max_label_length)
    colnames(cor_matrix) <- truncate_labels(colnames(cor_matrix), max_label_length)
    
    row_mean <- rowMeans(cor_matrix, na.rm = TRUE)
    ord <- order(row_mean, decreasing = TRUE)  # or order(row_max[keep_vars], decreasing = TRUE)
    cor_matrix <- cor_matrix[ord, ord, drop = FALSE]
    
    df_long <- as.data.frame(as.table(cor_matrix))
    colnames(df_long) <- c("var1", "var2", "abs_diff")
    
    plot_list$heatmap <- ggplot(df_long, aes(x = var1, y = var2, fill = abs_diff)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_gradientn(
        colours = c("white", "darkblue", "red"),
        limits = c(0,ceiling(max(df_long$abs_diff) / 0.2) * 0.2),
        name = "|\u0394 correlation|"
      ) +
      coord_equal() +
      labs(
        x = NULL,
        y = NULL,
        title = "Absolute correlation difference",
        subtitle = "Real vs synthetic data"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid = element_blank()
      )
  }
  
  # return one plot or both
  if (length(plot_list) == 1) return(plot_list[[1]])
  return(plot_list)
  
}
