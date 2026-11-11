# Privacy evaluators: distance to closest record (DCR) and nearest neighbor
# distance ratio (NNDR). The shared distance helpers are in utils.R.

# DCR

#' Distance to closest record (DCR)
#'
#' For every synthetic row, the distance to the closest real row. As a
#' baseline, the distance of every real row to the closest other real row is
#' also computed.
#'
#' @param real_df,syn_df data.frames with the real and the synthetic data
#'   (same column names).
#' @param seed Random seed.
#' @param scale_numeric Standardise the numeric columns (with the mean and sd
#'   of the real data).
#' @param encoding `"onehot"` (Euclidean distance on one-hot encoded data) or
#'   `"mixed"` (Gower distance).
#' @param digits Number of digits in the summary.
#' @param return_value `"summary"`, `"vector"` (the distances of the synthetic
#'   rows) or `"all"`.
#' @param max_unique_for_cat Numeric columns with at most this many unique
#'   values are treated as categorical (see [detect_types()]).
#' @param print_outcome Print the summary.
#' @return Depending on `return_value`: a summary data.frame, a numeric
#'   vector, or a list with both, the baseline and the exact match rate.
#' @examples
#' \dontrun{
#' DCR(real_df, syn_df)
#' }
#' @export
DCR <- function(real_df, syn_df,
                seed = 1,
                scale_numeric = TRUE,
                encoding = c("onehot", "mixed"),
                digits = 3,
                return_value = c("summary", "vector", "all"), max_unique_for_cat = 2, print_outcome = TRUE)
{
  return_value <- match.arg(return_value)
  encoding     <- match.arg(encoding)
  
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' is required. Please install it with install.packages('FNN').")
  }
  if (encoding == "mixed"  && !requireNamespace("gower", quietly = TRUE))
    stop("Package 'gower' is required (install.packages('gower')).")
  
  set.seed(seed)
  
  prep <- .prepare_distance_data(real_df, syn_df, scale_numeric = scale_numeric,
                                 encoding = encoding, max_unique_for_cat = max_unique_for_cat)
  
  # synthetic → nearest real
  dcr_syn  <- as.numeric(.knn_query(query = prep$X_syn, ref = prep$X_real, k = 1, prep$encoding)$dist[, 1])
  exact_match_rate <- mean(dcr_syn == 0) 
  # real → nearest other real (baseline)
  dcr_real <- as.numeric(.knn_self(data = prep$X_real, k = 1, prep$encoding)$dist[, 1])
  
  out_summary <- .summary_stats(values_list = list(DCR_syn = dcr_syn, DCR_real_baseline = dcr_real), digits = digits)
  
  if(print_outcome) {
    cat("Distance to Closest Record (DCR) with encoding:", encoding, "\n")
    print(out_summary, row.names = FALSE)
  }
  
  switch(return_value,
         summary = out_summary,
         vector  = dcr_syn,
         all     = list(summary           = out_summary,
                        dcr_syn           = dcr_syn,
                        dcr_real_baseline = dcr_real,
                        p05_syn           = unname(stats::quantile(dcr_syn,  0.05)),
                        p05_real          = unname(stats::quantile(dcr_real, 0.05)),
                        exact_match_rate  = exact_match_rate
         ))
}

# NNDR

#' Nearest neighbor distance ratio (NNDR)
#'
#' For every synthetic row, the distance to the closest real row divided by
#' the distance to the second closest real row. As a baseline, the same
#' ratio is computed for every real row against the other real rows.
#'
#' @inheritParams DCR
#' @param k `real_df` must have more than `k` rows. The ratio itself always
#'   uses the two nearest neighbors.
#' @param return_value `"summary"`, `"vector"` (the ratios of the synthetic
#'   rows) or `"all"`.
#' @inherit DCR return
#' @examples
#' \dontrun{
#' NNDR(real_df, syn_df)
#' }
#' @export
NNDR <- function(real_df, syn_df,
                 seed = 1,
                 k = 5,
                 scale_numeric = TRUE,
                 encoding = c("onehot", "mixed"),
                 digits = 3,
                 return_value = c("summary", "vector", "all"), max_unique_for_cat = 2, print_outcome = TRUE)
{
  return_value <- match.arg(return_value)
  encoding     <- match.arg(encoding)
  
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' is required. Please install it with install.packages('FNN').")
  }
  if (encoding == "mixed"  && !requireNamespace("gower", quietly = TRUE))
    stop("Package 'gower' is required (install.packages('gower')).")
  if (!is.numeric(k) || length(k) != 1 || k < 1) {
    stop("k must be a single integer >= 1.")
  }
  
  set.seed(seed)
  
  prep <- .prepare_distance_data(real_df, syn_df,
                                 scale_numeric      = scale_numeric,
                                 encoding           = encoding,
                                 max_unique_for_cat = max_unique_for_cat)
  
  if (nrow(prep$X_real) <= k || nrow(prep$X_real) < 3) {
    stop(paste0("real_df must contain at least ", min(k, 3), " rows to compute NNDR."))
  }
  
  # synthetic → 2 nearest real datapoins
  nn_syn   <- .knn_query(query = prep$X_syn, ref = prep$X_real, k = 2, encoding = prep$encoding)
  d1_syn   <- nn_syn$dist[, 1]; d2_syn <- nn_syn$dist[, 2]
  nndr_syn <- ifelse(d2_syn  == 0, 1, d1_syn  / d2_syn)
  exact_match_rate <- mean(d1_syn == 0) # Number of synthetic rows that are exact matches with the real dataset.
  
  # real → 2 nearest OTHER real (baseline)
  nn_real   <- .knn_self(data = prep$X_real, k = 2, encoding = prep$encoding)
  d1_real   <- nn_real$dist[, 1]; d2_real <- nn_real$dist[, 2]
  nndr_real <- ifelse(d2_real == 0, 1, d1_real / d2_real)
  
  out_summary <- .summary_stats(
    list(NNDR_syn = nndr_syn, NNDR_real_baseline = nndr_real),
    digits = digits
  )
  
  if(print_outcome) {
    cat("Nearest Neighbor Distance Ratio (NNDR) with encoding:", encoding, "\n")
    print(out_summary, row.names = FALSE)
  }
  
  switch(return_value,
         summary = out_summary,
         vector  = nndr_syn,
         all     = list(summary = out_summary,
                        nndr_syn = nndr_syn,
                        nndr_real_baseline = nndr_real,
                        p05_syn = unname(stats::quantile(nndr_syn,  0.05)),
                        p05_real = unname(stats::quantile(nndr_real, 0.05)),
                        exact_match_rate = exact_match_rate
         ))
  
}
