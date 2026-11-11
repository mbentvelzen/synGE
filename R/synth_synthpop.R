#' Generate synthetic data with synthpop (CART)
#' 
#' Wrapper around [synthpop::syn()]. Character columns are converted to factors first.
#' @param real_df A data.frame with the real data.
#' @param seed Optional random seed.
#' @param ... Further arguments that are passed to [synthpop::syn()].
#' @return A data.frame with the synthetic data.
#' @examples
#' \dontrun{
#' syn_df <- synth_sp(iris, seed = 123)
#' }
#' @export

synth_sp <- function(real_df, seed = NULL, ...) {
  stopifnot(is.data.frame(real_df))
  
  # 1) mutate characters to factors
  df <- real_df %>%
    mutate(across(where(is.character), as.factor))
  
  # 2) set seed
  if (!is.null(seed)) set.seed(seed)
  
  # 3) create synthetic data
  syn_fit <- synthpop::syn(df, ...)
  as.data.frame(syn_fit$syn)
}
