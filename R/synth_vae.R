#' Generate synthetic data with a tabular VAE (TVAE)
#'
#' Fits a Tabular Variational Autoencoder using the Python package `ctgan`
#' (through reticulate), optionally with a grid search over hyperparameters.
#'
#' @param real_df A data.frame with the real data.
#' @param discrete_columns Character vector with the names of the discrete (categorical) columns.
#' @param tune If `TRUE`, it runs a grid search and returns the best model.
#' @param tune_grid A data.frame of settings for the grid search to try, one row per run (see 
#' Details). If `NULL`, a default grid is used.
#' @param return_tune If `TRUE` (and `tune = TRUE`), also return the full
#'   tuning results.
#' @param tune_eval Metric used to pick the best run: `"MWS"` (Wasserstein
#'   distance, see [eval_MWS()]) or `"quality"` (mean KS / TVD). Lower is better.
#' @param seed Optional random seed for R, NumPy, Python and torch.
#' @param embedding_dim,compress_dims,decompress_dims,l2scale,batch_size,epochs,loss_factor,enable_gpu,verbose
#'   TVAE settings for a single run (ignored when `tune = TRUE`). See the
#'   `ctgan` Python documentation.
#' @param return_model If `TRUE`, also return the fitted Python ctgan model. A list with the
#'   synthetic data (`synth`) and the fitted Python model (`model`) is returned.
#'  
#' @details `tune_grid` needs the columns `embedding_dim`, `compress_dims`, `decompress_dims`, `l2scale`,
#'   `batch_size`, `epochs`, `loss_factor`, `enable_gpu` and `verbose`. Use
#'   `I(list(...))` for `compress_dims` and `decompress_dims`, for example:
#'   \preformatted{
#'   grid <- expand.grid(embedding_dim = c(64L, 128L, 265L), compress_dims = list(c(256L, 256L), c(128L, 128L), c(64L, 64L)), 
#'     decompress_dims = list(c(256L, 256L), c(128L, 128L), c(64L, 64L)), l2scale = c(1e-4,1e-5,1e-6),
#'     batch_size = c(500L), epochs = c(300L), loss_factor = 2L, enable_gpu = TRUE, verbose = FALSE,  KEEP.OUT.ATTRS = FALSE )
#'   }
#'   
#' @return A data.frame, or a list with `syn`, and (depending on the
#'   arguments) `model`, `param` and `tune`.
#' @examples
#' \dontrun{
#' syn_df <- synth_vae(penguins, discrete_columns = c("species", "island", "sex"))
#' }
#' @export
synth_vae <- function(real_df,
                      discrete_columns,
                      tune = FALSE,
                      tune_grid = NULL,
                      return_tune = TRUE,
                      tune_eval = "MWS",
                      return_model = FALSE,
                      seed = NULL,
                      embedding_dim = 128L,
                      compress_dims = c(128L, 128L),
                      decompress_dims = c(128L, 128L),
                      l2scale = 1e-5,
                      batch_size = 500L,
                      epochs = 300L,
                      loss_factor = 2L,
                      enable_gpu = TRUE,
                      verbose = FALSE) {
  stopifnot(is.data.frame(real_df))
  stopifnot(is.character(discrete_columns))
  stopifnot(all(discrete_columns %in% names(real_df)))
  
  ctgan      <- .import_ctgan()
  py_builtin <- reticulate::import("builtins")
  .set_seed_all(seed)
  
  py_df         <- reticulate::r_to_py(real_df, convert = FALSE)
  discrete_cols <- as.list(discrete_columns)
  numeric_cols  <- setdiff(names(real_df), discrete_columns)
  
  # fit and sample one configuration
  .fit_and_sample <- function(embedding_dim, compress_dims, decompress_dims, l2scale,
                              batch_size, epochs, loss_factor, enable_gpu, verbose) {
    model <- ctgan$TVAE(
      embedding_dim   = as.integer(embedding_dim),
      compress_dims   = py_builtin$tuple(as.integer(compress_dims)),
      decompress_dims = py_builtin$tuple(as.integer(decompress_dims)),
      l2scale         = l2scale,
      batch_size      = as.integer(batch_size),
      epochs          = as.integer(epochs),
      loss_factor     = loss_factor,
      enable_gpu      = enable_gpu,
      verbose         = verbose
    )
    model$fit(py_df, discrete_cols)
    syn_df <- .as_r_df(model$sample(as.integer(nrow(real_df))))
    list(synth = syn_df, model = model)
  }
  
  # single run
  if (!tune) {
    res <- .fit_and_sample(embedding_dim, compress_dims, decompress_dims, l2scale,
                           batch_size, epochs, loss_factor, enable_gpu, verbose)
    if (return_model) return(list(syn = res$synth, model = res$model))
    return(res$synth)
  }
  
  # tuning
  if (is.null(tune_grid)) {
    tune_grid <- expand.grid(
      embedding_dim   = c(64L, 128L, 256L),
      compress_dims   = I(list(c(256L, 256L), c(128L, 128L), c(64L, 64L))),
      decompress_dims = I(list(c(256L, 256L), c(128L, 128L), c(64L, 64L))),
      l2scale         = c(1e-4, 1e-5, 1e-6),
      batch_size      = 500L,
      epochs          = 300L,
      loss_factor     = 2L,
      enable_gpu      = TRUE,
      verbose         = FALSE,
      KEEP.OUT.ATTRS  = FALSE)
  }
  
  fit_row <- function(p) {
    .fit_and_sample(embedding_dim   = as.integer(p$embedding_dim),
                    compress_dims   = p$compress_dims[[1]],
                    decompress_dims = p$decompress_dims[[1]],
                    l2scale         = p$l2scale,
                    batch_size      = as.integer(p$batch_size),
                    epochs          = as.integer(p$epochs),
                    loss_factor     = p$loss_factor,
                    enable_gpu      = p$enable_gpu,
                    verbose         = p$verbose)
  }
  
  best <- .tune_loop(tune_grid, fit_row, real_df, numeric_cols,
                     discrete_columns, tune_eval)
  .tune_output(best, return_model, return_tune)
}
