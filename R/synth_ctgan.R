#' Generate synthetic data with CTGAN
#'
#' Fits a Conditional Tabular GAN using the Python package `ctgan` (through
#' reticulate), optionally with a grid search over hyperparameters.
#'
#' @param real_df A data.frame with the real data.
#' @param discrete_columns Character vector with the names of the discrete (categorical) columns.
#' @param tune If `TRUE`, it runs a grid search and returns the best model.
#' @param tune_grid A data.frame of settings to try for a grid search, one row per run (see
#'   Details). If `NULL`, a default grid is used.
#' @param return_tune If `TRUE` (and `tune = TRUE`), also return the full
#'   tuning results.
#' @param tune_eval Metric used to pick the best run: `"MWS"` (Wasserstein
#'   distance, see [eval_MWS()]) or `"quality"` (mean KS / TVD). Lower is better.
#' @param seed Optional random seed for R, NumPy, Python and torch.
#' @param embedding_dim,generator_dim,discriminator_dim,generator_lr,generator_decay,discriminator_lr,discriminator_decay,batch_size,discriminator_steps,log_frequency,verbose,epochs,pac,enable_gpu
#'   CTGAN settings for a single run (ignored when `tune = TRUE`). See the
#'   `ctgan` Python documentation.
#' @param return_model If `TRUE`, also return the fitted Python ctgan model. A list with the
#'   synthetic data (`synth`) and the fitted Python model (`model`) is returned.
#'
#' @details `tune_grid` needs the columns `epochs`, `batch`, `glr`, `dlr`,
#'   `gdim`, `ddim`, `embed`, `g_decay`, `d_decay`, `d_steps`,
#'   `log_frequency`, `verbose`, `pac` and `enable_gpu`. Use `I(list(...))`
#'   for `gdim` and `ddim`, for example:
#'   \preformatted{
#'   grid <- expand.grid(epochs = c(300L, 500L), batch = 500L, glr = 2e-4,
#'     dlr = 2e-4, gdim = I(list(c(256L, 256L), c(128L, 64L))),
#'     ddim = I(list(c(256L, 256L), c(128L, 64L))), embed = 128L,
#'     g_decay = 1e-6, d_decay = 1e-6, d_steps = 1L, log_frequency = TRUE,
#'     verbose = FALSE, pac = 10L, enable_gpu = TRUE, KEEP.OUT.ATTRS = FALSE)
#'   }
#'
#' @return A data.frame, or a list with `syn`, and (depending on the
#'   arguments) `model`, `param` and `tune`.
#' @examples
#' \dontrun{
#' syn_df <- synth_ctgan(penguins, discrete_columns = c("species", "island", "sex"))
#' }
#' @export
synth_ctgan <- function(real_df,
                        discrete_columns,
                        tune = FALSE,
                        tune_grid = NULL,
                        return_tune = TRUE,
                        tune_eval = "MWS",
                        seed = NULL,
                        embedding_dim = 128L,
                        generator_dim = c(256L, 256L),
                        discriminator_dim = c(256L, 256L),
                        generator_lr = 2e-4,
                        generator_decay = 1e-6,
                        discriminator_lr = 2e-4,
                        discriminator_decay = 1e-6,
                        batch_size = 500L,
                        discriminator_steps = 1L,
                        log_frequency = TRUE,
                        verbose = FALSE,
                        epochs = 300L,
                        pac = 10L,
                        enable_gpu = TRUE,
                        return_model = FALSE) {
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
  .fit_and_sample <- function(epochs, batch_size, g_lr, d_lr, g_dim, d_dim, emb_dim,
                              g_decay, d_decay, d_steps, log_frequency, verbose,
                              pac, enable_gpu) {
    model <- ctgan$CTGAN(
      embedding_dim       = as.integer(emb_dim),
      generator_dim       = py_builtin$tuple(as.integer(g_dim)),
      discriminator_dim   = py_builtin$tuple(as.integer(d_dim)),
      generator_lr        = g_lr,
      generator_decay     = g_decay,
      discriminator_lr    = d_lr,
      discriminator_decay = d_decay,
      batch_size          = as.integer(batch_size),
      discriminator_steps = as.integer(d_steps),
      log_frequency       = log_frequency,
      verbose             = verbose,
      epochs              = as.integer(epochs),
      pac                 = as.integer(pac),
      enable_gpu          = enable_gpu
    )
    model$fit(py_df, discrete_cols)
    syn_df <- .as_r_df(model$sample(as.integer(nrow(real_df))))
    list(synth = syn_df, model = model)
  }
  
  # single run
  if (!tune) {
    res <- .fit_and_sample(epochs, batch_size, generator_lr, discriminator_lr,
                           generator_dim, discriminator_dim, embedding_dim,
                           generator_decay, discriminator_decay, discriminator_steps,
                           log_frequency, verbose, pac, enable_gpu)
    if (return_model) return(list(syn = res$synth, model = res$model))
    return(res$synth)
  }
  
  # tuning
  if (is.null(tune_grid)) {
    tune_grid <- expand.grid(
      epochs = 300L, batch = c(100L, 500L),
      glr = c(1e-4, 2e-4, 1e-3), dlr = c(1e-4, 2e-4, 1e-3),
      gdim = I(list(c(256L, 256L), c(128L, 128L), c(64L, 64L))),
      ddim = I(list(c(256L, 256L), c(128L, 128L), c(64L, 64L))),
      embed = c(64L, 128L), g_decay = 1e-6, d_decay = 1e-6, d_steps = 1L,
      log_frequency = TRUE, verbose = FALSE, pac = 10L, enable_gpu = TRUE,
      KEEP.OUT.ATTRS = FALSE)
  }
  
  fit_row <- function(p) {
    .fit_and_sample(as.integer(p$epochs), as.integer(p$batch), p$glr, p$dlr,
                    p$gdim[[1]], p$ddim[[1]], as.integer(p$embed),
                    p$g_decay, p$d_decay, p$d_steps, p$log_frequency,
                    p$verbose, p$pac, p$enable_gpu)
  }
  
  best <- .tune_loop(tune_grid, fit_row, real_df, numeric_cols,
                     discrete_columns, tune_eval)
  .tune_output(best, return_model, return_tune)
}
