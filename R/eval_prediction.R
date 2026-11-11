# Prediction based evaluators: how well can a classifier tell real and
# synthetic rows apart? An AUC close to 0.5 means they are hard to tell apart.

# Logistic regression (ridge)

#' Distinguish real from synthetic data with logistic regression
#'
#' Fits a penalised logistic regression (glmnet) that predicts whether a row
#' is real or synthetic.
#'
#' @param real_df,syn_df data.frames with the real and the synthetic data
#'   (same column names).
#' @param test_size Fraction of rows used as test set (only used when
#'   `repeatedcv = FALSE`).
#' @param seed Random seed.
#' @param alpha Elastic net mixing parameter (0 = ridge, 1 = lasso).
#' @param n_folds Number of cross-validation folds.
#' @param n_repeat Number of repeats of the cross-validation.
#' @param repeatedcv If `TRUE`, use repeated cross-validation (caret). If
#'   `FALSE`, use a single train/test split.
#' @param digits Number of digits that are printed.
#' @param measure Measure that is printed when `repeatedcv = FALSE`: `"auc"`,
#'   `"acc"`, `"deviance"`, `"class"`, `"mse"`, `"mae"` or `"all"`.
#' @param max_unique_for_cat Numeric columns with at most this many unique
#'   values are treated as categorical (see [detect_types()]).
#' @return If `repeatedcv = TRUE`, a list with `fit`, `auc_mean`, `auc_sd`,
#'   `acc_mean` and `acc_sd`. Otherwise the output of
#'   [glmnet::assess.glmnet()].
#' @examples
#' \dontrun{
#' eval_logreg(real_df, syn_df)
#' }
#' @export
eval_logreg <- function(real_df, syn_df,
                        test_size = 0.30,
                        seed = 123,
                        alpha = 0,      # ridge
                        n_folds = 10,
                        n_repeat = 10,
                        repeatedcv = TRUE,
                        digits = 3,
                        measure = c("auc"),
                        max_unique_for_cat = 2) # 
{
  if (!measure %in% c("auc", "acc", "deviance", "class", "mse", "mae" , "all")) {
    stop("Measure must be one of: auc, acc, deviance, class, mse, mae, all ")}
  
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
  
  # add labels and combines df's
  df <- rbind(
    transform(real_df, label = "real"),
    transform(syn_df,  label = "fake")
  )
  df$label <- factor(df$label, levels = c("real", "fake"))
  
  # design matrix once on full data (avoids train/test dummy mismatches)
  X_full <- model.matrix(label ~ . - 1, data = df)
  y_full <- factor(df$label, levels = c("real", "fake"))
  
  # shuffle & split
  set.seed(seed)
  idx <- sample.int(nrow(X_full))
  X_full <- X_full[idx, , drop = FALSE]
  y_full <- y_full[idx]
  
  n_train <- floor((1 - test_size) * nrow(X_full))
  X_train <- X_full[seq_len(n_train), , drop = FALSE]
  y_train <- y_full[seq_len(n_train)]
  X_test  <- X_full[-seq_len(n_train), , drop = FALSE]
  y_test  <- y_full[-seq_len(n_train)]
  
  # print metrics
  if (repeatedcv) {
    mySummary <- function(data, lev = NULL, model = NULL) {
      roc_stats <- caret::twoClassSummary(data, lev = lev, model = model)
      acc_stats <- caret::defaultSummary(data, lev = lev, model = model)
      c(ROC = roc_stats["ROC"], Accuracy = acc_stats["Accuracy"])
    }
    
    ctrl <- trainControl(
      method = "repeatedcv",
      number = n_folds,
      repeats = n_repeat,
      classProbs = TRUE,
      summaryFunction = mySummary,
      savePredictions = "final"
    )
    lambda_grid <- 10^seq(-4, 1, length.out = 50)
    tune_grid <- expand.grid(alpha = alpha, lambda = lambda_grid    )
    
    set.seed(seed)
    
    fit <- caret::train(
      label ~ .,
      data = df,
      method = "glmnet",
      family = "binomial",
      metric = "ROC",
      trControl = ctrl,
      tuneGrid = tune_grid
    )
    
    auc_mean <- mean(fit$resample$ROC)
    auc_sd   <- sd(fit$resample$ROC)
    acc_mean <- mean(fit$resample$Accuracy)
    acc_sd   <- sd(fit$resample$Accuracy)
    
    cat("AUC (repeated CV): ", round(auc_mean, digits), "\u00B1",round(auc_sd, digits), "\n")
    cat("Acc (repeated CV): ", round(acc_mean, digits), "\u00B1",round(acc_sd, digits), "\n")
    
    return(list(fit = fit,
                auc_mean = auc_mean,
                auc_sd = auc_sd,
                acc_mean = acc_mean,
                acc_sd = acc_sd))
  } else {
    # cross val ridge logistic regression
    cvfit  <- cv.glmnet(X_train, y_train, family = "binomial", alpha = alpha, nfolds = n_folds)
    lambda <- cvfit$lambda.min
    # predict prob(real) and threshold at 0.5
    p    <- as.numeric(predict(cvfit, newx = X_test, s = lambda, type = "response"))
    yhat <- ifelse(p > 0.5, "fake", "real")
    acc <- mean(yhat == y_test)
    
    assess_fit <- assess.glmnet(cvfit, newx = X_test, newy = y_test)
    
    if (measure == "auc" || measure == "all") {cat("AUC is ", assess_fit$auc, "\n") }
    if (measure == "acc" || measure == "all") {cat("Accuracy:", round(acc, digits), "\n") }
    if (measure == "deviance" || measure == "all") {cat("Binomial Deviance is ", round(assess_fit$deviance, digits), "\n") }
    if (measure == "class" || measure == "all") {cat("Misclassification Errror is ", round(assess_fit$class, digits), "\n") }
    if (measure == "mse" || measure == "all") {cat("Mean-Squared Error is ", round(assess_fit$mse, digits), "\n") }
    if (measure == "mae" || measure == "all") {cat("Mean Absolute Error is ", round(assess_fit$mae, digits), "\n") }
    return(assess_fit)
  }
  
}

# Random forest

#' Distinguish real from synthetic data with a random forest
#'
#' Fits a random forest (caret, repeated cross-validation) that predicts
#' whether a row is real or synthetic.
#'
#' @param real_df,syn_df data.frames with the real and the synthetic data
#'   (same column names).
#' @param n_repeat Number of repeats of the cross-validation.
#' @param n_folds Number of cross-validation folds.
#' @param seed Random seed.
#' @param ntree Number of trees.
#' @param digits Number of digits that are printed.
#' @param max_unique_for_cat Numeric columns with at most this many unique
#'   values are treated as categorical (see [detect_types()]).
#' @return A list with a data.frame (`AUC` and `Accuracy` per resample) and
#'   the fitted caret model.
#' @examples
#' \dontrun{
#' eval_randomForest(real_df, syn_df)
#' }
#' @export
eval_randomForest <- function (real_df, 
                               syn_df, 
                               # test_size = 0.30,
                               # repeated = TRUE,
                               n_repeat = 10,
                               n_folds = 10,
                               seed = 123,
                               ntree = 500,
                               digits = 3,
                               max_unique_for_cat = 2) {
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
  
  # combine & add labels
  df <- rbind(
    transform(real_df, label = "real"),
    transform(syn_df,  label = "fake")
  )
  df$label <- factor(df$label, levels = c("real","fake"))
  
  # shuffle & split
  set.seed(seed)

  train_cv <- caret::trainControl(
    method = "repeatedcv",
    number = n_folds,
    repeats = n_repeat,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final"
  )
  
  # random forest
  rf_model <- caret::train(
    label ~ .,
    data = df,
    method = "rf",
    ntree = ntree,
    metric = "ROC",
    trControl = train_cv,
    tuneGrid = expand.grid(mtry = ceiling(c(sqrt(ncol(df)) * c(0.5, 1, 2))))
  )
  
  best <- rf_model$bestTune
  preds <- rf_model$pred
  
  auc_by_fold <- preds %>%
    dplyr::filter(mtry == best$mtry) %>%
    dplyr::group_by(Resample) %>%
    dplyr::summarise(
      AUC = as.numeric(pROC::roc(obs, fake, levels = c("real", "fake"), direction = "<", quiet = TRUE)$auc),
      Accuracy = mean(pred == obs)
    )
  
  auc_mean <- mean(auc_by_fold$AUC)
  auc_sd   <- sd(auc_by_fold$AUC)
  
  acc_mean <- mean(auc_by_fold$Accuracy)
  acc_sd   <- sd(auc_by_fold$Accuracy)
  
  cat("RF AUC:", round(auc_mean, digits), "\u00B1", round(auc_sd, digits), "\n")
  cat("RF Accuracy:", round(acc_mean, digits), "\u00B1", round(acc_sd, digits), "\n")
  
  return(list(data.frame(AUC = auc_by_fold$AUC, Accuracy = auc_by_fold$Accuracy), rf_model))
  
}
