#' Fit a user-supplied supervised model on the evidence feature matrix
#'
#' Wraps `glmnet::cv.glmnet(family = "binomial")` to fit a sparse
#' logistic model that ranks genes by feature-space evidence. The
#' fit is **owned entirely by the user**: labels must be supplied
#' explicitly. The package ships no pretrained model -- this is a
#' deliberate constraint for PDAC biomarker discovery, where
#' validated non-circular ground truth is unavailable.
#'
#' Use case: a user has their own positive set (e.g. their lab's
#' validated panel) and wants a continuous prioritisation score
#' across the rest of the atlas, with feature-level
#' interpretability that survives Bioconductor review.
#'
#' @section Reproducibility:
#' RNG is scoped via [withr::local_seed()] inside the fit. Two
#' calls with the same `(features, labels, alpha, nfolds, seed)`
#' produce byte-identical models.
#'
#' @section Cross-validation:
#' 10-fold CV by `glmnet::cv.glmnet()` over the lambda path; AUC
#' is recomputed per fold for the model card. Class imbalance is
#' not reweighted by default (use the `weights` arg if needed).
#'
#' @param features `data.table` from [make_evidence_features()] (or
#'   any frame with `gene_symbol` plus numeric feature columns).
#'   Z-scored input (`scale = "z"`) is recommended for stable
#'   coefficients.
#' @param labels Integer / logical vector of length `nrow(features)`,
#'   one entry per row, with `1` / `TRUE` for positives and
#'   `0` / `FALSE` for negatives.
#' @param method Currently only `"elastic_net"` is implemented.
#'   Argument exists so a future release can add `"rf"` /
#'   `"gbm"` without API churn.
#' @param alpha Elastic-net mixing parameter in `[0, 1]`. `0` is
#'   ridge, `1` is lasso. Default `0.5` (mid-point) keeps non-zero
#'   coefficients while damping correlated features.
#' @param nfolds Number of CV folds. Default `10`.
#' @param weights Optional per-row observation weights passed
#'   through to `glmnet::cv.glmnet()`. Default `NULL`
#'   (uniform weights).
#' @param seed Integer RNG seed for reproducibility. Default `1L`.
#' @return An S3 object of class `pdactrace_user_model` (a list)
#'   with elements: `cvfit` (the underlying `cv.glmnet` object),
#'   `coef_table` (data.table: `feature`, `coef_value`),
#'   `lambda_min`, `lambda_1se`, `cv_auc_mean`, `cv_auc_sd`,
#'   `n_train`, `n_features`, `n_positives`, `top_pos`,
#'   `top_neg`, `alpha`, `seed`, `method`, `feature_set_version`,
#'   `card`.
#' @examples
#' \dontrun{
#'   if (requireNamespace("glmnet", quietly = TRUE)) {
#'     feats <- make_evidence_features(scale = "z",
#'                                       drop_na_rows = TRUE)
#'     y <- as.integer(feats$gene_symbol %in%
#'                       c("LGALS3BP", "TIMP1", "LRG1"))
#'     fit <- fit_user_evidence_model(feats, y)
#'     explain_user_evidence_model(fit, top_n = 5)
#'   }
#' }
#' @seealso [predict_user_evidence_model()],
#'   [explain_user_evidence_model()], [model_card()],
#'   [make_evidence_features()].
#' @export
fit_user_evidence_model <- function(features,
                                     labels,
                                     method  = "elastic_net",
                                     alpha   = 0.5,
                                     nfolds  = 10L,
                                     weights = NULL,
                                     seed    = 1L) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("`glmnet` is required for fit_user_evidence_model(). ",
         "Install via install.packages(\"glmnet\").",
         call. = FALSE)
  }
  if (!identical(method, "elastic_net")) {
    stop("Only `method = 'elastic_net'` is supported in v0.99.6.",
         call. = FALSE)
  }
  if (is.null(labels)) {
    stop("`labels` must be supplied; no shipped supervised model.",
         call. = FALSE)
  }
  if (!is.data.frame(features)) {
    stop("`features` must be a data.frame / data.table.",
         call. = FALSE)
  }

  feat_dt <- data.table::as.data.table(features)
  if (!"gene_symbol" %in% names(feat_dt)) {
    stop("`features` must include a `gene_symbol` column.",
         call. = FALSE)
  }
  feat_cols <- setdiff(names(feat_dt), "gene_symbol")
  X <- as.matrix(feat_dt[, feat_cols, with = FALSE])
  y <- as.integer(as.logical(labels))
  if (length(y) != nrow(X)) {
    stop("length(labels) must equal nrow(features); got ",
         length(y), " vs ", nrow(X), ".", call. = FALSE)
  }
  if (!all(y %in% c(0L, 1L))) {
    stop("`labels` must coerce to {0, 1}; got values: ",
         paste(unique(y), collapse = ", "), ".", call. = FALSE)
  }
  n_pos <- sum(y == 1L)
  if (n_pos == 0L || n_pos == length(y)) {
    stop("`labels` must contain both positives and negatives.",
         call. = FALSE)
  }

  # complete-case rows only (cv.glmnet does not accept NAs)
  cc <- stats::complete.cases(X)
  if (sum(cc) < 20L) {
    stop("Fewer than 20 complete-case rows after NA-drop ",
         "(", sum(cc), "); refusing to fit.", call. = FALSE)
  }
  Xcc <- X[cc, , drop = FALSE]
  ycc <- y[cc]
  wcc <- if (!is.null(weights)) as.numeric(weights)[cc] else NULL
  if (sum(ycc == 1L) < 2L || sum(ycc == 0L) < 2L) {
    stop("After NA-drop, need >= 2 positives and >= 2 negatives.",
         call. = FALSE)
  }

  cvfit <- withr::with_seed(seed, {
    if (is.null(wcc)) {
      glmnet::cv.glmnet(Xcc, ycc, family = "binomial",
                        alpha = alpha, nfolds = nfolds,
                        type.measure = "auc")
    } else {
      glmnet::cv.glmnet(Xcc, ycc, family = "binomial",
                        alpha = alpha, nfolds = nfolds,
                        type.measure = "auc", weights = wcc)
    }
  })

  # CV AUC at lambda.min
  cvm  <- cvfit$cvm
  cvsd <- cvfit$cvsd
  idx_min <- which(cvfit$lambda == cvfit$lambda.min)
  cv_auc_mean <- if (length(idx_min)) cvm[idx_min] else NA_real_
  cv_auc_sd   <- if (length(idx_min)) cvsd[idx_min] else NA_real_

  # Coefficients at lambda.min
  coefs_sp <- stats::coef(cvfit, s = "lambda.min")
  coef_dt <- data.table::data.table(
    feature    = rownames(coefs_sp),
    coef_value = as.numeric(coefs_sp))
  coef_dt <- coef_dt[feature != "(Intercept)"]
  coef_dt <- coef_dt[order(-abs(coef_value))]

  pos <- coef_dt[coef_value > 0][order(-coef_value)]
  neg <- coef_dt[coef_value < 0][order(coef_value)]

  out <- list(
    cvfit                = cvfit,
    coef_table           = coef_dt,
    lambda_min           = cvfit$lambda.min,
    lambda_1se           = cvfit$lambda.1se,
    cv_auc_mean          = cv_auc_mean,
    cv_auc_sd            = cv_auc_sd,
    n_train              = nrow(Xcc),
    n_features           = ncol(Xcc),
    n_positives          = sum(ycc == 1L),
    top_pos              = head(pos, 10L),
    top_neg              = head(neg, 10L),
    alpha                = alpha,
    seed                 = seed,
    method               = method,
    feature_set_version  = attr(features,
                                  "feature_set_version") %||% "v1.0")
  out$card <- model_card(source = "user_model", model =
                            structure(out,
                                      class = "pdactrace_user_model"),
                          verbose = FALSE)
  class(out) <- "pdactrace_user_model"
  out
}

#' Predict from a `pdactrace_user_model`
#'
#' @param model A `pdactrace_user_model` from
#'   [fit_user_evidence_model()].
#' @param new_features `data.table` from [make_evidence_features()]
#'   (or a frame with `gene_symbol` and the same numeric feature
#'   columns the model was trained on).
#' @param s `lambda` value at which to evaluate. Default
#'   `"lambda.min"`. Passed through to `predict.cv.glmnet()`.
#' @return A `data.table` with `gene_symbol` and `predicted_prob`
#'   (P(label=1) under the fitted model). Rows with any `NA`
#'   feature are returned with `predicted_prob = NA`.
#' @examples
#' \dontrun{
#'   feats <- make_evidence_features(scale = "z",
#'                                     drop_na_rows = TRUE)
#'   y <- as.integer(feats$gene_symbol %in%
#'                     c("LGALS3BP", "TIMP1", "LRG1"))
#'   fit <- fit_user_evidence_model(feats, y)
#'   head(predict_user_evidence_model(fit, feats))
#' }
#' @export
predict_user_evidence_model <- function(model,
                                         new_features,
                                         s = "lambda.min") {
  if (!inherits(model, "pdactrace_user_model")) {
    stop("`model` must be a `pdactrace_user_model`.", call. = FALSE)
  }
  if (!"gene_symbol" %in% names(new_features)) {
    stop("`new_features` must include `gene_symbol`.",
         call. = FALSE)
  }
  feats <- data.table::as.data.table(new_features)
  feat_cols <- setdiff(names(feats), "gene_symbol")

  # Align column order with what the model was trained on
  trained_cols <- rownames(stats::coef(model$cvfit,
                                         s = "lambda.min"))
  trained_cols <- setdiff(trained_cols, "(Intercept)")
  miss <- setdiff(trained_cols, feat_cols)
  if (length(miss) > 0L) {
    stop("`new_features` is missing trained feature columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  X <- as.matrix(feats[, trained_cols, with = FALSE])
  cc <- stats::complete.cases(X)

  prob <- rep(NA_real_, nrow(X))
  if (any(cc)) {
    pr <- stats::predict(model$cvfit, newx = X[cc, , drop = FALSE],
                          s = s, type = "response")
    prob[cc] <- as.numeric(pr)
  }
  data.table::data.table(gene_symbol    = feats$gene_symbol,
                         predicted_prob = prob)[]
}

#' Explain a `pdactrace_user_model` via its top contributing features
#'
#' Prints (and invisibly returns) the top-N positive and top-N
#' negative coefficients of the fitted elastic-net model, with
#' feature names and signed coefficient values. This is the
#' user-facing rationale for the supervised score: which evidence
#' axes pushed the user's positive set up.
#'
#' @param model A `pdactrace_user_model` from
#'   [fit_user_evidence_model()].
#' @param top_n Number of positive AND number of negative
#'   coefficients to surface. Default `10`.
#' @param verbose Logical. If `TRUE` (default), prints to console.
#' @return Invisibly, a list with `top_pos`, `top_neg`,
#'   `cv_auc_mean`, `cv_auc_sd`, `lambda_min`, `n_features`,
#'   `n_train`, `n_positives`.
#' @examples
#' \dontrun{
#'   feats <- make_evidence_features(scale = "z",
#'                                     drop_na_rows = TRUE)
#'   y <- as.integer(feats$gene_symbol %in%
#'                     c("LGALS3BP", "TIMP1", "LRG1"))
#'   fit <- fit_user_evidence_model(feats, y)
#'   explain_user_evidence_model(fit, top_n = 5)
#' }
#' @seealso [model_card()] for the full reviewer-facing card.
#' @export
explain_user_evidence_model <- function(model,
                                         top_n = 10L,
                                         verbose = TRUE) {
  if (!inherits(model, "pdactrace_user_model")) {
    stop("`model` must be a `pdactrace_user_model`.", call. = FALSE)
  }
  pos <- head(model$top_pos, top_n)
  neg <- head(model$top_neg, top_n)
  if (isTRUE(verbose)) {
    cat("User evidence model -- explain (top ", top_n,
        " each side)\n", sep = "")
    cat(strrep("-", 60), "\n", sep = "")
    cat(sprintf("  CV AUC: %.3f +/- %.3f at lambda.min = %.4f\n",
                model$cv_auc_mean, model$cv_auc_sd,
                model$lambda_min))
    cat(sprintf("  fit on n=%d genes (%d positives, %d features)\n\n",
                model$n_train, model$n_positives, model$n_features))
    cat("Main positive contributors (push score UP):\n")
    if (nrow(pos) == 0L) cat("  (none -- all coefficients <= 0)\n")
    for (i in seq_len(nrow(pos))) {
      cat(sprintf("  +%.3f  %s\n",
                  pos$coef_value[i], pos$feature[i]))
    }
    cat("\nMain negative contributors (push score DOWN):\n")
    if (nrow(neg) == 0L) cat("  (none -- all coefficients >= 0)\n")
    for (i in seq_len(nrow(neg))) {
      cat(sprintf("  %.3f  %s\n",
                  neg$coef_value[i], neg$feature[i]))
    }
  }
  invisible(list(top_pos     = pos,
                 top_neg     = neg,
                 cv_auc_mean = model$cv_auc_mean,
                 cv_auc_sd   = model$cv_auc_sd,
                 lambda_min  = model$lambda_min,
                 n_features  = model$n_features,
                 n_train     = model$n_train,
                 n_positives = model$n_positives))
}

#' @export
print.pdactrace_user_model <- function(x, ...) {
  cat("<pdactrace_user_model>\n")
  cat("  method:      ", x$method, "\n", sep = "")
  cat("  n_train:     ", x$n_train, "  (positives = ",
      x$n_positives, ")\n", sep = "")
  cat("  n_features:  ", x$n_features, "\n", sep = "")
  cat(sprintf("  CV AUC:      %.3f +/- %.3f at lambda.min = %.4f\n",
              x$cv_auc_mean, x$cv_auc_sd, x$lambda_min))
  cat("  alpha:       ", x$alpha,
      "  seed: ", x$seed, "\n", sep = "")
  cat("\nUse explain_user_evidence_model() to see top contributors\n",
      "or predict_user_evidence_model() to score new genes.\n",
      sep = "")
  invisible(x)
}
