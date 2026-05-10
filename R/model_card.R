#' Model card metadata for the optional ML layer
#'
#' Single point of documentation reviewers can audit when a
#' descriptive ([score_anchor_similarity()]) or supervised
#' ([fit_user_evidence_model()]) score is in use. Returns a
#' structured list and prints a multi-line summary to the console.
#'
#' Two sources are supported:
#' \describe{
#'   \item{`"anchor_similarity"`}{The descriptive layer. No labels
#'     are used to train a predictor; the bundled
#'     [pdactrace_external_anchors] table is used **only** to
#'     define a reference centroid in the z-scored feature space.
#'     Always available, no `model` argument needed.}
#'   \item{`"user_model"`}{The user-supplied supervised fit
#'     produced by [fit_user_evidence_model()]. Pass the fit
#'     object via `model = ...`.}
#' }
#'
#' The package **ships no pretrained predictor**. This deliberate
#' constraint is the central defense against the
#' "what was the training label?" reviewer question for PDAC
#' biomarker discovery, where validated non-circular ground truth
#' is not available.
#'
#' @param source One of `"anchor_similarity"` or `"user_model"`.
#' @param model Optional `pdactrace_user_model` from
#'   [fit_user_evidence_model()] (required when
#'   `source = "user_model"`).
#' @param verbose Logical. If `TRUE` (default), prints the card to
#'   the console. The structured list is returned invisibly either
#'   way.
#' @return Invisibly, a list with at least: `source`,
#'   `feature_set_version`, `intended_use`, `not_intended_for`,
#'   `labels_used_for_training`, `limitations`. For
#'   `"user_model"`, also includes `n_features`, `n_train`,
#'   `n_positives`, `cv_auc_mean`, `cv_auc_sd`, `top_pos`,
#'   `top_neg`, `seed`, `alpha`, `lambda_min`.
#' @examples
#' model_card("anchor_similarity")
#' \dontrun{
#'   if (requireNamespace("glmnet", quietly = TRUE)) {
#'     feats <- make_evidence_features(scale = "z",
#'                                       drop_na_rows = TRUE)
#'     y <- as.integer(feats$gene_symbol %in%
#'                       c("LGALS3BP", "TIMP1", "LRG1"))
#'     fit <- fit_user_evidence_model(feats, y)
#'     model_card(source = "user_model", model = fit)
#'   }
#' }
#' @seealso [score_anchor_similarity()],
#'   [fit_user_evidence_model()].
#' @export
model_card <- function(source = c("anchor_similarity", "user_model"),
                        model = NULL,
                        verbose = TRUE) {
  source <- match.arg(source)

  card <- if (source == "anchor_similarity") {
    .mc_anchor_card()
  } else {
    if (is.null(model) || !inherits(model, "pdactrace_user_model")) {
      stop("`model` must be a `pdactrace_user_model` from ",
           "fit_user_evidence_model() when source='user_model'.",
           call. = FALSE)
    }
    .mc_user_card(model)
  }

  if (isTRUE(verbose)) .mc_print(card)
  invisible(card)
}

# ---- builders ----------------------------------------------------------

.mc_anchor_card <- function() {
  list(
    source                   = "anchor_similarity",
    feature_set_version      = "v1.0",
    intended_use             = paste0(
      "Descriptive feature-space proximity to a curated set of ",
      "external positive anchors, for candidate prioritisation."),
    not_intended_for         = paste0(
      "Clinical diagnosis, regulatory submissions, or as a ",
      "supervised classifier output."),
    labels_used_for_training = paste0(
      "None. Anchor labels are used only to define the reference ",
      "centroid in the z-scored feature space; no supervised ",
      "model is fit and no predictor is shipped."),
    leakage_controls         = paste0(
      "Anchors are evaluation-only by package convention; the same ",
      "anchor set is also consumed by ",
      "evaluate_anchor_enrichment() in the post-hoc evaluator. ",
      "If you also use this similarity score for ranking you must ",
      "report enrichment computed on a held-out anchor split."),
    validation_summary       = paste0(
      "T1_validated anchors should land in the top decile by ",
      "anchor_similarity (verified by test-anchor-similarity.R)."),
    limitations              = c(
      "Single bundled anchor table (n=30) limits centroid stability.",
      "Equal weighting per feature; no learned feature weights.",
      "Cosine similarity is direction-only; large effect-size genes ",
      "may not rank above weak-signal anchors with similar shape."))
}

.mc_user_card <- function(m) {
  list(
    source                   = "user_model",
    feature_set_version      = m$feature_set_version %||% "v1.0",
    method                   = m$method %||% "elastic_net",
    intended_use             = paste0(
      "Candidate prioritisation under labels supplied by the user. ",
      "The fit is produced and owned entirely by the user; the ",
      "package ships no pretrained model."),
    not_intended_for         = paste0(
      "Generalisation beyond the user's labelled cohort, clinical ",
      "diagnosis, or regulatory submissions."),
    labels_used_for_training = paste0(
      "User-supplied binary labels (n=", m$n_train, "; ",
      m$n_positives, " positives, ", m$n_train - m$n_positives,
      " negatives)."),
    n_features               = m$n_features,
    n_train                  = m$n_train,
    n_positives              = m$n_positives,
    cv_auc_mean              = m$cv_auc_mean,
    cv_auc_sd                = m$cv_auc_sd,
    alpha                    = m$alpha,
    lambda_min               = m$lambda_min,
    seed                     = m$seed,
    top_pos                  = m$top_pos,
    top_neg                  = m$top_neg,
    leakage_controls         = paste0(
      "10-fold cross-validation by glmnet::cv.glmnet(). The user ",
      "must verify their label set is independent from any held-out ",
      "evaluation set used downstream."),
    limitations              = c(
      paste0("Linear (logistic) decision boundary; no learned ",
             "feature interactions."),
      "Class imbalance is not reweighted by default.",
      paste0("Fit is reproducible only at the same `seed` and ",
             "`alpha`.")))
}

# ---- printer -----------------------------------------------------------

.mc_print <- function(card) {
  cat("Model card -- ", card$source, "\n", sep = "")
  cat(strrep("-", 60), "\n", sep = "")
  ord <- c("source", "feature_set_version", "method",
           "intended_use", "not_intended_for",
           "labels_used_for_training",
           "n_features", "n_train", "n_positives",
           "cv_auc_mean", "cv_auc_sd",
           "alpha", "lambda_min", "seed",
           "leakage_controls", "validation_summary",
           "top_pos", "top_neg",
           "limitations")
  for (k in ord) {
    if (!k %in% names(card)) next
    v <- card[[k]]
    if (is.null(v)) next
    if (k %in% c("limitations") && length(v) > 1L) {
      cat(sprintf("  %s:\n", k))
      for (li in v) cat("    - ", li, "\n", sep = "")
      next
    }
    if (k %in% c("top_pos", "top_neg") &&
        data.table::is.data.table(v)) {
      cat(sprintf("  %s:\n", k))
      print(v)
      next
    }
    if (is.numeric(v) && length(v) == 1L) {
      cat(sprintf("  %-25s %s\n", paste0(k, ":"),
                  formatC(v, format = "g", digits = 4)))
    } else if (length(v) == 1L) {
      cat(sprintf("  %-25s %s\n", paste0(k, ":"), as.character(v)))
    }
  }
  invisible(NULL)
}

# null-coalesce (re-defined locally to avoid load-order dependency on
# format_provenance.R's `%||%`).
`%||%` <- function(a, b) if (is.null(a)) b else a
