test_that("model_card('anchor_similarity') returns required fields", {
  card <- model_card("anchor_similarity", verbose = FALSE)
  expect_type(card, "list")
  for (k in c("source", "feature_set_version", "intended_use",
              "not_intended_for", "labels_used_for_training",
              "leakage_controls", "limitations")) {
    expect_true(k %in% names(card), info = k)
  }
  expect_equal(card$source, "anchor_similarity")
  expect_match(card$labels_used_for_training,
               "no supervised model is fit", fixed = TRUE)
})

test_that("model_card('user_model') errors without a fit", {
  expect_error(model_card("user_model", model = NULL,
                           verbose = FALSE),
               regexp = "pdactrace_user_model")
})

test_that("model_card('user_model') populates CV summary from the fit", {
  skip_if_not_installed("glmnet")
  feats <- make_evidence_features(scale = "z", impute = "mean")
  data("pdactrace_external_anchors", package = "pdactrace")
  t1 <- pdactrace_external_anchors$gene[
    pdactrace_external_anchors$evidence_tier == "T1_validated"]
  y <- as.integer(feats$gene_symbol %in% t1)
  fit <- suppressWarnings(
    fit_user_evidence_model(feats, y, alpha = 0.5, seed = 1L))
  card <- model_card("user_model", model = fit, verbose = FALSE)
  for (k in c("cv_auc_mean", "cv_auc_sd", "lambda_min", "alpha",
              "seed", "n_train", "n_features", "n_positives",
              "top_pos", "top_neg", "leakage_controls",
              "limitations")) {
    expect_true(k %in% names(card), info = k)
  }
  expect_equal(card$source, "user_model")
})

test_that("printing the card produces output without error", {
  out <- capture.output(model_card("anchor_similarity"),
                         type = "output")
  joined <- paste(out, collapse = "\n")
  expect_match(joined, "Model card", fixed = TRUE)
  expect_match(joined, "anchor_similarity", fixed = TRUE)
  expect_match(joined, "no supervised model is fit", fixed = TRUE)
})
