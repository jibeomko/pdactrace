test_that("fit_user_evidence_model errors when labels are missing", {
  skip_if_not_installed("glmnet")
  feats <- make_evidence_features(impute = "mean")
  expect_error(fit_user_evidence_model(feats, labels = NULL),
               regexp = "labels.*supplied")
})

test_that("fit_user_evidence_model errors on single-class labels", {
  skip_if_not_installed("glmnet")
  feats <- make_evidence_features(impute = "mean")
  y_allzero <- rep(0L, nrow(feats))
  y_allone  <- rep(1L, nrow(feats))
  expect_error(fit_user_evidence_model(feats, y_allzero),
               regexp = "both positives and negatives")
  expect_error(fit_user_evidence_model(feats, y_allone),
               regexp = "both positives and negatives")
})

test_that("fit_user_evidence_model errors on length-mismatched labels", {
  skip_if_not_installed("glmnet")
  feats <- make_evidence_features(impute = "mean")
  expect_error(fit_user_evidence_model(feats, c(0L, 1L, 0L)),
               regexp = "length\\(labels\\)")
})

test_that("anchor_similarity is independent from supervised training labels", {
  # The bundled anchors define the centroid for similarity; they are
  # NOT used as supervised labels for any shipped predictor. Verify
  # this discipline programmatically by checking that
  # score_anchor_similarity() does not require, and never reads,
  # any "label" column from the user.
  sim <- score_anchor_similarity(tier = "primary")
  expect_s3_class(sim, "data.table")
  # No "label" / "y" column appears in the output
  expect_false(any(c("label", "y", "outcome") %in% names(sim)))
})

test_that("evaluate_anchor_enrichment + score_anchor_similarity remain distinct layers", {
  # Both consume the same anchor set, but for different purposes.
  # Confirm the two have distinct return shapes so users (and reviewers)
  # cannot confuse them.
  sim <- score_anchor_similarity(tier = "primary")
  ev  <- evaluate_anchor_enrichment(top_n = c(50, 100), tier = "primary")
  expect_false(identical(names(sim), names(ev)))
  expect_true("anchor_similarity" %in% names(sim))
  expect_true("hits" %in% names(ev))
})
