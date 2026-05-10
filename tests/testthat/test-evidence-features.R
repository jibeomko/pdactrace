test_that("make_evidence_features returns a wide numeric data.table", {
  feats <- make_evidence_features()
  expect_s3_class(feats, "data.table")
  expect_true("gene_symbol" %in% names(feats))
  num_cols <- setdiff(names(feats), "gene_symbol")
  expect_gte(length(num_cols), 18L)
  expect_true(all(vapply(feats[, num_cols, with = FALSE],
                          is.numeric, logical(1L))))
})

test_that("no duplicated gene_symbol rows", {
  feats <- make_evidence_features()
  expect_equal(anyDuplicated(feats$gene_symbol), 0L)
})

test_that("trajectory_delta_rho equals rho_best - rho_runner_up from atlas", {
  data("pdactrace_reference", package = "pdactrace")
  feats <- make_evidence_features()
  for (g in c("LGALS3BP", "LTBP1", "TIMP1")) {
    target <- g
    expected <- pdactrace_reference[gene_symbol == target,
                                     rna_pattern_rho - rna_pattern_rho_runner_up]
    got <- feats[gene_symbol == target, trajectory_delta_rho]
    expect_equal(got, expected, tolerance = 1e-8)
  }
})

test_that("scale='z' produces column means ~0 and sds ~1 on complete cases", {
  feats <- make_evidence_features(scale = "z", impute = "mean")
  num_cols <- setdiff(names(feats), "gene_symbol")
  for (j in num_cols) {
    v <- feats[[j]]
    if (length(unique(stats::na.omit(v))) < 2L) next
    expect_lt(abs(mean(v, na.rm = TRUE)), 0.05, label = j)
    expect_lt(abs(stats::sd(v, na.rm = TRUE) - 1), 0.05, label = j)
  }
})

test_that("impute='mean' eliminates NAs in feature columns", {
  feats <- make_evidence_features(impute = "mean")
  num_cols <- setdiff(names(feats), "gene_symbol")
  expect_false(any(is.na(feats[, num_cols, with = FALSE])))
})

test_that("impute='zero' fills NAs with zero", {
  feats <- make_evidence_features(impute = "zero")
  num_cols <- setdiff(names(feats), "gene_symbol")
  expect_false(any(is.na(feats[, num_cols, with = FALSE])))
  # at least one column that was sparse should now contain literal zeros
  expect_true(any(feats$rna_protein_cosine == 0))
})

test_that("genes argument subsets correctly", {
  feats <- make_evidence_features(genes = c("LGALS3BP", "LTBP1"))
  expect_equal(nrow(feats), 2L)
  expect_setequal(feats$gene_symbol, c("LGALS3BP", "LTBP1"))
})

test_that("attr(feature_set_version) is set on the returned object", {
  feats <- make_evidence_features()
  expect_equal(attr(feats, "feature_set_version"), "v1.0")
  expect_equal(attr(feats, "scale"), "none")
  expect_equal(attr(feats, "impute"), "none")
})
