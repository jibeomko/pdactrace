test_that("explain_score returns the expected breakdown structure", {
  res <- explain_score("LGALS3BP", verbose = FALSE)
  expect_type(res, "list")
  expect_named(res, c("gene", "audit_class", "audit_score",
                       "positive_score", "axes", "gates", "explanation"),
               ignore.order = TRUE)
  expect_s3_class(res$axes, "data.table")
  expect_equal(nrow(res$axes), 3L)
  expect_equal(sum(res$axes$weight), 1, tolerance = 1e-6)
  expect_equal(nrow(res$gates), 2L)
  expect_match(res$explanation, res$audit_class, fixed = TRUE)
})

test_that("explain_score errors clearly on missing gene", {
  expect_error(explain_score("NOT_A_REAL_GENE_XYZ", verbose = FALSE),
               regexp = "not in the bundled atlas")
})

test_that("compare_candidates returns one row per input gene", {
  cmp <- compare_candidates(c("LGALS3BP", "LTBP1", "ALB", "GAPDH"))
  expect_s3_class(cmp, "data.table")
  expect_equal(nrow(cmp), 4L)
  expect_true("redundancy_with" %in% names(cmp))
  # audit_score sort: NA last, descending otherwise
  scores <- cmp$audit_score
  expect_true(all(diff(stats::na.omit(scores)) <= 0))
})

test_that("compare_candidates pads NA rows for genes outside the atlas", {
  cmp <- compare_candidates(c("LGALS3BP", "NOT_A_REAL_GENE_XYZ"))
  expect_equal(nrow(cmp), 2L)
  miss <- cmp[gene_symbol == "NOT_A_REAL_GENE_XYZ"]
  expect_equal(nrow(miss), 1L)
  expect_true(is.na(miss$audit_class))
})

test_that("as_summarized_experiment produces a valid SE with the right shape", {
  skip_if_not_installed("SummarizedExperiment")
  se <- as_summarized_experiment()
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(ncol(se), 4L)
  expect_equal(SummarizedExperiment::assayNames(se)[1L], "rna_beta")
  ass <- SummarizedExperiment::assay(se, "rna_beta")
  expect_equal(colnames(ass), c("Normal", "Early", "Mid", "Late"))
  # Reference column is identically zero by construction
  expect_true(all(ass[, "Normal"] == 0))
  rd <- SummarizedExperiment::rowData(se)
  expect_true("audit_class" %in% colnames(rd))
})
