test_that("project_user_cohort runs end-to-end on toy data", {
  skip_if_not_installed("DESeq2")
  res <- project_user_cohort(
    rna        = toy_counts,
    coldata    = toy_coldata,
    stage_col  = "stage",
    cohort_col = "cohort")
  expect_s3_class(res, "pdactrace_user_projection")
  expect_named(res,
               c("rna_fit", "rna_pattern", "prot_fit", "prot_pattern",
                 "evidence", "audit", "summary"),
               ignore.order = TRUE)
  expect_s3_class(res$rna_fit,  "data.table")
  expect_s3_class(res$audit,    "data.table")
  expect_null(res$prot_fit)
  expect_false(res$summary$has_protein)
  expect_gt(res$summary$n_genes_audit, 0L)
})

test_that("project_user_cohort accepts a SummarizedExperiment", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = toy_counts),
    colData = S4Vectors::DataFrame(toy_coldata))
  res <- project_user_cohort(
    rna        = se,
    stage_col  = "stage",
    cohort_col = "cohort")
  expect_s3_class(res, "pdactrace_user_projection")
})

test_that("project_user_cohort layers in a protein cohort when supplied", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("limma")
  res <- project_user_cohort(
    rna        = toy_counts,
    coldata    = toy_coldata,
    stage_col  = "stage",
    cohort_col = "cohort",
    protein    = toy_protein)
  expect_true(res$summary$has_protein)
  expect_s3_class(res$prot_fit, "data.table")
})

test_that("project_user_cohort errors when coldata is missing for matrix input", {
  expect_error(
    project_user_cohort(rna = toy_counts, stage_col = "stage"),
    regexp = "coldata.*required")
})
