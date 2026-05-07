test_that("fit_stage_de S4 dispatches to SummarizedExperiment method", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("DESeq2")

  # Synthetic 500-gene x 24-sample stage-progressive count matrix
  set.seed(42)
  n_g <- 500L
  n_s <- 24L
  stage <- rep(c("Normal", "Early", "Mid", "Late"), each = n_s / 4L)
  cohort <- rep(c("A", "B"), times = n_s / 2L)
  base <- matrix(stats::rpois(n_g * n_s, lambda = 50),
                  nrow = n_g, ncol = n_s)
  rownames(base) <- paste0("G", seq_len(n_g))
  # Inject a stage-progressive signal in the first 100 rows
  for (j in seq_len(n_s)) {
    s <- match(stage[j], c("Normal", "Early", "Mid", "Late")) - 1L
    base[seq_len(100), j] <- base[seq_len(100), j] + 5L * s
  }

  # ── Path 1: matrix interface ────────────────────────────────
  fit_mat <- fit_stage_de(base, stage, cohort)

  # ── Path 2: SummarizedExperiment interface ──────────────────
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = base),
    colData = S4Vectors::DataFrame(stage = stage, cohort = cohort))
  fit_se <- fit_stage_de(se, stage_col = "stage", cohort_col = "cohort")

  expect_s3_class(fit_mat, "data.table")
  expect_s3_class(fit_se,  "data.table")
  # The two paths must produce identical numeric output
  expect_equal(fit_mat$beta_E, fit_se$beta_E)
  expect_equal(fit_mat$lrt_padj, fit_se$lrt_padj)
})

test_that("fit_stage_de_protein S4 dispatches to SE method", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("limma")

  set.seed(42)
  n_g <- 200L
  n_s <- 24L
  stage <- rep(c("Normal", "Early", "Mid", "Late"), each = n_s / 4L)
  cohort <- rep(c("A", "B"), times = n_s / 2L)
  intensity <- matrix(stats::rnorm(n_g * n_s, mean = 12, sd = 1.5),
                       nrow = n_g, ncol = n_s)
  rownames(intensity) <- paste0("P", seq_len(n_g))
  for (j in seq_len(n_s)) {
    s <- match(stage[j], c("Normal", "Early", "Mid", "Late")) - 1L
    intensity[seq_len(40), j] <- intensity[seq_len(40), j] + 0.5 * s
  }

  fit_mat <- fit_stage_de_protein(intensity, stage, cohort)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(intensity = intensity),
    colData = S4Vectors::DataFrame(stage = stage, cohort = cohort))
  fit_se <- fit_stage_de_protein(
    se, stage_col = "stage", cohort_col = "cohort",
    assay_name = "intensity")

  expect_equal(fit_mat$beta_E, fit_se$beta_E)
  expect_equal(fit_mat$lrt_padj, fit_se$lrt_padj)
})

test_that("fit_stage_de errors on missing colData column", {
  skip_if_not_installed("SummarizedExperiment")

  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = matrix(1L, 5L, 5L,
                                     dimnames = list(LETTERS[1:5], NULL))),
    colData = S4Vectors::DataFrame(grade = letters[1:5]))
  expect_error(
    fit_stage_de(se, stage_col = "stage"),
    regexp = "stage")
  expect_error(
    fit_stage_de(se, stage_col = "grade", cohort_col = "missing"),
    regexp = "missing")
})
