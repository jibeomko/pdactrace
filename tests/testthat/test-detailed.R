test_that("query_gene_detailed returns expected slots for LTBP1", {
  d <- query_gene_detailed("LTBP1")
  expect_true(inherits(d, "pdactrace_gene_detailed"))
  expect_setequal(setdiff(names(d), c("summary", "provenance")),
                    c("per_stage", "per_cohort", "per_celltype",
                      "filter_diag", "serum_per_cohort", "meta"))
})

test_that("per_stage table has 4 stages with valid CI math", {
  d <- query_gene_detailed("LTBP1")
  ps <- d$per_stage
  expect_equal(nrow(ps), 4L)
  expect_setequal(as.character(ps$stage),
                    c("Normal", "Early", "Mid", "Late"))
  # CI math: ci_lo = beta - 1.96*SE, ci_hi = beta + 1.96*SE
  for (i in seq_len(nrow(ps))) {
    if (!is.na(ps$lfcSE[i]) && ps$lfcSE[i] > 0) {
      expect_equal(ps$ci_lo[i], ps$log2FC[i] - 1.96 * ps$lfcSE[i],
                    tolerance = 1e-9)
      expect_equal(ps$ci_hi[i], ps$log2FC[i] + 1.96 * ps$lfcSE[i],
                    tolerance = 1e-9)
    }
  }
})

test_that("per_cohort exposes 4 RNA cohorts with attributes", {
  d <- query_gene_detailed("LTBP1")
  pc <- d$per_cohort
  expect_setequal(pc$cohort,
                    c("TCGA", "CPTAC", "GSE224564", "GSE79668"))
  # Stouffer attributes are exposed
  expect_true(!is.null(attr(pc, "stouffer_p")))
  expect_true(!is.null(attr(pc, "agreement_pct")))
})

test_that("per_celltype is sorted descending and includes tau", {
  d <- query_gene_detailed("LTBP1")
  ct <- d$per_celltype
  expect_true(nrow(ct) >= 5L)
  expect_true(all(diff(ct$mean_expression) <= 0))   # descending
  expect_true(!is.null(attr(ct, "specificity_tau")))
  # LTBP1 myCAF should be top
  expect_equal(ct$celltype[1], "myCAF")
})

test_that("filter_diag has 7 steps with underlying metrics", {
  d <- query_gene_detailed("LTBP1")
  fd <- d$filter_diag
  expect_equal(nrow(fd), 7L)
  expect_setequal(fd$step,
                    c("signal_peptide", "serum_measurable",
                      "serum_significant", "pancreatitis_pdac",
                      "pancreatitis_hc", "direction_match", "final"))
  expect_true(all(grepl(".+", fd$underlying_metric)))
})

test_that("SERPINA1 passes all 7 filters (phase60_final)", {
  d <- query_gene_detailed("SERPINA1")
  fd <- d$filter_diag
  expect_true(all(fd$pass == TRUE))
})

test_that("LTBP1 fails 6 of 7 phase60 filters (only SignalP passes)", {
  d <- query_gene_detailed("LTBP1")
  fd <- d$filter_diag
  expect_equal(sum(fd$pass == TRUE, na.rm = TRUE), 1L)
  expect_true(fd[step == "signal_peptide", pass])
  expect_false(fd[step == "final", pass])
})

test_that("summarize_gene_evidence(detail=TRUE) adds per-stage block", {
  default_out <- summarize_gene_evidence("LTBP1")
  detail_out  <- summarize_gene_evidence("LTBP1", detail = TRUE)
  expect_true(nchar(detail_out) > nchar(default_out))
  expect_match(detail_out, "Per-stage:")
  expect_match(detail_out, "Per-cohort:")
})

test_that("4 detailed plot functions return ggplot objects", {
  expect_true(inherits(plot_stage_effect("LTBP1"), "ggplot"))
  expect_true(inherits(plot_per_cohort("LTBP1"), "ggplot"))
  expect_true(inherits(plot_filter_diagnostics("LTBP1"), "ggplot"))
  expect_true(inherits(plot_celltype_full("LTBP1"), "ggplot"))
})
