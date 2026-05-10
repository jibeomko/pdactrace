# Helpers --------------------------------------------------------

.make_synth_serum <- function(seed = 1, g = 50, n = 24) {
  set.seed(seed)
  X <- matrix(stats::rnorm(g * n, mean = 5), nrow = g, ncol = n,
              dimnames = list(paste0("GENE", seq_len(g)),
                               paste0("S", seq_len(n))))
  X[1:10, 1:8] <- X[1:10, 1:8] + 2L
  cd <- data.frame(
    sample = colnames(X),
    group  = rep(c("PDAC", "HC", "Pancreatitis"), each = 8),
    stringsAsFactors = FALSE)
  list(intensity = X, coldata = cd)
}

# Tests ----------------------------------------------------------

test_that("project_user_serum_cohort returns the expected schema", {
  d <- .make_synth_serum()
  ss <- suppressWarnings(project_user_serum_cohort(
    d$intensity, d$coldata, pan_label = "Pancreatitis"))
  expect_s3_class(ss, "data.table")
  expect_named(ss, c("gene_symbol",
                      "serum_log2fc_PDAC_vs_HC",
                      "serum_padj_PDAC_vs_HC",
                      "serum_log2fc_Pan_vs_HC",
                      "serum_padj_Pan_vs_HC",
                      "translation_class",
                      "serum_detected"),
               ignore.order = TRUE)
})

test_that("project_user_serum_cohort log2FC matches manual mean diff", {
  d <- .make_synth_serum()
  ss <- suppressWarnings(project_user_serum_cohort(
    d$intensity, d$coldata, pan_label = "Pancreatitis"))
  manual <- unname(rowMeans(d$intensity[, 1:8]) -
                    rowMeans(d$intensity[, 9:16]))
  expect_equal(ss$serum_log2fc_PDAC_vs_HC, manual,
               tolerance = 1e-10)
})

test_that("planted PDAC-up signal becomes serum_detected", {
  d <- .make_synth_serum()
  ss <- suppressWarnings(project_user_serum_cohort(
    d$intensity, d$coldata, pan_label = "Pancreatitis",
    test = "wilcox"))
  hits <- ss[serum_detected == TRUE, gene_symbol]
  expect_true(any(c("GENE1","GENE5","GENE10") %in% hits))
})

test_that("missing PDAC / HC labels error with helpful message", {
  d <- .make_synth_serum()
  d$coldata$group <- "Cancer"   # no matches
  expect_error(project_user_serum_cohort(d$intensity, d$coldata),
               regexp = "must contain group labels")
})

test_that("nrow(coldata) != ncol(intensity) errors", {
  d <- .make_synth_serum()
  expect_error(project_user_serum_cohort(d$intensity,
                                          d$coldata[1:3, ]),
               regexp = "nrow\\(coldata\\)")
})

test_that("missing group_col errors with helpful message", {
  d <- .make_synth_serum()
  d$coldata$group <- NULL
  expect_error(project_user_serum_cohort(d$intensity, d$coldata),
               regexp = "no column 'group'")
})

test_that("link_to_atlas = FALSE leaves translation_class as NA", {
  d <- .make_synth_serum()
  ss <- suppressWarnings(project_user_serum_cohort(
    d$intensity, d$coldata, pan_label = "Pancreatitis",
    link_to_atlas = FALSE))
  expect_true(all(is.na(ss$translation_class)))
})

test_that("translation_class joined from atlas matches sign comparison", {
  d <- .make_synth_serum()
  rownames(d$intensity)[1:6] <- c("LTBP1", "LGALS3BP", "TIMP1",
                                    "ALB", "GAPDH", "SERPINA1")
  ss <- suppressWarnings(project_user_serum_cohort(
    d$intensity, d$coldata, pan_label = "Pancreatitis"))
  hit <- ss[gene_symbol %in% c("LTBP1","LGALS3BP","TIMP1","ALB",
                                 "GAPDH","SERPINA1")]
  # All 6 had positive PDAC-vs-HC log2FC in the synth signal:
  # tissue Early_Burst_Up genes (LTBP1, LGALS3BP, SERPINA1) -> A;
  # Early_Trough (ALB) -> B; NA pattern (GAPDH) -> C.
  expect_equal(hit[gene_symbol == "LTBP1",   translation_class], "A")
  expect_equal(hit[gene_symbol == "LGALS3BP",translation_class], "A")
  expect_equal(hit[gene_symbol == "ALB",     translation_class], "B")
  expect_equal(hit[gene_symbol == "GAPDH",   translation_class], "C")
})

test_that("fit_stage_de friendly error fires on N vs T input", {
  data("toy_counts",  package = "pdactrace")
  data("toy_coldata", package = "pdactrace")
  binary_stage <- ifelse(toy_coldata$stage == "Normal",
                         "Normal", "Tumor")
  expect_error(
    fit_stage_de(toy_counts, stage = binary_stage),
    regexp = "stage-aware data")
})
