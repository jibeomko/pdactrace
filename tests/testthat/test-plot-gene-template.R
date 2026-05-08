test_that("plot_gene_template returns ggplot for an Early-surface gene", {
  p <- plot_gene_template("LGALS3BP", "rna")
  expect_s3_class(p, "ggplot")
  # Title contains gene name and rho annotation
  expect_match(p$labels$title, "LGALS3BP")
  expect_match(p$labels$title, "rho")
})

test_that("plot_gene_template falls back to non-Early best-match", {
  # GAPDH: housekeeping — atlas rna_pattern is NA but argmax recovers
  # a non-Early template. Confirm the subtitle flag is set.
  p <- plot_gene_template("GAPDH", "rna")
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$subtitle %||% "",
               "non-Early best-match", fixed = TRUE)
})

test_that("plot_gene_template writes PDF when output_file is set", {
  tmp <- tempfile("gene-tpl-", fileext = ".pdf")
  fp <- plot_gene_template("LTBP1", "rna", output_file = tmp)
  expect_true(file.exists(fp))
  expect_match(fp, "\\.pdf$")
  unlink(fp)
})

test_that("plot_gene_template works on the protein layer", {
  p <- plot_gene_template("LGALS3BP", "protein")
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$title, "LGALS3BP")
})

test_that("plot_gene_template errors on missing gene", {
  expect_error(plot_gene_template("NOT_A_REAL_GENE_XYZ", "rna"),
               regexp = "not in")
})

# tiny null-coalesce helper for the subtitle expectation above (NULL
# subtitle when not set)
`%||%` <- function(a, b) if (is.null(a)) b else a
