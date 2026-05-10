test_that("viz_gene default returns pdactrace_viz_panels (split mode)", {
  p <- viz_gene("LGALS3BP")
  expect_s3_class(p, "pdactrace_viz_panels")
  expect_length(p, 6L)
  expect_setequal(names(p),
                   c("rna", "protein", "cell", "serum",
                     "filter", "hexagon"))
})

test_that("viz_gene(layout='compact') returns a patchwork object", {
  skip_if_not_installed("patchwork")
  p <- viz_gene("LGALS3BP", layout = "compact")
  expect_s3_class(p, "patchwork")
})

test_that("viz_gene(layout='compact') supports 2x3 / 3x2 / 6x1 layouts", {
  skip_if_not_installed("patchwork")
  p2 <- viz_gene("LGALS3BP", layout = "compact", ncol = 2L)
  p1 <- viz_gene("LGALS3BP", layout = "compact", ncol = 1L)
  expect_s3_class(p2, "patchwork")
  expect_s3_class(p1, "patchwork")
})

test_that("viz_gene errors clearly on missing gene", {
  expect_error(viz_gene("NOT_A_REAL_GENE_XYZ"),
               regexp = "not in the bundled atlas")
})

test_that("viz_gene(layout='compact') composes the six expected sub-panels", {
  skip_if_not_installed("patchwork")
  # patchwork stores assembled plots in $patches$plots; the final
  # plot lives on the patchwork object itself, hence + 1L.
  p <- viz_gene("LGALS3BP", layout = "compact")
  panels <- p$patches$plots
  expect_equal(length(panels) + 1L, 6L)
})

test_that("viz_gene gracefully degrades for genes with sparse layers", {
  for (g in c("ALB", "GAPDH", "SERPINA1")) {
    p <- viz_gene(g)
    expect_s3_class(p, "pdactrace_viz_panels")
    expect_length(p, 6L)
  }
})

test_that("plot_serum_direction returns a ggplot for genes with serum data", {
  p <- plot_serum_direction("LTBP1")
  expect_true(inherits(p, "ggplot"))
})

test_that("plot_serum_direction errors on empty input", {
  expect_error(plot_serum_direction(character(0L)),
               regexp = "non-empty character vector")
})

test_that("viz_gene(layout='split') returns named list of 6 ggplots", {
  panels <- viz_gene("LGALS3BP", layout = "split")
  expect_s3_class(panels, "pdactrace_viz_panels")
  expect_type(panels, "list")
  expect_length(panels, 6L)
  expect_setequal(names(panels),
                   c("rna", "protein", "cell", "serum",
                     "filter", "hexagon"))
  for (nm in names(panels)) {
    expect_true(inherits(panels[[nm]], "ggplot"),
                info = nm)
  }
})

test_that("print.pdactrace_viz_panels renders all six panels", {
  panels <- viz_gene("LGALS3BP")
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  out <- print(panels)
  expect_identical(out, panels)
})

test_that("viz_gene(layout='split', output_dir=...) writes 6 PDFs", {
  tmpd <- tempfile("viz_split_")
  res <- suppressMessages(viz_gene("LGALS3BP",
                                     layout = "split",
                                     output_dir = tmpd))
  files <- attr(res, "files")
  expect_equal(length(files), 6L)
  expect_true(all(file.exists(files)))
  expect_setequal(basename(files), sprintf(
    "viz_gene_LGALS3BP_%s.pdf",
    c("rna", "protein", "cell", "serum", "filter", "hexagon")))
})

test_that("viz_gene(layout='split') works for sparse-evidence genes", {
  panels <- viz_gene("ALB", layout = "split")
  expect_length(panels, 6L)
  for (nm in names(panels))
    expect_true(inherits(panels[[nm]], "ggplot"), info = nm)
})
