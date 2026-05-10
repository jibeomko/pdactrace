test_that("viz_gene returns a patchwork object for a real gene", {
  skip_if_not_installed("patchwork")
  p <- viz_gene("LGALS3BP")
  expect_s3_class(p, "patchwork")
})

test_that("viz_gene supports 2x2 (default) and 4x1 (vertical) layouts", {
  skip_if_not_installed("patchwork")
  p2 <- viz_gene("LGALS3BP", ncol = 2L)
  p1 <- viz_gene("LGALS3BP", ncol = 1L)
  expect_s3_class(p2, "patchwork")
  expect_s3_class(p1, "patchwork")
})

test_that("viz_gene errors clearly on missing gene", {
  skip_if_not_installed("patchwork")
  expect_error(viz_gene("NOT_A_REAL_GENE_XYZ"),
               regexp = "not in the bundled atlas")
})

test_that("viz_gene composes the six expected sub-panels", {
  skip_if_not_installed("patchwork")
  # patchwork stores assembled plots in $patches$plots; the final
  # plot lives on the patchwork object itself, hence + 1L.
  p <- viz_gene("LGALS3BP")
  panels <- p$patches$plots
  expect_equal(length(panels) + 1L, 6L)
})

test_that("viz_gene gracefully degrades for genes with sparse layers", {
  skip_if_not_installed("patchwork")
  for (g in c("ALB", "GAPDH", "SERPINA1")) {
    p <- viz_gene(g)
    expect_s3_class(p, "patchwork")
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
