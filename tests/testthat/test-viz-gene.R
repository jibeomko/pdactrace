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

test_that("viz_gene composes the four expected sub-panels", {
  skip_if_not_installed("patchwork")
  # patchwork stores assembled plots in $patches$plots
  p <- viz_gene("LGALS3BP")
  panels <- p$patches$plots
  expect_equal(length(panels) + 1L, 4L)  # last plot is in p itself
})
