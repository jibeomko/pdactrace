test_that("compare_genes long form returns a long data.table with expected columns", {
  tab <- compare_genes(c("LGALS3BP", "LTBP1", "TIMP1"))
  expect_s3_class(tab, "data.table")
  expect_named(tab, c("gene", "axis", "metric", "value"),
               ignore.order = TRUE)
  # All three genes should be represented, with > 1 row each
  expect_setequal(unique(tab$gene), c("LGALS3BP", "LTBP1", "TIMP1"))
  expect_true(all(table(tab$gene) > 10L))
})

test_that("compare_genes wide form has one row per gene with axis.metric columns", {
  w <- compare_genes(c("LGALS3BP", "LTBP1"), wide = TRUE)
  expect_s3_class(w, "data.table")
  expect_equal(nrow(w), 2L)
  expect_true("trajectory_fit.delta_rho" %in% names(w))
  expect_true("rna_protein_coupling.cosine" %in% names(w))
  expect_true("filter_survival.passed" %in% names(w))
})

test_that("compare_genes axes filter restricts the metric set", {
  tab <- compare_genes(c("LGALS3BP", "LTBP1"),
                        axes = "trajectory_fit")
  expect_setequal(unique(tab$axis), "trajectory_fit")
  # rna_pattern + rho_best + rho_runner_up + delta_rho + note = 5 metrics
  expect_equal(nrow(tab), 2L * 5L)
})

test_that("compare_genes flags missing genes inline rather than erroring", {
  tab <- compare_genes(c("LGALS3BP", "NOT_A_REAL_GENE_XYZ"))
  expect_s3_class(tab, "data.table")
  miss <- tab[gene == "NOT_A_REAL_GENE_XYZ"]
  expect_equal(nrow(miss), 1L)
  expect_equal(miss$metric, "missing")
  expect_match(miss$value, "not in bundled atlas")
})

test_that("compare_genes errors on unknown axis name", {
  expect_error(compare_genes("LGALS3BP", axes = "not_an_axis"),
               regexp = "Unknown axis name")
})
