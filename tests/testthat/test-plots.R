test_that("plot_gene_evidence returns a patchwork", {
  p <- plot_gene_evidence("LTBP1")
  expect_true(inherits(p, "patchwork"))
})

test_that("plot_gene_evidence returns NULL for unknown gene", {
  expect_message(p <- plot_gene_evidence("FOOBAR_NOT_A_GENE"),
                  "No evidence")
  expect_null(p)
})

test_that("plot_filter_trace returns patchwork by default", {
  p <- plot_filter_trace(c("SERPINA1", "LTBP1", "SPARC"))
  expect_true(inherits(p, "patchwork"))
})

test_that("plot_filter_trace returns ggplot when show_routes=FALSE", {
  p <- plot_filter_trace("LTBP1", show_routes = FALSE)
  expect_true(inherits(p, "ggplot"))
})

test_that("plot_panel_heatmap returns ggplot", {
  p <- plot_panel_heatmap(c("LTBP1", "SERPINA1", "CDH13"))
  expect_true(inherits(p, "ggplot"))
})

test_that("plot_candidate_landscape returns ggplot", {
  p <- plot_candidate_landscape()
  expect_true(inherits(p, "ggplot"))
})
