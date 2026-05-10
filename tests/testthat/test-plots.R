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

test_that("plot_filter_trace serum strip can be toggled off", {
  p1 <- plot_filter_trace("LTBP1", show_routes = FALSE,
                           show_serum = FALSE)
  p2 <- plot_filter_trace("LTBP1", show_routes = FALSE,
                           show_serum = TRUE)
  # Both render; with show_serum=TRUE the result is a patchwork
  # (ggplot subclass), without it just a ggplot.
  expect_true(inherits(p1, "ggplot"))
  expect_true(inherits(p2, "ggplot"))
})

test_that("plot_stage_effect supports layer = 'protein'", {
  p_rna  <- plot_stage_effect("LTBP1", layer = "rna")
  p_prot <- plot_stage_effect("LTBP1", layer = "protein")
  expect_true(inherits(p_rna, "ggplot"))
  expect_true(inherits(p_prot, "ggplot"))
})

test_that("plot_stage_effect protein gracefully handles genes absent from protein atlas", {
  # Pick a gene we know is in pdactrace_reference but NOT in
  # pdactrace_protein_betas.
  data("pdactrace_reference",     package = "pdactrace")
  data("pdactrace_protein_betas", package = "pdactrace")
  rna_only <- setdiff(pdactrace_reference$gene_symbol,
                       pdactrace_protein_betas$gene_symbol)
  if (length(rna_only) == 0L) skip("no RNA-only gene in atlas")
  p <- plot_stage_effect(rna_only[1L], layer = "protein")
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
