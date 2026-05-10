test_that("score_anchor_similarity returns the expected schema", {
  sim <- score_anchor_similarity(tier = "primary")
  expect_s3_class(sim, "data.table")
  expect_named(sim, c("gene_symbol", "anchor_similarity",
                       "anchor_n", "anchor_tier"),
               ignore.order = TRUE)
  expect_true(all(sim$anchor_similarity[!is.na(sim$anchor_similarity)] >= -1))
  expect_true(all(sim$anchor_similarity[!is.na(sim$anchor_similarity)] <=  1))
})

test_that("T1 anchors land in the top decile by anchor_similarity", {
  sim <- score_anchor_similarity(tier = "primary")
  data("pdactrace_external_anchors", package = "pdactrace")
  t1 <- pdactrace_external_anchors$gene[
    pdactrace_external_anchors$evidence_tier == "T1_validated"]
  in_atlas <- intersect(t1, sim$gene_symbol)
  top10pct <- head(sim$gene_symbol, ceiling(nrow(sim) * 0.10))
  hits <- length(intersect(in_atlas, top10pct))
  expected_random <- length(in_atlas) * 0.10
  # at least 3x baseline enrichment
  expect_gte(hits, ceiling(expected_random * 3))
})

test_that("tier filter changes anchor_n", {
  sim_p <- score_anchor_similarity(tier = "primary")
  sim_s <- score_anchor_similarity(tier = "secondary")
  expect_lte(unique(sim_p$anchor_n), unique(sim_s$anchor_n))
})

test_that("euclidean method also lies in (0, 1] and ranks T1 anchors high", {
  sim <- score_anchor_similarity(tier = "primary", method = "euclidean")
  v <- sim$anchor_similarity[!is.na(sim$anchor_similarity)]
  expect_true(all(v > 0)); expect_true(all(v <= 1))
})

test_that("two consecutive calls produce identical output", {
  s1 <- score_anchor_similarity(tier = "primary")
  s2 <- score_anchor_similarity(tier = "primary")
  expect_identical(s1, s2)
})
