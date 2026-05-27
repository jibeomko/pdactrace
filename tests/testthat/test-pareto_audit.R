test_that(".fast_nondominated_sort_3d handles 5-point hierarchy", {
  # Hand-derived layers: gene 1 strictly dominates the rest;
  # genes 2 and 3 are mutually non-dominating after 1 is removed
  # (gene 2 strong on axes 2,3; gene 3 strong on axis 1);
  # gene 4 dominated by 2 and 3; gene 5 dominated by 4.
  X <- matrix(c(
    1.0, 1.0, 1.0,   # 1: layer 1 (strict apex)
    0.8, 0.9, 0.9,   # 2: dominated by 1 -> layer 2
    0.9, 0.8, 0.5,   # 3: dominated by 1, non-dominated by 2 -> layer 2
    0.7, 0.7, 0.4,   # 4: dominated by 2 and 3 -> layer 3
    0.5, 0.5, 0.3    # 5: dominated by 4 -> layer 4
  ), ncol = 3L, byrow = TRUE)
  layers <- pdactrace:::.fast_nondominated_sort_3d(X)
  expect_equal(layers, c(1L, 2L, 2L, 3L, 4L))
})

test_that(".fast_nondominated_sort_3d finds non-trivial layer 1", {
  # Three mutually non-dominating high-corner points + a clearly
  # dominated interior. Gene 4 is dominated by all three corners
  # because each corner is >= on every axis and strictly > on one.
  X <- matrix(c(
    1.0, 0.5, 0.5,
    0.5, 1.0, 0.5,
    0.5, 0.5, 1.0,
    0.4, 0.4, 0.4
  ), ncol = 3L, byrow = TRUE)
  layers <- pdactrace:::.fast_nondominated_sort_3d(X)
  expect_equal(layers, c(1L, 1L, 1L, 2L))
})

test_that(".fast_nondominated_sort_3d handles trivial inputs", {
  expect_equal(pdactrace:::.fast_nondominated_sort_3d(matrix(numeric(0),
                                                              ncol = 3L)),
               integer(0))
  expect_equal(pdactrace:::.fast_nondominated_sort_3d(matrix(c(1, 1, 1),
                                                              ncol = 3L)),
               1L)
  # Two ties on every axis -> both in layer 1 (neither dominates).
  X <- matrix(c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5), ncol = 3L, byrow = TRUE)
  expect_equal(pdactrace:::.fast_nondominated_sort_3d(X), c(1L, 1L))
})

test_that(".crowding_distance_3d gives Inf to boundary points", {
  X <- matrix(c(
    1.0, 0.5, 0.5,   # boundary on axis 1 (high)
    0.5, 1.0, 0.5,   # boundary on axis 2 (high)
    0.5, 0.5, 1.0,   # boundary on axis 3 (high)
    0.7, 0.7, 0.7    # interior
  ), ncol = 3L, byrow = TRUE)
  layers <- pdactrace:::.fast_nondominated_sort_3d(X)
  expect_true(all(layers == 1L))
  cd <- pdactrace:::.crowding_distance_3d(X, layers)
  expect_true(all(is.infinite(cd[1:3])))
  expect_true(is.finite(cd[4]))
})

test_that("compute_pareto_layers preserves frozen audit_score", {
  skip_on_cran()
  out <- compute_pareto_layers(top_n = 500L)
  expect_true(all(c("gene_symbol", "pareto_layer", "pareto_rank",
                    "pareto_excluded_by_gate") %in% names(out)))
  # At least one gene in layer 1 of a non-empty top-N pool.
  expect_gte(sum(out$pareto_layer == 1L, na.rm = TRUE), 1L)
  # Excluded genes get NA layer + NA rank.
  excluded <- out[out$pareto_excluded_by_gate]
  expect_true(all(is.na(excluded$pareto_layer)))
  expect_true(all(is.na(excluded$pareto_rank)))
  # All in-pool ranks are unique 1..N.
  in_pool <- out[!out$pareto_excluded_by_gate]
  expect_equal(sort(in_pool$pareto_rank), seq_len(nrow(in_pool)))
})

test_that("compute_pareto_layers rejects malformed atlas", {
  bad <- data.table::data.table(gene_symbol = "X")
  expect_error(compute_pareto_layers(atlas = bad),
               "missing required columns")
})

test_that("evaluate_pareto_stability returns expected schema", {
  skip_on_cran()
  out <- evaluate_pareto_stability(n_draws = 20L, top_n = 200L,
                                    seed = 1L)
  expect_named(out,
               c("gene_symbol", "pareto_layer_median",
                 "pareto_stability_top1", "pareto_layer_lo95",
                 "pareto_layer_hi95", "pareto_top10_pct_stability"))
  # Stability is a proportion in [0, 1].
  expect_true(all(out$pareto_stability_top1 >= 0 &
                    out$pareto_stability_top1 <= 1, na.rm = TRUE))
  expect_true(all(out$pareto_top10_pct_stability >= 0 &
                    out$pareto_top10_pct_stability <= 1, na.rm = TRUE))
})

test_that("evaluate_anchor_enrichment accepts pareto_rank score_col", {
  skip_on_cran()
  ref <- pdactrace:::.get_reference(NULL)
  skip_if_not("pareto_rank" %in% names(ref),
              "atlas does not carry pareto_rank yet")
  res <- evaluate_anchor_enrichment(top_n = 100, tier = "secondary",
                                     score_col = "pareto_rank")
  expect_true(nrow(res) == 1L)
  expect_true(is.numeric(res$fold))
  expect_true(is.numeric(res$pval))
  expect_true(res$pval >= 0 && res$pval <= 1)
})
