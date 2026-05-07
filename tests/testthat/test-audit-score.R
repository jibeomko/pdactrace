test_that("compute_audit_score returns 3-axis + 2-gate + 4-class output", {
  out <- compute_audit_score(c("LTBP1", "GAPDH", "ALB", "THBS2", "LGALS3BP"))
  expect_setequal(out$gene_symbol,
                  c("LTBP1", "GAPDH", "ALB", "THBS2", "LGALS3BP"))
  expect_true(all(c("evidence_strength", "biological_coherence",
                    "translational_relevance", "leakage_gate",
                    "heterogeneity_gate", "positive_score",
                    "audit_score", "audit_class") %in% names(out)))
  expect_equal(out[gene_symbol == "GAPDH", leakage_gate], 0)
  expect_equal(out[gene_symbol == "GAPDH", audit_score], 0)
  expect_equal(out[gene_symbol == "GAPDH", audit_class], "excluded")
  expect_equal(out[gene_symbol == "ALB", leakage_gate], 0.5)
  expect_equal(out[gene_symbol == "ALB", audit_class], "penalized")
  expect_equal(out[gene_symbol == "LTBP1", audit_class],
               "supported_uncertain")
  expect_equal(out[gene_symbol == "THBS2", audit_class],
               "high_confidence")
  expect_equal(out[gene_symbol == "LGALS3BP", audit_class],
               "high_confidence")
  expect_true(out[gene_symbol == "THBS2", audit_score] > 0.65)
})

test_that("propagate_uncertainty returns stored MC case-study classes", {
  out <- propagate_uncertainty(c("LTBP1", "GAPDH", "THBS2"))
  expect_equal(out[gene_symbol == "GAPDH", confidence_class], "excluded")
  expect_equal(out[gene_symbol == "GAPDH", audit_score_lo95], 0)
  expect_equal(out[gene_symbol == "GAPDH", audit_score_hi95], 0)
  expect_equal(out[gene_symbol == "THBS2", confidence_class], "stable_high")
  expect_equal(out[gene_symbol == "LTBP1", confidence_class],
               "high_uncertain")
  expect_true(out[gene_symbol == "LTBP1", uncertainty_width] > 0.5)
})

test_that("evaluate_anchor_enrichment reproduces secondary top-100 signal", {
  out <- evaluate_anchor_enrichment(top_n = 100, tier = "secondary")
  expect_true(out$hits >= 5L)
  expect_true(out$fold > 20)
  expect_true(out$pval < 1e-5)
})

test_that("build_evidence_graph and extract_graph_features are auditable", {
  g <- build_evidence_graph("LTBP1")
  expect_true(inherits(g, "pdactrace_evidence_graph"))
  expect_true(all(c("nodes", "edges", "score") %in% names(g)))
  expect_true(nrow(g$nodes) >= 8L)
  expect_true("rescue_signal" %in% g$nodes$node_id)
  expect_equal(g$nodes[node_id == "rescue_signal", value_numeric], 1)

  feat <- extract_graph_features("LTBP1")
  expect_true(feat$rescue_eligible)
  expect_true(feat$positive_score > 0.55)
})

test_that("case_study combines score, uncertainty, and call labels", {
  cs <- case_study(c("THBS2", "LTBP1", "GAPDH", "ALB"))
  expect_equal(nrow(cs), 4L)
  expect_equal(cs[gene_symbol == "THBS2", pdactrace_call],
               "confirmed_anchor")
  expect_equal(cs[gene_symbol == "LTBP1", pdactrace_call],
               "rescued_candidate")
  expect_equal(cs[gene_symbol == "GAPDH", pdactrace_call],
               "rejected_artifact")
  expect_equal(cs[gene_symbol == "ALB", pdactrace_call],
               "rejected_artifact")
})
