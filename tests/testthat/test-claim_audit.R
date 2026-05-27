test_that("classify_claim_tier returns expected bundled-gene tiers", {
  out <- classify_claim_tier(genes = c("LGALS3BP", "LTBP1", "ALB",
                                       "GAPDH"))
  expect_named(out,
               c("gene_symbol", "claim_tier", "claim_strength",
                 "tissue_supported", "exportable_plausible",
                 "serum_observed", "serum_signal",
                 "serum_concordant", "translation_status",
                 "confounder_risk", "tracd_class", "tracd_confidence",
                 "audit_class", "audit_score", "claim_reason"))
  expect_equal(out[gene_symbol == "GAPDH", claim_tier], "excluded")
  expect_equal(out[gene_symbol == "ALB", claim_tier], "confounded")
  expect_true(out[gene_symbol == "LGALS3BP", claim_strength] >=
                out[gene_symbol == "LTBP1", claim_strength])
})

test_that("classify_claim_tier separates tissue support from blood claim", {
  toy <- data.table::data.table(
    gene_symbol = c("A_serum", "B_export", "C_unknown", "D_hk",
                    "E_plasma", "F_none"),
    rna_pattern = c("Early_Burst_Up", "Early_Burst_Up",
                    "Early_Burst_Up", "Early_Burst_Up",
                    "Early_Burst_Up", NA_character_),
    prot_pattern = c("Early_Burst_Up", "Early_Burst_Up",
                     "Early_Burst_Up", "Early_Burst_Up",
                     "Early_Burst_Up", NA_character_),
    cell_origin_top = c("ductal", "ductal", "ductal", "ductal",
                        "ductal", NA_character_),
    max_abs_beta_meta = c(1, 1, 1, 1, 1, 0),
    max_I2_meta = c(10, 10, 10, 10, 10, NA),
    rna_pattern_rho = c(0.95, 0.95, 0.95, 0.95, 0.95, NA),
    rna_lrt_padj = c(1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1),
    rna_cohort_agreement = c(1, 1, 1, 1, 1, 0),
    serum_detected = c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE),
    serum_n_cohorts_detected = c(1L, 0L, 0L, 1L, 1L, 0L),
    serum_log2fc_PDAC_vs_HC = c(0.8, NA, NA, 0.8, 0.8, NA),
    serum_log2fc_Pan_vs_HC = c(0, NA, NA, 0, 0, NA),
    translation_class = c("A", NA_character_, NA_character_, "B",
                          "A", NA_character_),
    flt_signal_peptide = c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE),
    flt_direction_match = c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE),
    audit_score_layer = c(1, 0.75, 0.5, 1, 1, 0),
    audit_score_direction = c(1, 1, 1, 1, 1, 0),
    audit_score_early = c(1, 1, 1, 1, 1, 0),
    audit_score_serum = c(1, 0.3, 0, 1, 1, 0),
    audit_score_rescue = c(0, 0, 0, 0, 0, 0),
    audit_positive_score = c(1, 0.75, 0.5, 1, 1, 0),
    audit_leakage_mult = c(1, 1, 1, 0, 0.5, 1),
    audit_heterogeneity_mult = c(1, 1, 1, 1, 1, 1),
    audit_score_raw = c(1, 0.75, 0.5, 0, 0.5, 0),
    audit_score = c(1, 0.75, 0.5, 0, 0.5, 0),
    audit_is_housekeeping = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
    audit_is_plasma_high_abundance = c(FALSE, FALSE, FALSE, FALSE,
                                       TRUE, FALSE),
    audit_rescue_eligible = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
  )
  out <- classify_claim_tier(atlas = toy)
  expect_equal(out$claim_tier,
               c("serum_concordant", "exportable_plausible",
                 "translation_unknown", "excluded", "confounded",
                 "insufficient_tissue_evidence"))
})

test_that("summarize_translation_gap returns tier counts", {
  out <- summarize_translation_gap()
  expect_true(all(c("claim_tier", "n", "pct") %in% names(out)))
  expect_equal(sum(out$n), nrow(pdactrace:::.get_reference()))
  expect_equal(sum(out$pct), 1, tolerance = 1e-8)
})

test_that("run_weight_robustness returns top-N stability proportions", {
  out <- run_weight_robustness(genes = c("LGALS3BP", "LTBP1"),
                               n_draws = 20L,
                               top_n = c(50L, 100L),
                               seed = 1L)
  expect_equal(nrow(out), 2L)
  expect_true(all(c("weight_rank_median", "weight_top50_stability",
                    "weight_top100_stability") %in% names(out)))
  expect_true(all(out$weight_top50_stability >= 0 &
                    out$weight_top50_stability <= 1))
})

test_that("compute_pareto_class reports manuscript-facing classes", {
  out <- compute_pareto_class(genes = c("LGALS3BP", "LTBP1", "GAPDH"),
                              top_n = nrow(pdactrace:::.get_reference()))
  expect_true(all(c("gene_symbol", "pareto_class", "pareto_layer",
                    "pareto_rank") %in% names(out)))
  expect_true(all(out$pareto_class %in%
                    c("pareto_front", "pareto_supported",
                      "pareto_lower", "gate_excluded",
                      "outside_top_n")))
  expect_equal(out[gene_symbol == "GAPDH", pareto_class],
               "gate_excluded")
})
