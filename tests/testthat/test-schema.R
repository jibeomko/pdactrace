test_that("schema_spec returns 113 columns aligned with reference DB (v0.4.0)", {
  s <- schema_spec()
  expect_equal(nrow(s), 113L)
  expect_setequal(unique(s$layer),
                    c("identifier", "rna", "protein", "scrna",
                      "serum", "clinical", "filter_status",
                      "annotation", "provenance",
                      "meta_analysis", "tier", "audit", "t3_ready"))

  ref <- pdactrace:::.get_reference()
  expect_setequal(s$name, names(ref))
})

test_that("Random-effects meta-analysis columns populated", {
  ref <- pdactrace:::.get_reference()
  for (col in c("meta_NvE_beta", "meta_MvE_I2", "meta_LvE_padj",
                "max_abs_beta_meta", "max_I2_meta",
                "meta_cohort_divergent"))
    expect_true(col %in% names(ref))
  # v0.4.0: confidence_tier / early_onset_score / heterogeneity_factor
  # were removed because v0.3.0 audit_class supersedes them.
  for (col in c("confidence_tier", "early_onset_score",
                "heterogeneity_factor"))
    expect_false(col %in% names(ref))
})

test_that("v0.1.1 detail columns are populated", {
  ref <- pdactrace:::.get_reference()
  expect_true("rna_padj_E" %in% names(ref))
  expect_true("rna_lfcSE_E" %in% names(ref))
  expect_true(any(!is.na(ref$rna_lfcSE_E)))
  # v0.4.0: per-contrast Wald padj distinct from LRT padj copies
  for (col in c("rna_wald_padj_E", "rna_wald_padj_M",
                "rna_wald_padj_L"))
    expect_true(col %in% names(ref))
  expect_true(any(!is.na(ref$rna_wald_padj_E)))
  expect_true("rna_per_cohort_trend" %in% names(ref))
  expect_true(is.list(ref$rna_per_cohort_trend))
  expect_true("rna_stouffer_p" %in% names(ref))
  for (col in c("ann_pan_vs_hc_logfc", "ann_pan_vs_hc_pval",
                "ann_pdac_mean", "ann_pan_mean", "ann_hc_mean",
                "ann_pdac_vs_pan_pval", "cell_specificity_tau"))
    expect_true(col %in% names(ref))
})

test_that("4-Early atlas universe is 1,356 (RNA × Protein concordant, v0.4.0)", {
  ref <- pdactrace:::.get_reference()
  n_4early <- ref[rnaprot_concordant == TRUE &
                    rna_pattern %in% early_pattern_names(),
                    data.table::uniqueN(gene_symbol)]
  expect_equal(n_4early, 1356L)
})

test_that("12-template exclusion flags (v0.4.0): Mid/Late/Monotonic", {
  ref <- pdactrace:::.get_reference()
  for (col in c("excluded_mid_pattern", "excluded_late_pattern",
                "excluded_monotonic_pattern"))
    expect_true(col %in% names(ref))
  expect_true(sum(ref$excluded_mid_pattern,        na.rm = TRUE) > 0)
  expect_true(sum(ref$excluded_late_pattern,       na.rm = TRUE) > 0)
  expect_true(sum(ref$excluded_monotonic_pattern,  na.rm = TRUE) > 0)
  # Any excluded pattern call -> rna_pattern is NA (Early × 4 surface only)
  excluded_any <- ref$excluded_mid_pattern |
                   ref$excluded_late_pattern |
                   ref$excluded_monotonic_pattern
  expect_true(all(is.na(ref$rna_pattern[excluded_any])))
})

test_that("phase60 7-step coverage = 1,449 / final pass count > 0", {
  ref <- pdactrace:::.get_reference()
  n_p60 <- ref[!is.na(flt_signal_peptide), .N]
  expect_equal(n_p60, 1449L)
  n_final <- ref[isTRUE(flt_final) | (!is.na(flt_final) & flt_final == TRUE), .N]
  expect_true(n_final > 0)
  expect_true(n_final < n_p60)  # strict funnel
})

test_that("LTBP1 = phase77_classB, SERPINA1 = phase60_final routes", {
  tf <- trace_filters(c("LTBP1", "SERPINA1"))
  expect_equal(tf[gene_symbol == "LTBP1", class_route],
                "phase77_classB")
  expect_equal(tf[gene_symbol == "SERPINA1", class_route],
                "phase60_final")
  expect_equal(tf[gene_symbol == "SERPINA1", n_phase60_pass], 7L)
  expect_equal(tf[gene_symbol == "LTBP1", n_phase60_pass], 1L)
})

test_that("Direction T3 columns are NA placeholders", {
  ref <- pdactrace:::.get_reference()
  expect_true(all(is.na(ref$serum_direction_label)))
  expect_true(all(is.na(ref$direction_model_trainable)))
  expect_true(all(is.na(ref$direction_model_card_ref)))
})

test_that("v0.3.0 audit columns (3 axes + 2 gates + 4 classes) are populated", {
  ref <- pdactrace:::.get_reference()
  for (col in c("audit_score", "audit_class",
                "audit_evidence_strength", "audit_biological_coherence",
                "audit_translational_relevance",
                "audit_leakage_gate", "audit_heterogeneity_gate",
                "audit_score_median", "audit_score_lo95",
                "audit_score_hi95", "audit_confidence_class"))
    expect_true(col %in% names(ref))
  expect_equal(ref[gene_symbol == "GAPDH", audit_score], 0)
  expect_equal(ref[gene_symbol == "GAPDH", audit_class], "excluded")
  expect_equal(ref[gene_symbol == "GAPDH", audit_confidence_class], "excluded")
  expect_equal(ref[gene_symbol == "ALB", audit_class], "penalized")
  expect_equal(ref[gene_symbol == "LTBP1", audit_class], "supported_uncertain")
  expect_true(ref[gene_symbol == "LTBP1", audit_score_median] > 0.5)
  expect_equal(ref[gene_symbol == "THBS2", audit_class], "high_confidence")
  expect_equal(ref[gene_symbol == "THBS2", audit_confidence_class],
               "stable_high")
  expect_equal(ref[gene_symbol == "LGALS3BP", audit_class], "high_confidence")
})

test_that("4-Early templates surface in atlas correctly", {
  expect_setequal(early_pattern_names(),
    c("Early_Burst_Up", "Early_Loss_Down",
      "Early_Peak", "Early_Trough"))
  expect_setequal(mid_pattern_names_excluded(),
    c("Mid_Peak", "Mid_Trough",
      "Mid_Plateau_Up", "Mid_Plateau_Down"))
})
