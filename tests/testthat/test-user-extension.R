test_that("compute_audit_score(evidence = NULL) matches bundled atlas", {
  default_out <- compute_audit_score(c("LGALS3BP", "LTBP1", "ALB", "GAPDH"))
  null_out <- compute_audit_score(c("LGALS3BP", "LTBP1", "ALB", "GAPDH"),
                                    evidence = NULL)
  expect_equal(default_out, null_out)
})

test_that("assemble_user_evidence builds a scoreable evidence table", {
  rna <- data.table::data.table(
    gene_symbol = c("FAKE_UP_1", "FAKE_DOWN_1", "FAKE_FLAT"),
    rna_pattern = c("Early_Burst_Up", "Early_Loss_Down", NA_character_),
    rna_pattern_rho = c(0.95, 0.92, NA_real_),
    lrt_padj = c(1e-4, 1e-3, NA_real_))

  prot <- data.table::data.table(
    gene_symbol = c("FAKE_UP_1", "FAKE_DOWN_1"),
    prot_pattern = c("Early_Burst_Up", "Early_Loss_Down"))

  ev <- assemble_user_evidence(rna_fit = rna, prot_fit = prot,
                                signal_peptide = c("FAKE_UP_1"))
  expect_true("rna_pattern" %in% names(ev))
  expect_true("prot_pattern" %in% names(ev))
  expect_true("flt_signal_peptide" %in% names(ev))
  expect_true("max_I2_meta" %in% names(ev))
  expect_equal(nrow(ev), 3L)

  scored <- compute_audit_score(evidence = ev)
  expect_true(all(c("evidence_strength", "biological_coherence",
                    "translational_relevance", "leakage_gate",
                    "heterogeneity_gate", "audit_score",
                    "audit_class") %in% names(scored)))
  expect_equal(nrow(scored), 3L)
})

test_that("evaluate_anchor_enrichment accepts user evidence + anchors", {
  rna <- data.table::data.table(
    gene_symbol = paste0("GENE", 1:50),
    rna_pattern = c(rep("Early_Burst_Up", 10),
                    rep("Early_Loss_Down", 10),
                    rep(NA_character_, 30)),
    rna_pattern_rho = c(rep(0.95, 20), rep(NA_real_, 30)),
    lrt_padj = c(rep(1e-4, 20), rep(NA_real_, 30)))
  prot <- data.table::data.table(
    gene_symbol = paste0("GENE", 1:20),
    prot_pattern = c(rep("Early_Burst_Up", 10),
                     rep("Early_Loss_Down", 10)))
  ev <- assemble_user_evidence(rna_fit = rna, prot_fit = prot)

  user_anchors <- data.table::data.table(
    gene = c("GENE1", "GENE2", "GENE3"),
    include_primary_eval = c(TRUE, TRUE, FALSE),
    include_secondary_eval = c(TRUE, TRUE, TRUE))

  res <- evaluate_anchor_enrichment(top_n = 10, tier = "secondary",
                                     evidence = ev,
                                     anchors = user_anchors)
  expect_true(nrow(res) >= 1L)
  expect_true(all(c("hits", "fold", "pval") %in% names(res)))
})

test_that("plot_gene_hexagon returns a ggplot for one or more genes", {
  p1 <- plot_gene_hexagon("LTBP1")
  expect_true(inherits(p1, "ggplot"))
  p2 <- plot_gene_hexagon(c("LGALS3BP", "LTBP1", "ALB", "GAPDH"))
  expect_true(inherits(p2, "ggplot"))
  p3 <- plot_gene_hexagon("LTBP1",
                            comparison = "high_confidence_mean")
  expect_true(inherits(p3, "ggplot"))
})

test_that("fit_stage_de_protein returns expected schema", {
  skip_if_not_installed("limma")
  skip_on_os("mac")
  set.seed(42)
  intensity <- matrix(rnorm(100 * 24), nrow = 100,
                       dimnames = list(paste0("G", 1:100), NULL))
  stage <- rep(c("Normal", "Early", "Mid", "Late"), each = 6)
  cohort <- rep(c("A", "B", "C"), times = 8)

  out <- fit_stage_de_protein(intensity, stage, cohort)
  expect_true(all(c("gene_symbol", "beta_N", "beta_E", "beta_M",
                    "beta_L", "lrt_padj", "lrt_significant")
                  %in% names(out)))
  expect_equal(out$beta_N, rep(0, nrow(out)))
})

# ── align_patient_profile() — single-patient trajectory alignment ──

test_that("align_patient_profile is deterministic given (input, atlas)", {
  ref <- pdactrace:::.get_reference()
  pool <- ref[!is.na(rna_beta_L) & audit_score >= 0.3][seq_len(200)]
  set.seed(42)
  patient <- setNames(pool$rna_beta_L + rnorm(nrow(pool), sd = 0.05),
                       pool$gene_symbol)
  a1 <- align_patient_profile(patient, top_n_genes = 150)
  a2 <- align_patient_profile(patient, top_n_genes = 150)
  expect_identical(a1$rna, a2$rna)
  expect_equal(a1$summary, a2$summary)
})

test_that("align_patient_profile recovers the synthetic Late axis", {
  ref <- pdactrace:::.get_reference()
  pool <- ref[!is.na(rna_beta_L) & audit_score >= 0.3][seq_len(200)]
  set.seed(7)
  patient <- setNames(0.1 * pool$rna_beta_L +
                          rnorm(nrow(pool), sd = 0.05),
                       pool$gene_symbol)
  out <- align_patient_profile(patient, top_n_genes = 150)
  rho_emL <- out$rna$cor_to_stage_axis[2:4]
  expect_equal(which.max(rho_emL), 3L)        # Late wins (E=2, M=3, L=4)
  expect_gt(out$rna[stage == "Late", cor_to_stage_axis], 0.5)
  expect_s3_class(out, "pdactrace_patient_alignment")
})

test_that("align_patient_profile errors on too few overlapping genes", {
  patient <- setNames(rnorm(20), paste0("FAKE_GENE_", seq_len(20)))
  expect_error(align_patient_profile(patient), regexp = ">=10")
})

test_that("align_patient_profile warns on low coverage", {
  ref <- pdactrace:::.get_reference()
  pool <- ref[!is.na(rna_beta_L) & audit_score >= 0.3][seq_len(30)]
  set.seed(11)
  patient <- setNames(pool$rna_beta_L + rnorm(30, sd = 0.05),
                       pool$gene_symbol)
  expect_warning(
    align_patient_profile(patient, top_n_genes = NULL,
                            min_genes = 50L),
    regexp = "min_genes")
})

test_that("align_patient_profile keeps protein layer separate, not combined", {
  ref <- pdactrace:::.get_reference()
  pool <- ref[!is.na(rna_beta_L) & audit_score >= 0.3][seq_len(200)]
  set.seed(99)
  patient_rna <- setNames(pool$rna_beta_L + rnorm(nrow(pool), sd = 0.1),
                           pool$gene_symbol)

  out_rna_only <- align_patient_profile(patient_rna, top_n_genes = 150)
  expect_null(out_rna_only$prot)

  prot_pool <- ref[!is.na(prot_pattern) &
                     prot_pattern %in% c("Early_Burst_Up",
                                          "Early_Loss_Down")]
  patient_prot <- setNames(
    ifelse(prot_pool$prot_pattern == "Early_Burst_Up", 1.0, -1.0) +
      rnorm(nrow(prot_pool), sd = 0.05),
    prot_pool$gene_symbol)
  out_both <- align_patient_profile(patient_rna,
                                      prot_logfc = patient_prot,
                                      top_n_genes = 150)
  expect_s3_class(out_both, "pdactrace_patient_alignment")
  expect_false(is.null(out_both$prot))
  expect_true(all(c("prot_pattern", "n_genes",
                     "n_concordant", "concordance")
                  %in% names(out_both$prot)))
})
