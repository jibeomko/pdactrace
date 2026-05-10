test_that("evidence_math returns a list with all 8 axes named", {
  m <- evidence_math("LGALS3BP")
  expect_type(m, "list")
  expect_named(m, c("gene", "trajectory_fit", "effect_magnitude",
                     "cohort_consistency", "rna_protein_coupling",
                     "serum_bridge", "cell_specificity",
                     "filter_survival", "clinical_role"),
               ignore.order = TRUE)
})

test_that("delta_rho matches rna_pattern_rho - rna_pattern_rho_runner_up from atlas", {
  data("pdactrace_reference", package = "pdactrace")
  ref <- pdactrace_reference
  for (g in c("LGALS3BP", "LTBP1", "TIMP1")) {
    m <- evidence_math(g)
    target <- g
    expected <- ref[gene_symbol == target,
                    rna_pattern_rho - rna_pattern_rho_runner_up]
    expect_equal(m$trajectory_fit$delta_rho, expected,
                 tolerance = 1e-8)
  }
})

test_that("rna_beta_norm equals sqrt(sum(beta_E,M,L^2))", {
  data("pdactrace_reference", package = "pdactrace")
  ref <- pdactrace_reference
  g <- "LGALS3BP"
  target <- g
  betas <- as.numeric(ref[gene_symbol == target,
                          c(rna_beta_E, rna_beta_M, rna_beta_L)])
  m <- evidence_math(g)
  expect_equal(m$effect_magnitude$rna_beta_norm,
               sqrt(sum(betas^2)), tolerance = 1e-8)
})

test_that("RNA-protein cosine for LGALS3BP is high and positive", {
  m <- evidence_math("LGALS3BP")
  expect_true(isTRUE(m$rna_protein_coupling$prot_in_atlas))
  expect_gte(m$rna_protein_coupling$cosine, 0.8)
  expect_lte(m$rna_protein_coupling$cosine, 1.0)
})

test_that("RNA-protein cosine returns NA when gene missing from protein betas", {
  data("pdactrace_reference", package = "pdactrace")
  data("pdactrace_protein_betas", package = "pdactrace")
  rna_only <- setdiff(pdactrace_reference$gene_symbol,
                       pdactrace_protein_betas$gene_symbol)
  if (length(rna_only) == 0L) skip("no RNA-only genes in this atlas")
  g <- rna_only[1L]
  m <- evidence_math(g)
  expect_false(isTRUE(m$rna_protein_coupling$prot_in_atlas))
  expect_true(is.na(m$rna_protein_coupling$cosine))
})

test_that("filter_survival passed count matches sum(flt_*)", {
  data("pdactrace_reference", package = "pdactrace")
  ref <- pdactrace_reference
  g <- "LTBP1"
  target <- g
  fls <- ref[gene_symbol == target,
             c(flt_signal_peptide, flt_serum_measurable,
               flt_serum_significant, flt_pancreatitis_pdac,
               flt_pancreatitis_hc, flt_direction_match, flt_final)]
  m <- evidence_math(g)
  expect_equal(m$filter_survival$passed,
               as.integer(sum(fls, na.rm = TRUE)))
  expect_equal(m$filter_survival$total, 7L)
})

test_that("evidence_math errors clearly on missing gene", {
  expect_error(evidence_math("NOT_A_REAL_GENE_XYZ"),
               regexp = "not in the bundled atlas")
})

test_that("evidence_math is reproducible across two calls", {
  m1 <- evidence_math("LGALS3BP")
  m2 <- evidence_math("LGALS3BP")
  attr(m1, "reference_version") <- NULL
  attr(m2, "reference_version") <- NULL
  expect_identical(m1, m2)
})
