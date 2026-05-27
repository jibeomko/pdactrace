test_that("compute_trace_d returns the expected schema", {
  out <- compute_trace_d()
  expect_named(out,
               c("gene_symbol", "tracd_tissue_dir", "tracd_serum_dir",
                 "tracd_class", "tracd_confidence",
                 "tracd_pancreatitis_overlap_score",
                 "tracd_pancreatitis_specificity",
                 "tracd_tissue_weight", "tracd_decision_path"))
  expect_true(all(out$tracd_class %in% c("A", "B", "C", NA_character_)))
  expect_true(all(out$tracd_pancreatitis_specificity %in%
                    c("pdac_specific", "shared_inflammation",
                      "ambiguous")))
  expect_true(all(out$tracd_tissue_weight %in% c(0L, 1L, 2L)))
  conf <- out$tracd_confidence
  expect_true(all(is.na(conf) | (conf >= 0 & conf <= 1)))
})

test_that("compute_trace_d direction consensus + class assignment", {
  toy <- data.table::data.table(
    gene_symbol = c("A_clean", "B_inverse", "C_decoupled",
                    "NA_no_tissue", "A_shared_inflam",
                    "A_pdac_opp", "discordant"),
    rna_pattern = c("Early_Burst_Up", "Early_Burst_Up",
                    "Early_Burst_Up", NA_character_,
                    "Early_Burst_Up", "Early_Burst_Up",
                    "Early_Burst_Up"),
    prot_pattern = c("Early_Burst_Up", "Early_Burst_Up",
                     NA_character_, NA_character_,
                     "Early_Burst_Up", NA_character_,
                     "Early_Loss_Down"),
    max_abs_beta_meta = c(1.5, 1.2, 0.8, 0.1, 1.0, 0.9, 0.8),
    serum_log2fc_PDAC_vs_HC = c(0.8, -0.5, NA, 0.3, 0.6, 0.4, 0.5),
    serum_log2fc_Pan_vs_HC  = c(0.1, 0.0, NA, 0.0, 0.55, -0.2, 0.5),
    serum_detected = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
  )
  out <- compute_trace_d(atlas = toy)
  expect_equal(out$tracd_class,
               c("A", "B", "C", NA_character_, "A", "A",
                 NA_character_))
  expect_equal(out$tracd_tissue_weight,
               c(2L, 2L, 1L, 0L, 2L, 1L, 0L))
  expect_equal(out$tracd_pancreatitis_specificity[5L],
               "shared_inflammation")
  expect_equal(out$tracd_pancreatitis_specificity[6L],
               "pdac_specific")
  expect_equal(out$tracd_pancreatitis_specificity[3L], "ambiguous")
  expect_equal(out$tracd_serum_dir,
               c("UP", "DOWN", NA_character_, "UP", "UP", "UP", "UP"))
})

test_that("compute_trace_d confidence formula respects weights", {
  toy <- data.table::data.table(
    gene_symbol = "X",
    rna_pattern = "Early_Burst_Up",
    prot_pattern = "Early_Burst_Up",
    max_abs_beta_meta = 1.0,
    serum_log2fc_PDAC_vs_HC = 1.0,
    serum_log2fc_Pan_vs_HC = 0.0,
    serum_detected = TRUE
  )
  out <- compute_trace_d(atlas = toy)
  expect_equal(out$tracd_class, "A")
  # tissue_score = 1.0, serum_score = 1.0, overlap = 0 -> spec_score = 1.
  # conf = 0.50*1 + 0.30*1 + 0.20*1 = 1.0
  expect_equal(out$tracd_confidence, 1.0)
})

test_that("compute_trace_d threshold parameters are honoured", {
  toy <- data.table::data.table(
    gene_symbol = c("strong_tissue", "weak_tissue"),
    rna_pattern = c("Early_Burst_Up", "Early_Burst_Up"),
    prot_pattern = c("Early_Burst_Up", "Early_Burst_Up"),
    max_abs_beta_meta = c(0.6, 0.3),
    serum_log2fc_PDAC_vs_HC = c(0.5, 0.5),
    serum_log2fc_Pan_vs_HC = c(0.0, 0.0),
    serum_detected = c(TRUE, TRUE)
  )
  default <- compute_trace_d(atlas = toy)
  expect_equal(default$tracd_class, c("A", NA_character_))
  permissive <- compute_trace_d(atlas = toy, tau_tissue = 0.2)
  expect_equal(permissive$tracd_class, c("A", "A"))
})


test_that("compute_trace_d keeps strict and legacy fallback modes separate", {
  toy <- data.table::data.table(
    gene_symbol = c("legacy_A", "legacy_B", "strict_only"),
    rna_pattern = c("Early_Burst_Up", "Early_Burst_Up", "Early_Burst_Up"),
    prot_pattern = c("Early_Burst_Up", "Early_Burst_Up", "Early_Burst_Up"),
    max_abs_beta_meta = c(0.2, 0.2, 1.0),
    serum_log2fc_PDAC_vs_HC = c(NA_real_, NA_real_, 0.8),
    serum_log2fc_Pan_vs_HC = c(NA_real_, NA_real_, 0.0),
    serum_detected = c(TRUE, TRUE, TRUE),
    translation_class = c("A", "B", NA_character_)
  )
  strict <- compute_trace_d(atlas = toy)
  expect_equal(strict$tracd_class, c(NA_character_, NA_character_, "A"))

  compat <- compute_trace_d(atlas = toy, legacy_translation = "fallback")
  expect_equal(compat$tracd_class, c("A", "B", "A"))
  expect_match(compat$tracd_decision_path[1], "legacy_translation")
  expect_match(compat$tracd_decision_path[2], "legacy_translation")
  expect_false(grepl("legacy_translation", compat$tracd_decision_path[3]))
  expect_true(all(is.na(compat$tracd_confidence[1:2]) |
                    compat$tracd_confidence[1:2] <= 0.5))
})

test_that("compute_trace_d rejects invalid arguments", {
  expect_error(compute_trace_d(tau_tissue = -1),
               "non-negative")
  expect_error(compute_trace_d(weights = c(tissue = 0.5)),
               "named numeric")
  expect_error(compute_trace_d(weights = c(tissue = 0.4, serum = 0.4,
                                            specificity = 0.4)),
               "sum to 1")
  bad <- data.table::data.table(gene_symbol = "X")
  expect_error(compute_trace_d(atlas = bad),
               "missing required columns")
})

test_that("evaluate_anchor_enrichment accepts tracd_confidence", {
  skip_on_cran()
  ref <- pdactrace:::.get_reference(NULL)
  skip_if_not("tracd_confidence" %in% names(ref),
              "atlas does not carry tracd_confidence yet")
  res <- evaluate_anchor_enrichment(top_n = 100, tier = "secondary",
                                     score_col = "tracd_confidence")
  expect_true(nrow(res) == 1L)
  expect_true(is.numeric(res$fold))
  expect_true(res$pval >= 0 && res$pval <= 1)
})
