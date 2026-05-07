test_that("query_gene returns LTBP1 with canonical evidence", {
  q <- query_gene("LTBP1")
  expect_true(inherits(q, "pdactrace_gene_evidence"))
  expect_equal(q$rna$pattern, "Early_Burst_Up")
  expect_true(q$rna$pattern_rho > 0.95)
  expect_equal(q$serum$translation_class, "B")
  expect_true(q$serum$phase77_strict)
  expect_true(q$clinical$panel_member)
  expect_equal(q$clinical$evidence_scope, "panel_validated")
  expect_match(q$provenance, "phase33")
  expect_match(q$provenance, "phase77")
  expect_match(q$provenance, "phase80")
})

test_that("query_gene returns NULL with message for unknown gene", {
  expect_message(out <- query_gene("FOOBAR_NOT_A_GENE"),
                  "No evidence")
  expect_null(out)
})

test_that("query_panel joined returns one row per matched gene", {
  expect_message(qp <- query_panel(c("LTBP1", "SERPINA1", "FOOBAR")),
                  "skipped")
  expect_equal(nrow(qp), 2L)
  expect_setequal(qp$gene_symbol, c("LTBP1", "SERPINA1"))
  expect_true(all(c("rna_pattern", "translation_class", "panel_member")
                    %in% names(qp)))
})

test_that("list_candidates new multi-param signature", {
  # Use min_audit_class = "ALL" to test the filter logic without the
  # default v0.2.0 SILVER-tier gate; a separate test below verifies
  # the tier default itself.
  # Translation class aliases
  # v0.4.0 12-template change: GAPDH and CDH13 reclassified as
  # Monotonic_Up → dropped from Early × 4 surface (rna_pattern = NA).
  expect_equal(nrow(list_candidates(translation_class = "inverse",
                                       min_audit_class = "ALL")), 11L)
  expect_equal(nrow(list_candidates(translation_class = "B",
                                       min_audit_class = "ALL")), 11L)
  expect_equal(nrow(list_candidates(translation_class = "concordant",
                                       min_audit_class = "ALL")), 9L)
  expect_equal(nrow(list_candidates(translation_class = "A",
                                       min_audit_class = "ALL")), 9L)
  # Atlas-internal filters
  expect_equal(nrow(list_candidates(phase77_strict = TRUE,
                                       min_audit_class = "ALL")), 20L)
  expect_equal(nrow(list_candidates(panel_member = TRUE,
                                       min_audit_class = "ALL")), 4L)
  # Tissue direction split
  up_n   <- nrow(list_candidates(tissue_direction = "Up",
                                   min_audit_class = "ALL"))
  down_n <- nrow(list_candidates(tissue_direction = "Down",
                                   min_audit_class = "ALL"))
  expect_true(up_n > 100 && down_n > 100)
  # Compound filter
  out <- list_candidates(signal_peptide = TRUE,
                          panel_member = TRUE,
                          min_audit_class = "ALL")
  # v0.4.0: CDH13 reclassified as Monotonic_Up under 12-template, dropped.
  expect_setequal(out$gene_symbol,
                    c("LTBP1", "SERPINA1", "CP", "FGB"))
})

test_that("list_candidates min_audit_class filters by v0.3.0 audit_class", {
  # v0.4.0: confidence_tier removed. New filter is by v0.3.0 audit_class.
  # Default min_audit_class = "ALL" returns all 4 panel members.
  all_panel <- list_candidates(panel_member = TRUE)
  expect_setequal(all_panel$gene_symbol,
                    c("LTBP1", "SERPINA1", "CP", "FGB"))

  # min_audit_class = "high_confidence" drops everything below
  hc_only <- list_candidates(min_audit_class = "high_confidence")
  if (nrow(hc_only) > 0) {
    expect_true(all(hc_only$audit_class == "high_confidence"))
  }
})

test_that("list_candidates(top_n) caps result", {
  out <- list_candidates(top_n = 50)
  expect_equal(nrow(out), 50L)
  expect_true("provenance" %in% names(out))
})

test_that("trace_filters surfaces phase60 + class_route", {
  tf <- trace_filters(c("LTBP1", "SPARC", "SERPINA1", "CDH13"))
  expect_equal(nrow(tf), 4L)
  expect_true("class_route" %in% names(tf))
  # SERPINA1 is the phase60_final exemplar
  expect_equal(tf[gene_symbol == "SERPINA1", class_route],
                "phase60_final")
  expect_equal(tf[gene_symbol == "SERPINA1", n_phase60_pass], 7L)
  # LTBP1 / SPARC / CDH13 are phase77_classB
  expect_equal(tf[gene_symbol == "LTBP1", class_route],
                "phase77_classB")
  expect_equal(tf[gene_symbol == "SPARC", class_route],
                "phase77_classB")
})

test_that("Mid-pattern transparency: excluded_mid_pattern flag works", {
  ref <- pdactrace:::.get_reference()
  n_excluded <- sum(ref$excluded_mid_pattern, na.rm = TRUE)
  expect_true(n_excluded > 0)
  # All Mid-excluded rows must have NA rna_pattern
  expect_true(all(is.na(ref$rna_pattern[ref$excluded_mid_pattern])))
})
