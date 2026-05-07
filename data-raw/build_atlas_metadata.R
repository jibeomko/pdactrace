#!/usr/bin/env Rscript
# ===========================================================
# Build atlas_metadata.rda  (v0.4.0 — 12-template + 3+2 audit)
#
# Cohort x layer x sample-n table — points users at source CSVs +
# documents the atlas snapshot.
# ===========================================================
suppressPackageStartupMessages({ library(data.table) })

PKG <- rprojroot::find_package_root_file()

atlas_metadata <- list(
  version       = "0.4.0",
  audit_rule_version = "v0.3.0 (frozen scoring rule)",
  snapshot_date = as.Date("2026-05-07"),
  scope         = "PDAC stage-aware tissue-to-serum biomarker atlas",
  classifier    = paste0(
    "12-template competitive trajectory matching ",
    "(Early x 4 + Mid x 4 + Late x 2 + Monotonic x 2); ",
    "atlas surface restricts visible rna_pattern / prot_pattern ",
    "calls to Early x 4 only. Non-Early best-matches are flagged ",
    "via excluded_mid_pattern, excluded_late_pattern, ",
    "excluded_monotonic_pattern."),
  pattern_scope = "Surfaced: Early x 4. Competing alternatives: Mid x 4 + Late x 2 + Monotonic x 2.",
  rna_cohorts = data.table(
    cohort     = c("TCGA", "CPTAC", "GSE71729", "GSE79668", "GSE93326",
                   "GSE119794", "GSE164665", "GSE224564", "GSE225767",
                   "GSE226762", "GSE253260"),
    samples    = NA_integer_,
    role       = "RNA-seq stage trajectory"),
  protein_cohorts = data.table(
    cohort  = c("CPTAC PDC000270 TMT", "CPTAC PDC000341 DIA",
                "KU PDC000248", "PXD015744", "PXD043111", "PXD048644"),
    samples = NA_integer_,
    role    = "Tissue proteomics stage trajectory"),
  scrna_cohorts = data.table(
    cohort  = "all_cohorts_scvi (10 cohort, ~372K cells)",
    samples = 372000L,
    role    = "scRNA cell origin"),
  serum_cohorts = data.table(
    cohort  = c("PXD048034", "PXD053603", "PXD066048",
                "PXD039273", "PXD046438", "PXD065581"),
    samples = NA_integer_,
    role    = "Serum proteomics MS detectability"),
  source_csvs = data.table(
    file = c("phase33_deseq2_coef_12template.csv",
             "phase34_protein_pooled_12template.csv",
             "multi_cohort_consistency.csv",
             "phase42_pancreatitis_check.csv",
             "phase77_strict_RNAprotConvergent_serum.csv",
             "phase29_resectable_markers.csv",
             "phase80_ltbp1_pancreatitis_predeclared_panels.csv",
             "phase2c_celltype_specificity.csv"),
    role = c("RNA 12-template best-match + LRT padj + per-stage Wald padj + b coefficients",
             "Tissue protein 12-template best-match (limma F-test)",
             "Multi-cohort Stouffer consistency (4 cohort)",
             "Serum log2FC PDAC/Pan/HC",
             "22 strict tissue+serum candidates",
             "25 phase29 resectable markers (s_pat 8-tpl discovery)",
             "LTBP1+SERPINA1 LOOCV panel (Class A x Class B)",
             "scRNA cell-type specificity"),
    note = c("Legacy alias *_8template.csv produced for back-compat",
             "Legacy alias *_8template.csv produced for back-compat",
             "", "", "", "", "", "")),
  audit_classes = c("high_confidence", "supported_uncertain",
                     "penalized", "excluded", "low"),
  caveats = c(
    paste0("v0.4.0 widens the catalog from 8 to 12 templates ",
            "(adds Late x 2 + Monotonic x 2). The atlas surface ",
            "restricts visible rna_pattern / prot_pattern calls ",
            "to Early x 4 only by design."),
    paste0("excluded_mid_pattern / excluded_late_pattern / ",
            "excluded_monotonic_pattern flag genes whose ",
            "12-template best-match was Mid / Late / Monotonic."),
    "resectable_pattern_phase29 retains the original phase29 s_pat label (8-tpl discovery).",
    "translation_class B is rare and manually curated, not predicted.",
    "Direction prediction is NOT included; reserved for v0.5+.",
    paste0("phase29 self-validation: discovery (s_pat, 8-tpl ",
            "single-cohort) 76% Early; production (coef_pat, ",
            "12-tpl cohort-adjusted) 53% Early in-atlas. Both ",
            "cited together to avoid reviewer audit risk.")))

out <- file.path(PKG, "data", "atlas_metadata.rda")
save(atlas_metadata, file = out, compress = "xz")
cat(sprintf("Saved: %s\n", out))
cat(sprintf("Version: %s, audit_rule: %s\n",
            atlas_metadata$version,
            atlas_metadata$audit_rule_version))
cat(sprintf("Classifier: 12-template (E x 4 + M x 4 + L x 2 + Mono x 2)\n"))
cat(sprintf("Atlas surface: Early x 4 only\n"))
cat(sprintf("Audit labels (%d): %s\n",
            length(atlas_metadata$audit_classes),
            paste(atlas_metadata$audit_classes, collapse = ", ")))
