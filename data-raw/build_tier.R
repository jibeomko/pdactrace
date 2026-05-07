#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════
# data-raw/build_tier.R  (DEPRECATED in v0.4.0)
#
# ⚠️  WARNING: This script computes columns that were REMOVED from
# the bundled atlas in v0.4.0:
#   - confidence_tier      (5-label GOLD/SILVER/BRONZE/OTHER/ARTIFACT)
#   - early_onset_score    (composite tier score)
#   - heterogeneity_factor (only used by tier-system internals)
#
# In v0.4.0 the canonical per-gene classification became `audit_class`
# (5-label, v0.3.0 frozen 3+2 framework: high_confidence /
# supported_uncertain / penalized / excluded / low). The dual
# v0.2.0-tier vs v0.3.0-class system was creating reviewer-grade
# audit ambiguity (e.g., LGALS3BP scoring `audit_class =
# high_confidence` for strong multi-layer convergence while
# `confidence_tier = ARTIFACT` for low single-layer effect size).
#
# If you re-run this script today it will recompute the legacy tier
# columns BUT they will be dropped immediately by
# `data-raw/drop_v0_2_tier_columns.R` (which runs after this in the
# v0.4.0 rebuild chain). The ONLY columns this script still
# legitimately produces are the two audit-framework dependencies:
#   - max_abs_beta_meta  (used by audit_score)
#   - max_I2_meta        (used by heterogeneity_gate)
#
# Recommended v0.4.0 rebuild path: skip this script. The two
# audit-framework columns are populated by `build_meta_analysis.R`
# directly, and the tier columns are no longer needed.
#
# This file is retained for the audit trail of how v0.2.0 tier was
# computed. Source-of-truth tier rules (legacy):
#
#   EarlyOnsetScore = rho × max|β_meta| × √(1 − max_I²/100) × cohort_agreement
#
#   GOLD     — rho≥0.85, max|β|≥1.0,   max_I²<50, agreement≥0.75
#   SILVER   — rho≥0.85, max|β|≥0.585, (max_I²<70 OR agreement≥0.5)
#   BRONZE   — rho≥0.75, max|β|≥0.585
#   ARTIFACT — rho≥0.85, max|β|<0.585, Stouffer p>0.5
#   OTHER    — borderline / unclassified
# ═══════════════════════════════════════════════════════════════
suppressPackageStartupMessages({ library(data.table) })
PKG <- rprojroot::find_package_root_file()

cat("=== Loading reference + meta_analysis ===\n")
load(file.path(PKG, "data", "pdactrace_reference.rda"))
load(file.path(PKG, "data", "meta_analysis.rda"))
cat(sprintf("  reference: %d genes × %d cols\n",
            nrow(pdactrace_reference), ncol(pdactrace_reference)))
cat(sprintf("  meta:      %d genes × %d cols\n",
            nrow(meta_analysis), ncol(meta_analysis)))

# Subset meta to columns we want to expose at atlas level
meta_keep <- meta_analysis[, .(
  gene_symbol,
  meta_NvE_beta = meta_Normal_vs_Early_beta,
  meta_NvE_se   = meta_Normal_vs_Early_se,
  meta_NvE_pval = meta_Normal_vs_Early_pval,
  meta_NvE_padj = meta_Normal_vs_Early_padj,
  meta_NvE_I2   = meta_Normal_vs_Early_I2,
  meta_NvE_k    = meta_Normal_vs_Early_k,
  meta_MvE_beta = meta_Mid_vs_Early_beta,
  meta_MvE_pval = meta_Mid_vs_Early_pval,
  meta_MvE_padj = meta_Mid_vs_Early_padj,
  meta_MvE_I2   = meta_Mid_vs_Early_I2,
  meta_MvE_k    = meta_Mid_vs_Early_k,
  meta_LvE_beta = meta_Late_vs_Early_beta,
  meta_LvE_pval = meta_Late_vs_Early_pval,
  meta_LvE_padj = meta_Late_vs_Early_padj,
  meta_LvE_I2   = meta_Late_vs_Early_I2,
  meta_LvE_k    = meta_Late_vs_Early_k,
  meta_cohort_divergent)]

# Join into reference
ref <- merge(pdactrace_reference, meta_keep, by = "gene_symbol",
             all.x = TRUE, sort = FALSE)
data.table::setkey(ref, gene_symbol)

# ── Composite score ─────────────────────────────────────────
ref[, max_abs_beta_meta := pmax(abs(meta_NvE_beta), abs(meta_MvE_beta),
                                  abs(meta_LvE_beta), na.rm = TRUE)]
ref[, max_I2_meta := pmax(meta_NvE_I2, meta_MvE_I2, meta_LvE_I2,
                           na.rm = TRUE)]
ref[, heterogeneity_factor := sqrt(pmax(0, 1 - max_I2_meta / 100))]
ref[, early_onset_score := rna_pattern_rho * max_abs_beta_meta *
                            heterogeneity_factor * rna_cohort_agreement]

# ── Tier classification ─────────────────────────────────────
ref[, confidence_tier := data.table::fcase(
  # GOLD: rho>=0.85, |β|>=1.0, I²<50, agreement>=75%
  !is.na(rna_pattern) & !is.na(max_abs_beta_meta) & !is.na(max_I2_meta) &
    rna_pattern_rho >= 0.85 & max_abs_beta_meta >= 1.0 &
    max_I2_meta < 50 & rna_cohort_agreement >= 0.75, "GOLD",
  # SILVER: rho>=0.85, |β|>=0.585, (I²<70 OR agree>=50%)
  !is.na(rna_pattern) & !is.na(max_abs_beta_meta) &
    rna_pattern_rho >= 0.85 & max_abs_beta_meta >= 0.585 &
    (max_I2_meta < 70 | rna_cohort_agreement >= 0.5), "SILVER",
  # BRONZE: rho>=0.75, |β|>=0.585
  !is.na(rna_pattern) & !is.na(max_abs_beta_meta) &
    rna_pattern_rho >= 0.75 & max_abs_beta_meta >= 0.585, "BRONZE",
  # ARTIFACT: rho>=0.85, |β|<0.585, Stouffer p>0.5
  !is.na(rna_pattern) & !is.na(max_abs_beta_meta) &
    rna_pattern_rho >= 0.85 & max_abs_beta_meta < 0.585 &
    rna_stouffer_p > 0.5, "ARTIFACT",
  # OTHER: borderline / has pattern but not in above categories
  !is.na(rna_pattern), "OTHER",
  default = NA_character_)]

cat("\n=== Composite score distribution ===\n")
cat(sprintf("  with score: %d / %d\n",
            sum(!is.na(ref$early_onset_score)), nrow(ref)))
print(summary(ref$early_onset_score))

cat("\n=== Tier breakdown ===\n")
print(table(ref$confidence_tier, useNA = "ifany"))

cat("\n=== Headline genes — tier check ===\n")
hl <- ref[gene_symbol %in% c("LTBP1","SERPINA1","CDH13","SPARC","CP","FGB",
                              "ACTB","GAPDH","GUSB","F5","C7","C2","C9","HBB",
                              "PON1","KLKB1","APOA1","AMBP","VNN1","CTSD","A1BG"),
          .(gene = gene_symbol,
            pat = rna_pattern,
            rho = round(rna_pattern_rho,2),
            maxBeta_meta = round(max_abs_beta_meta,2),
            max_I2 = round(max_I2_meta),
            agree = rna_cohort_agreement,
            score = round(early_onset_score,2),
            tier = confidence_tier)][order(tier, -score)]
print(hl)

# ── Save updated reference ──────────────────────────────────
pdactrace_reference <- ref
out_path <- file.path(PKG, "data", "pdactrace_reference.rda")
save(pdactrace_reference, file = out_path, compress = "xz")
cat(sprintf("\nSaved: %s\n", out_path))
cat(sprintf("Updated reference: %d genes × %d cols\n",
            nrow(pdactrace_reference), ncol(pdactrace_reference)))
