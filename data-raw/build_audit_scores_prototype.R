#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════
# data-raw/build_audit_scores_prototype.R   (pdactrace v0.3.0 X)
#
# Hand-engineered transparent audit_score (NO supervised training).
# Formula:
#   audit_score = positive_score × leakage_multiplier × heterogeneity_multiplier
# Weights locked in audit_score_weights_v1.csv (supplement Table S1).
# ═══════════════════════════════════════════════════════════════
suppressPackageStartupMessages({ library(data.table) })

PKG <- rprojroot::find_package_root_file()
load(file.path(PKG, "data", "pdactrace_reference.rda"))
ref <- as.data.table(pdactrace_reference)

# ── Reference sets ─────────────────────────────────────────────
hk_list     <- fread(file.path(PKG, "data-raw",
                    "external_negatives_hpa_housekeeping.csv"))$gene_symbol
plasma_list <- fread(file.path(PKG, "data-raw",
                    "external_negatives_hard_plasma.csv"))$gene_symbol
anchor_pos  <- fread(file.path(PKG, "data-raw",
                    "external_positives_tier1_literature.csv"))$gene_symbol

# Drop dual-role from negative side
plasma_list <- setdiff(plasma_list, c("APOA2","TTR"))

# ── Helpers ───────────────────────────────────────────────────
EARLY <- c("Early_Burst_Up","Early_Loss_Down","Early_Peak","Early_Trough")
direction_of <- function(p) {
  data.table::fcase(
    p %in% c("Early_Burst_Up","Early_Peak"), "UP",
    p %in% c("Early_Loss_Down","Early_Trough"), "DOWN",
    default = NA_character_)
}

# ── 1. Layer presence (weight 0.20) ───────────────────────────
ref[, layer_rna   := !is.na(rna_pattern)]
ref[, layer_prot  := !is.na(prot_pattern)]
ref[, layer_scrna := !is.na(cell_origin_top)]
ref[, layer_serum := !is.na(serum_detected) & serum_detected == TRUE]
ref[, layer_count := as.integer(layer_rna) + as.integer(layer_prot) +
                       as.integer(layer_scrna) + as.integer(layer_serum)]
ref[, score_layer := layer_count / 4]   # 0-1

# ── 2. Direction agreement (weight 0.20) ──────────────────────
ref[, dir_rna  := direction_of(rna_pattern)]
ref[, dir_prot := direction_of(prot_pattern)]
ref[, cross_layer_concord := !is.na(dir_rna) & !is.na(dir_prot) & dir_rna == dir_prot]
ref[, cross_cohort_concord := pmax(0, pmin(1,
       data.table::fcoalesce(rna_cohort_agreement, 0)))]
ref[, score_direction := 0.5 * cross_cohort_concord +
                            0.5 * as.integer(cross_layer_concord)]

# ── 3. Early-pattern concentration (weight 0.20) ──────────────
ref[, is_early := !is.na(rna_pattern) & rna_pattern %in% EARLY]
ref[, lrt_sig_factor := pmin(1, -log10(pmax(rna_lrt_padj, 1e-10, na.rm=TRUE)) / 4)]
ref[, score_early := ifelse(is_early & !is.na(rna_pattern_rho),
                              rna_pattern_rho * data.table::fcoalesce(lrt_sig_factor, 0),
                              0)]
ref[, score_early := pmax(0, pmin(1, score_early))]

# ── 4. Serum bridge (weight 0.10) ─────────────────────────────
ref[, score_serum := 0.4 * as.integer(layer_serum) +
                        0.3 * as.integer(!is.na(flt_signal_peptide) & flt_signal_peptide == TRUE) +
                        0.3 * as.integer(!is.na(flt_direction_match) & flt_direction_match == TRUE)]

# ── 5. Rescue signal (weight 0.10, strict eligibility) ────────
# RNA weak/artifact + ≥2 other layer + cross-layer direction concord
ref[, rna_weak := confidence_tier %in% c("ARTIFACT","OTHER") | is.na(confidence_tier)]
ref[, other_layer_count := as.integer(layer_prot) + as.integer(layer_scrna) +
                              as.integer(layer_serum)]
ref[, rescue_eligible := rna_weak & other_layer_count >= 2 & cross_layer_concord]
ref[, score_rescue := as.integer(rescue_eligible)]

# ── Combine positive components ───────────────────────────────
ref[, positive_score := 0.20 * score_layer +
                          0.20 * score_direction +
                          0.20 * score_early +
                          0.10 * score_serum +
                          0.10 * score_rescue]
# range 0-0.80 by construction

# ── Multipliers ───────────────────────────────────────────────
ref[, is_hk         := gene_symbol %in% hk_list]
ref[, is_plasma_hi  := gene_symbol %in% plasma_list]
ref[, leakage_mult  := data.table::fcase(
       is_hk, 0.00,
       is_plasma_hi, 0.50,
       default = 1.00)]

ref[, het_mult := data.table::fcase(
       is.na(max_I2_meta), 1.00,
       max_I2_meta < 50, 1.00,
       max_I2_meta < 70, 1.00,
       max_I2_meta < 90, 0.70,
       default = 0.30)]

# ── Final ─────────────────────────────────────────────────────
ref[, audit_score_raw := positive_score * leakage_mult * het_mult]
mx <- max(ref$audit_score_raw, na.rm = TRUE)
ref[, audit_score := audit_score_raw / mx]

# ═══════════════════════════════════════════════════════════════
#  Reporting
# ═══════════════════════════════════════════════════════════════
cat("\n========================================\n")
cat("=== Atlas-wide audit_score distribution ===\n")
cat("========================================\n")
print(summary(ref$audit_score))
cat(sprintf("audit_score == 0 : %d (%.1f%%)\n",
    sum(ref$audit_score == 0, na.rm = TRUE),
    100 * sum(ref$audit_score == 0, na.rm = TRUE) / nrow(ref)))

cat("\n========================================\n")
cat("=== 6 case study genes ===\n")
cat("========================================\n")
cs <- c("LTBP1","GAPDH","ALB","C6","AMBP","SERPINA1","ANPEP","TIMP1","THBS2")
print(ref[gene_symbol %in% cs,
          .(gene_symbol,
            pos = round(positive_score, 3),
            leak = leakage_mult,
            het = het_mult,
            score = round(audit_score, 3),
            tier_v0.2 = confidence_tier)][order(-score)])

cat("\n========================================\n")
cat("=== Top 30 audit_score (anchor enriched?) ===\n")
cat("========================================\n")
top30 <- ref[order(-audit_score)][1:30,
              .(gene_symbol,
                score = round(audit_score, 3),
                anchor = gene_symbol %in% anchor_pos,
                hk = is_hk,
                plasma = is_plasma_hi,
                tier_v0.2 = confidence_tier)]
print(top30)

cat("\n========================================\n")
cat("=== Anchor enrichment evaluation ===\n")
cat("========================================\n")
n_atlas      <- nrow(ref)
n_anchor_in  <- sum(ref$gene_symbol %in% anchor_pos)
for (top_n in c(50, 100, 200, 500, 1000)) {
  top_set <- ref[order(-audit_score)][1:top_n, gene_symbol]
  hits    <- sum(top_set %in% anchor_pos)
  expect  <- top_n * (n_anchor_in / n_atlas)
  fold    <- hits / expect
  pval    <- phyper(hits - 1, n_anchor_in, n_atlas - n_anchor_in, top_n,
                    lower.tail = FALSE)
  cat(sprintf("Top %4d : %2d hits (expect %.2f, fold=%.1fx, hyper p=%.2e)\n",
              top_n, hits, expect, fold, pval))
}

cat("\n========================================\n")
cat("=== Negative reference reject evaluation ===\n")
cat("========================================\n")
hk_in_atlas     <- intersect(hk_list, ref$gene_symbol)
plasma_in_atlas <- intersect(plasma_list, ref$gene_symbol)
neg_all <- unique(c(hk_in_atlas, plasma_in_atlas))
for (top_n in c(100, 500, 1000)) {
  top_set <- ref[order(-audit_score)][1:top_n, gene_symbol]
  hk_pollution     <- sum(top_set %in% hk_in_atlas)
  plasma_pollution <- sum(top_set %in% plasma_in_atlas)
  cat(sprintf("Top %4d : %d HK + %d plasma high-abund   (negatives in top : %d / %d)\n",
              top_n, hk_pollution, plasma_pollution,
              hk_pollution + plasma_pollution, length(neg_all)))
}

# v0.2 false-positive 18 reduction
cat("\n========================================\n")
cat("=== v0.2 false-positive 18 reduction ===\n")
cat("========================================\n")
fp_v02 <- ref[(gene_symbol %in% hk_in_atlas | gene_symbol %in% plasma_in_atlas) &
                confidence_tier %in% c("GOLD","SILVER"), gene_symbol]
cat(sprintf("v0.2 false positives (HK/plasma in GOLD/SILVER tier): %d genes\n",
            length(fp_v02)))
cat(sprintf("  list: %s\n", paste(fp_v02, collapse=", ")))
fp_v02_in_top500 <- sum(fp_v02 %in% ref[order(-audit_score)][1:500, gene_symbol])
cat(sprintf("Of these %d v0.2 FP, how many in NEW audit_score top 500: %d\n",
            length(fp_v02), fp_v02_in_top500))
fp_v02_in_top100 <- sum(fp_v02 %in% ref[order(-audit_score)][1:100, gene_symbol])
cat(sprintf("                                              top 100: %d\n",
            fp_v02_in_top100))

# Save attached cols
ATTACH_COLS <- c("score_layer","score_direction","score_early","score_serum",
                  "score_rescue","positive_score","leakage_mult","het_mult",
                  "audit_score_raw","audit_score","is_hk","is_plasma_hi",
                  "rescue_eligible")
out_path <- file.path(PKG, "data-raw", "audit_score_atlas_legacy_7feature.csv")
fwrite(ref[, c("gene_symbol", ATTACH_COLS), with = FALSE], out_path)
cat(sprintf("\nSaved: %s (atlas + audit_score columns)\n", out_path))
