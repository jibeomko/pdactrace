#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════
# Feature category ablation
#
# Each ablation: zero out one feature category (or set one
# multiplier to neutral) and recompute audit_score across atlas.
# Report: top-100 anchor enrichment + GAPDH/ALB/LTBP1 score change
# + Spearman vs full-rule + Top-100 Jaccard.
# ═══════════════════════════════════════════════════════════════
suppressPackageStartupMessages({ library(data.table) })
PKG <- rprojroot::find_package_root_file()
set.seed(42)

# ── Reference data ────────────────────────────────────────────
load(file.path(PKG, "data", "pdactrace_reference.rda"))
ref <- as.data.table(pdactrace_reference)

hk_list     <- fread(file.path(PKG, "data-raw",
                  "external_negatives_hpa_housekeeping.csv"))$gene_symbol
plasma_list <- setdiff(fread(file.path(PKG, "data-raw",
                  "external_negatives_hard_plasma.csv"))$gene_symbol,
                       c("APOA2","TTR"))
anc_secondary <- fread(file.path(PKG, "data-raw",
        "external_positives_anchors_v2.csv"))[include_secondary_eval == TRUE, gene]

EARLY <- c("Early_Burst_Up","Early_Loss_Down","Early_Peak","Early_Trough")
dir_of <- function(p) data.table::fcase(
  p %in% c("Early_Burst_Up","Early_Peak"), "UP",
  p %in% c("Early_Loss_Down","Early_Trough"), "DOWN",
  default = NA_character_)

# Pre-compute base
ref[, layer_rna := !is.na(rna_pattern)]
ref[, layer_prot := !is.na(prot_pattern)]
ref[, layer_scrna := !is.na(cell_origin_top)]
ref[, layer_serum := !is.na(serum_detected) & serum_detected == TRUE]
ref[, score_layer := (as.integer(layer_rna) + as.integer(layer_prot) +
                        as.integer(layer_scrna) + as.integer(layer_serum)) / 4]
ref[, dir_rna := dir_of(rna_pattern)]
ref[, dir_prot := dir_of(prot_pattern)]
ref[, cross_layer_concord := !is.na(dir_rna) & !is.na(dir_prot) & dir_rna == dir_prot]
ref[, cross_cohort_concord := pmax(0, pmin(1, fcoalesce(rna_cohort_agreement, 0)))]
ref[, score_direction := 0.5 * cross_cohort_concord + 0.5 * as.integer(cross_layer_concord)]
ref[, is_early := !is.na(rna_pattern) & rna_pattern %in% EARLY]
ref[, lrt_sig := pmin(1, -log10(pmax(rna_lrt_padj, 1e-10)) / 4)]
ref[, score_early := pmax(0, pmin(1,
       ifelse(is_early & !is.na(rna_pattern_rho),
              rna_pattern_rho * fcoalesce(lrt_sig, 0), 0)))]
ref[, score_serum := 0.4 * as.integer(layer_serum) +
                       0.3 * as.integer(!is.na(flt_signal_peptide) & flt_signal_peptide == TRUE) +
                       0.3 * as.integer(!is.na(flt_direction_match) & flt_direction_match == TRUE)]
ref[, rna_weak := confidence_tier %in% c("ARTIFACT","OTHER") | is.na(confidence_tier)]
ref[, other_layer_count := as.integer(layer_prot) + as.integer(layer_scrna) +
                              as.integer(layer_serum)]
ref[, rescue_eligible := rna_weak & other_layer_count >= 2 & cross_layer_concord]
ref[, score_rescue := as.integer(rescue_eligible)]
ref[, is_hk := gene_symbol %in% hk_list]
ref[, is_plasma_hi := gene_symbol %in% plasma_list]

# ── Score computation function ────────────────────────────────
compute <- function(w_layer = 0.20, w_dir = 0.20, w_early = 0.20,
                    w_serum = 0.10, w_rescue = 0.10,
                    leak_on = TRUE, het_on = TRUE) {
  pos <- w_layer * ref$score_layer +
         w_dir * ref$score_direction +
         w_early * ref$score_early +
         w_serum * ref$score_serum +
         w_rescue * ref$score_rescue
  leak <- if (leak_on) {
    ifelse(ref$is_hk, 0.00,
      ifelse(ref$is_plasma_hi, 0.50, 1.00))
  } else rep(1.00, nrow(ref))
  het <- if (het_on) {
    data.table::fcase(
      is.na(ref$max_I2_meta), 1.00,
      ref$max_I2_meta < 70, 1.00,
      ref$max_I2_meta < 90, 0.70,
      default = 0.30)
  } else rep(1.00, nrow(ref))
  raw <- pos * leak * het
  raw / max(raw, na.rm = TRUE)
}

# ── Reference (full rule) ──────────────────────────────────────
base_score <- compute()
base_top100 <- ref$gene_symbol[order(-base_score)][1:100]

eval_top100 <- function(scores) {
  ord <- order(-scores)
  top100 <- ref$gene_symbol[ord[1:100]]
  top500 <- ref$gene_symbol[ord[1:500]]
  hits100 <- sum(top100 %in% anc_secondary)
  expect <- 100 * length(intersect(anc_secondary, ref$gene_symbol)) / nrow(ref)
  fold <- hits100 / max(expect, 1e-12)
  pval <- phyper(hits100 - 1,
                  length(intersect(anc_secondary, ref$gene_symbol)),
                  nrow(ref) - length(intersect(anc_secondary, ref$gene_symbol)),
                  100, lower.tail = FALSE)
  jaccard <- length(intersect(top100, base_top100)) /
              length(union(top100, base_top100))
  spear <- suppressWarnings(cor(scores, base_score, method = "spearman",
                                 use = "pairwise.complete.obs"))
  list(hits = hits100, fold = round(fold, 1), pval = pval,
       jaccard = round(jaccard, 2), spearman = round(spear, 3),
       hk_top500 = sum(top500 %in% hk_list),
       plasma_top500 = sum(top500 %in% plasma_list),
       gapdh = round(scores[ref$gene_symbol == "GAPDH"], 3),
       alb = round(scores[ref$gene_symbol == "ALB"], 3),
       ltbp1 = round(scores[ref$gene_symbol == "LTBP1"], 3),
       ltbp1_rank = which(ref$gene_symbol[order(-scores)] == "LTBP1")[1])
}

# ── Run ablations ─────────────────────────────────────────────
runs <- list(
  "FULL_RULE"           = compute(),
  "no_layer_presence"   = compute(w_layer = 0),
  "no_direction_agree"  = compute(w_dir = 0),
  "no_early_pattern"    = compute(w_early = 0),
  "no_serum_bridge"     = compute(w_serum = 0),
  "no_rescue_signal"    = compute(w_rescue = 0),
  "no_leakage_gate"     = compute(leak_on = FALSE),
  "no_heterogeneity_gate" = compute(het_on = FALSE),
  "no_both_gates"       = compute(leak_on = FALSE, het_on = FALSE),
  "only_positives_eq"   = compute(),  # placeholder, default
  "weights_uniform"     = compute(w_layer = 0.16, w_dir = 0.16,
                                    w_early = 0.16, w_serum = 0.16,
                                    w_rescue = 0.16))

results <- list()
for (nm in names(runs)) {
  m <- eval_top100(runs[[nm]])
  m$ablation <- nm
  results[[nm]] <- as.data.table(m)
}
res <- rbindlist(results)
setcolorder(res, "ablation")

cat("========================================\n")
cat("=== Feature ablation summary ===\n")
cat("========================================\n")
print(res)

cat("\n=== Interpretation ===\n")
cat("Driver ranking (largest fold drop when ablated):\n")
full_fold <- res[ablation == "FULL_RULE", fold]
res[, fold_drop := full_fold - fold]
res[ablation != "FULL_RULE",
    print(.SD[order(-fold_drop), .(ablation, fold, fold_drop, jaccard, spearman,
                                       hk_top500, plasma_top500, gapdh, alb, ltbp1)])]

fwrite(res, file.path(PKG, "data-raw", "feature_ablation_results.csv"))
cat("\nSaved: data-raw/feature_ablation_results.csv\n")
