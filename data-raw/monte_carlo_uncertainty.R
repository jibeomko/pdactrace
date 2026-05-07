#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════
# Monte Carlo uncertainty propagation for audit_score
#
# Scoring rule: FROZEN at v1 freeze (compute_audit_score logic).
# MC perturbs underlying *evidence* (rho, agreement, I²), NOT weights.
#
# Rules (사용자 lock):
#   - leakage multiplier: NOT perturbed (HK gate stays 0, plasma stays 0.5)
#   - weights: NOT perturbed
#   - heterogeneity → uncertainty width via I² band crossing
#   - GAPDH/HK: CI fixed [0, 0]
#
# Output columns:
#   audit_score (median), audit_score_lo95, audit_score_hi95
#   rank_median, rank_lo95, rank_hi95
#   uncertainty_width
#   confidence_class: excluded / stable_high / high_uncertain / medium / low
# ═══════════════════════════════════════════════════════════════
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
PKG <- rprojroot::find_package_root_file()
set.seed(42)
N_MC <- 500

# ── Load atlas + frozen leakage flags + base features ─────────
load(file.path(PKG, "data", "pdactrace_reference.rda"))
ref <- as.data.table(pdactrace_reference)

hk_list     <- fread(file.path(PKG, "data-raw",
                    "external_negatives_hpa_housekeeping.csv"))$gene_symbol
plasma_list <- setdiff(fread(file.path(PKG, "data-raw",
                    "external_negatives_hard_plasma.csv"))$gene_symbol,
                       c("APOA2","TTR"))

EARLY <- c("Early_Burst_Up","Early_Loss_Down","Early_Peak","Early_Trough")
dir_of <- function(p) data.table::fcase(
  p %in% c("Early_Burst_Up","Early_Peak"), "UP",
  p %in% c("Early_Loss_Down","Early_Trough"), "DOWN",
  default = NA_character_)

# Pre-compute deterministic features
ref[, layer_rna   := !is.na(rna_pattern)]
ref[, layer_prot  := !is.na(prot_pattern)]
ref[, layer_scrna := !is.na(cell_origin_top)]
ref[, layer_serum := !is.na(serum_detected) & serum_detected == TRUE]
ref[, score_layer := (as.integer(layer_rna) + as.integer(layer_prot) +
                        as.integer(layer_scrna) + as.integer(layer_serum)) / 4]
ref[, dir_rna  := dir_of(rna_pattern)]
ref[, dir_prot := dir_of(prot_pattern)]
ref[, cross_layer_concord := !is.na(dir_rna) & !is.na(dir_prot) & dir_rna == dir_prot]
ref[, is_early := !is.na(rna_pattern) & rna_pattern %in% EARLY]
ref[, lrt_sig := pmin(1, -log10(pmax(rna_lrt_padj, 1e-10)) / 4)]
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
ref[, leakage_mult := data.table::fcase(
       is_hk, 0.00,
       is_plasma_hi, 0.50,
       default = 1.00)]

n_genes <- nrow(ref)
score_mc <- matrix(NA_real_, nrow = n_genes, ncol = N_MC)
rank_mc  <- matrix(NA_integer_, nrow = n_genes, ncol = N_MC)

# ── MC loop (vectorized over genes) ───────────────────────────
cat(sprintf("Running %d MC iterations over %d genes ...\n", N_MC, n_genes))
t0 <- Sys.time()
rho_se <- 0.05  # ρ standard error proxy (n_cohorts ≈ 4)
i2_se  <- 8     # I² standard error approximation (Higgins ~5-10)

# Pre-extract base values
rho_b <- ref$rna_pattern_rho
agree_b <- ref$rna_cohort_agreement
i2_b <- ref$max_I2_meta
lrt_b <- ref$lrt_sig
is_early_b <- ref$is_early
cross_concord_b <- as.integer(ref$cross_layer_concord)
sl_b <- ref$score_layer
ss_b <- ref$score_serum
sr_b <- ref$score_rescue
leak_b <- ref$leakage_mult

for (iter in 1:N_MC) {
  # Perturb ρ via Fisher-z (ρ ∈ [-1,1])
  rho_p <- ifelse(is.na(rho_b), NA_real_,
                   tanh(atanh(pmin(0.999, pmax(-0.999, rho_b))) +
                        rnorm(n_genes, 0, rho_se)))
  rho_p <- pmax(-1, pmin(1, rho_p))

  # Perturb cohort_agreement via Beta(k+1, n-k+1), n_cohorts=4 assumption
  k <- pmax(0, pmin(4, round(data.table::fcoalesce(agree_b, 0) * 4)))
  agree_p <- rbeta(n_genes, k + 1, 4 - k + 1)

  # Perturb I² (Higgins SE approx)
  i2_p <- ifelse(is.na(i2_b), NA_real_,
                  pmax(0, pmin(100, i2_b + rnorm(n_genes, 0, i2_se))))

  # Recompute MC-dependent components
  se_iter <- ifelse(is_early_b & !is.na(rho_p),
                     pmax(0, pmin(1, rho_p * data.table::fcoalesce(lrt_b, 0))), 0)
  sd_iter <- 0.5 * agree_p + 0.5 * cross_concord_b

  pos <- 0.20 * sl_b +
         0.20 * sd_iter +
         0.20 * se_iter +
         0.10 * ss_b +
         0.10 * sr_b

  het_mult_iter <- data.table::fcase(
    is.na(i2_p), 1.00,
    i2_p < 50, 1.00,
    i2_p < 70, 1.00,
    i2_p < 90, 0.70,
    default = 0.30)

  raw <- pos * leak_b * het_mult_iter
  score_mc[, iter] <- raw / max(raw, na.rm = TRUE)
  rank_mc[, iter]  <- frank(-score_mc[, iter], na.last = "keep")
  if (iter %% 100 == 0) cat(sprintf("  iter %d / %d (elapsed %.1fs)\n",
                                     iter, N_MC, as.numeric(Sys.time() - t0)))
}
cat(sprintf("MC done. Elapsed: %.1fs\n", as.numeric(Sys.time() - t0)))

# ── Quantile summaries ────────────────────────────────────────
ref[, audit_score_median := apply(score_mc, 1, median, na.rm = TRUE)]
ref[, audit_score_lo95   := apply(score_mc, 1, quantile, probs = 0.025, na.rm = TRUE)]
ref[, audit_score_hi95   := apply(score_mc, 1, quantile, probs = 0.975, na.rm = TRUE)]
ref[, uncertainty_width  := audit_score_hi95 - audit_score_lo95]
ref[, rank_median := apply(rank_mc, 1, median, na.rm = TRUE)]
ref[, rank_lo95   := apply(rank_mc, 1, quantile, probs = 0.025, na.rm = TRUE)]
ref[, rank_hi95   := apply(rank_mc, 1, quantile, probs = 0.975, na.rm = TRUE)]

# Confidence class
ref[, confidence_class := data.table::fcase(
       is_hk | (audit_score_hi95 == 0 & audit_score_lo95 == 0), "excluded",
       audit_score_lo95 >= 0.5, "stable_high",
       audit_score_hi95 >= 0.5 & audit_score_lo95 < 0.5, "high_uncertain",
       audit_score_lo95 >= 0.3, "medium",
       default = "low")]

# ── Reports ───────────────────────────────────────────────────
cat("\n========================================\n")
cat("=== Confidence class breakdown ===\n")
cat("========================================\n")
print(table(ref$confidence_class))

cat("\n========================================\n")
cat("=== Case study CI ===\n")
cat("========================================\n")
cs <- c("LTBP1","GAPDH","ALB","C6","AMBP","SERPINA1","ANPEP","TIMP1","THBS2","MUC16")
print(ref[gene_symbol %in% cs,
          .(gene_symbol,
            score_med = round(audit_score_median, 3),
            lo95 = round(audit_score_lo95, 3),
            hi95 = round(audit_score_hi95, 3),
            width = round(uncertainty_width, 3),
            rank_med = rank_median,
            rank_lo = rank_lo95,
            rank_hi = rank_hi95,
            class = confidence_class)][order(-score_med)])

cat("\n========================================\n")
cat("=== Top 30 by audit_score_median (with CI) ===\n")
cat("========================================\n")
print(ref[order(-audit_score_median)][1:30,
          .(gene_symbol,
            score = round(audit_score_median, 3),
            ci = sprintf("[%.2f, %.2f]", audit_score_lo95, audit_score_hi95),
            width = round(uncertainty_width, 3),
            class = confidence_class)])

# ── Save atlas with MC results ────────────────────────────────
out_path <- file.path(PKG, "data-raw", "audit_score_mc_v1.csv")
fwrite(ref[, .(gene_symbol,
               audit_score_median, audit_score_lo95, audit_score_hi95,
               uncertainty_width,
               rank_median, rank_lo95, rank_hi95,
               confidence_class,
               is_hk, is_plasma_hi, leakage_mult,
               max_I2_meta, rna_pattern_rho, rna_cohort_agreement)],
       out_path)

# ── Figures ───────────────────────────────────────────────────
# Fig 1: Top 50 candidates score + 95% CI
top50 <- ref[order(-audit_score_median)][1:50]
top50[, gene_factor := factor(gene_symbol, levels = rev(gene_symbol))]
p1 <- ggplot(top50, aes(x = audit_score_median, y = gene_factor)) +
  geom_errorbarh(aes(xmin = audit_score_lo95, xmax = audit_score_hi95,
                      color = confidence_class), height = 0, linewidth = 0.4) +
  geom_point(aes(color = confidence_class), size = 1.5) +
  scale_color_manual(values = c("stable_high" = "#1B5E20", "high_uncertain" = "#F57C00",
                                  "medium" = "#1976D2", "low" = "grey60",
                                  "excluded" = "black")) +
  labs(title = "Top 50 audit_score with 95% CI (Monte Carlo, N=500)",
       x = "audit_score (median ± 95% CI)", y = NULL, color = "class") +
  theme_minimal(base_size = 8)
ggsave(file.path(PKG, "data-raw", "mc_top50_ci.png"), p1, width = 6, height = 9)

# Fig 2: Case study CI
cs_dt <- ref[gene_symbol %in% c("LTBP1","GAPDH","THBS2","TIMP1","AMBP","SERPINA1","ALB")]
cs_dt[, gene_factor := factor(gene_symbol,
                                 levels = c("THBS2","TIMP1","AMBP","SERPINA1",
                                            "LTBP1","ALB","GAPDH"))]
p2 <- ggplot(cs_dt, aes(x = audit_score_median, y = gene_factor)) +
  geom_errorbarh(aes(xmin = audit_score_lo95, xmax = audit_score_hi95,
                      color = confidence_class), height = 0.2, linewidth = 0.6) +
  geom_point(aes(color = confidence_class), size = 3) +
  scale_color_manual(values = c("stable_high" = "#1B5E20", "high_uncertain" = "#F57C00",
                                  "medium" = "#1976D2", "low" = "grey60",
                                  "excluded" = "black")) +
  labs(title = "Case-study CI: anchors + LTBP1 + GAPDH + ALB",
       x = "audit_score (median ± 95% CI)", y = NULL, color = "class") +
  theme_minimal(base_size = 10) +
  geom_vline(xintercept = c(0.3, 0.5), linetype = "dashed", color = "grey80")
ggsave(file.path(PKG, "data-raw", "mc_case_study_ci.png"), p2, width = 6, height = 4)

# Fig 3: uncertainty_width vs audit_score scatter
p3 <- ggplot(ref[!is.na(audit_score_median)],
              aes(x = audit_score_median, y = uncertainty_width)) +
  geom_point(aes(color = confidence_class), alpha = 0.4, size = 0.6) +
  scale_color_manual(values = c("stable_high" = "#1B5E20", "high_uncertain" = "#F57C00",
                                  "medium" = "#1976D2", "low" = "grey60",
                                  "excluded" = "black")) +
  labs(title = "Uncertainty width vs audit_score (high-score / high-uncertainty quadrant)",
       x = "audit_score (median)", y = "95% CI width", color = "class") +
  theme_minimal(base_size = 9) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 0.2, linetype = "dashed")
ggsave(file.path(PKG, "data-raw", "mc_uncertainty_scatter.png"), p3, width = 7, height = 5)

cat(sprintf("\nSaved: %s\n", out_path))
cat(sprintf("Figures: mc_top50_ci.png, mc_case_study_ci.png, mc_uncertainty_scatter.png\n"))
