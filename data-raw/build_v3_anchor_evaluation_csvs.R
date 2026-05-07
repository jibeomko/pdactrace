#!/usr/bin/env Rscript
# ===========================================================
# data-raw/build_v3_anchor_evaluation_csvs.R
#
# Regenerates the four v0.3.0 anchor-evaluation CSVs that
# manuscript numbers reference, but using whatever atlas
# (audit_score column) is currently bundled. Run after every
# atlas rebuild that changes audit_score values.
#
# Output:
#   data-raw/anchor_enrichment_3plus2_results.csv
#   data-raw/loo_anchor_enrichment_v3_results.csv
#   data-raw/bootstrap_anchor_enrichment_v3_results.csv
#   data-raw/top100_audit_v3_deterministic.csv
# ===========================================================
suppressPackageStartupMessages({
  library(data.table)
})

PKG <- rprojroot::find_package_root_file()
devtools::load_all(PKG, quiet = TRUE)

set.seed(42)

# ── Atlas + anchors ────────────────────────────────────────
load(file.path(PKG, "data", "pdactrace_reference.rda"))
load(file.path(PKG, "data", "pdactrace_external_anchors.rda"))
ref <- as.data.table(pdactrace_reference)
anc <- as.data.table(pdactrace_external_anchors)

primary_anc <- intersect(unique(anc[include_primary_eval == TRUE, gene]),
                            ref$gene_symbol)
secondary_anc <- intersect(unique(anc[include_secondary_eval == TRUE, gene]),
                              ref$gene_symbol)
cat(sprintf("Atlas: %d genes  Anchors: P=%d  S=%d (in-atlas)\n",
            nrow(ref), length(primary_anc), length(secondary_anc)))

# ── 1. anchor_enrichment_3plus2_results.csv ───────────────
# 15 rows = 5 top_n × {primary-deterministic, secondary-deterministic,
#                      secondary-mc_median}
top_n_vec <- c(50, 100, 200, 500, 1000)
n_atlas <- nrow(ref)

build_row <- function(score_vec, anchors, top_n, n_decl, tier, tier_key,
                       ranking) {
  ord <- order(-score_vec)
  top_set <- ref$gene_symbol[ord[seq_len(top_n)]]
  n_in <- length(anchors)
  hits <- sum(top_set %in% anchors)
  expect <- top_n * n_in / n_atlas
  fold <- hits / max(expect, 1e-12)
  pval <- phyper(hits - 1, n_in, n_atlas - n_in, top_n, lower.tail = FALSE)
  data.table(tier = tier, tier_key = tier_key, top_n = top_n,
              n_anchor_declared = n_decl, n_anchor_in_atlas = n_in,
              hits = hits, expect = expect, fold = fold, pval = pval,
              score_col = "audit_score", ranking = ranking)
}

all_res <- list()
n_decl_p <- length(unique(anc[include_primary_eval == TRUE, gene]))
n_decl_s <- length(unique(anc[include_secondary_eval == TRUE, gene]))

for (tn in top_n_vec) {
  all_res[[length(all_res) + 1]] <- build_row(
    ref$audit_score, primary_anc, tn, n_decl_p,
    "PRIMARY (T1_validated, direct)", "primary", "deterministic")
}
for (tn in top_n_vec) {
  all_res[[length(all_res) + 1]] <- build_row(
    ref$audit_score, secondary_anc, tn, n_decl_s,
    "SECONDARY (T1+T2)", "secondary", "deterministic")
}
for (tn in top_n_vec) {
  all_res[[length(all_res) + 1]] <- build_row(
    ref$audit_score_median, secondary_anc, tn, n_decl_s,
    "SECONDARY (T1+T2)", "secondary", "mc_median")
}
res_3plus2 <- rbindlist(all_res)
fwrite(res_3plus2,
        file.path(PKG, "data-raw", "anchor_enrichment_3plus2_results.csv"))
cat(sprintf("Wrote anchor_enrichment_3plus2_results.csv (%d rows)\n",
            nrow(res_3plus2)))

# ── 2. loo_anchor_enrichment_v3_results.csv ───────────────
# Drop one secondary anchor at a time, recompute top-100 enrichment
ord <- order(-ref$audit_score)
top100 <- ref$gene_symbol[ord[seq_len(100)]]

loo_rows <- list()
for (drop in secondary_anc) {
  k_remain <- secondary_anc[secondary_anc != drop]
  hits <- sum(top100 %in% k_remain)
  K <- length(k_remain)
  expect <- 100 * K / n_atlas
  fold <- hits / max(expect, 1e-12)
  pval <- phyper(hits - 1, K, n_atlas - K, 100, lower.tail = FALSE)
  loo_rows[[length(loo_rows) + 1]] <- data.table(
    dropped = drop, K_remaining = K, hits = hits,
    expect = expect, fold = fold, pval = pval)
}
loo_dt <- rbindlist(loo_rows)
fwrite(loo_dt,
        file.path(PKG, "data-raw", "loo_anchor_enrichment_v3_results.csv"))
cat(sprintf("Wrote loo_anchor_enrichment_v3_results.csv (%d rows)\n",
            nrow(loo_dt)))

# ── 3. bootstrap_anchor_enrichment_v3_results.csv ─────────
B <- 1000
boot_rows <- vector("list", B)
for (b in seq_len(B)) {
  resampled <- sample(secondary_anc, length(secondary_anc), replace = TRUE)
  uniq <- unique(resampled)
  K <- length(uniq)
  hits <- sum(top100 %in% uniq)
  expect <- 100 * K / n_atlas
  fold <- hits / max(expect, 1e-12)
  boot_rows[[b]] <- data.table(boot_iter = b, n_unique_anchors = K,
                                  hits = hits, fold = fold)
}
boot_dt <- rbindlist(boot_rows)
fwrite(boot_dt,
        file.path(PKG, "data-raw", "bootstrap_anchor_enrichment_v3_results.csv"))
cat(sprintf("Wrote bootstrap_anchor_enrichment_v3_results.csv (%d rows)\n",
            nrow(boot_dt)))
cat(sprintf("  median fold = %.2f, 95%% CI = [%.2f, %.2f]\n",
            median(boot_dt$fold),
            quantile(boot_dt$fold, 0.025),
            quantile(boot_dt$fold, 0.975)))

# ── 4. top100_audit_v3_deterministic.csv ──────────────────
top100_dt <- ref[ord[seq_len(100)],
                   .(gene_symbol,
                     evidence_strength       = audit_evidence_strength,
                     biological_coherence    = audit_biological_coherence,
                     translational_relevance = audit_translational_relevance,
                     leakage_gate            = audit_leakage_gate,
                     heterogeneity_gate      = audit_heterogeneity_gate,
                     positive_score          = audit_positive_score,
                     audit_score             = audit_score,
                     audit_class             = audit_class,
                     rank                    = seq_len(.N),
                     audit_score_median      = audit_score_median,
                     audit_score_lo95        = audit_score_lo95,
                     audit_score_hi95        = audit_score_hi95,
                     uncertainty_width       = audit_uncertainty_width,
                     confidence_class        = audit_confidence_class)]
top100_dt[, anchor_secondary := gene_symbol %in% secondary_anc]
top100_dt[, anchor_primary   := gene_symbol %in% primary_anc]
fwrite(top100_dt,
        file.path(PKG, "data-raw", "top100_audit_v3_deterministic.csv"))
cat(sprintf("Wrote top100_audit_v3_deterministic.csv (%d rows)\n",
            nrow(top100_dt)))
cat(sprintf("  Anchor recovery in top-100: P=%d, S=%d\n",
            sum(top100_dt$anchor_primary),
            sum(top100_dt$anchor_secondary)))

cat("Done.\n")
