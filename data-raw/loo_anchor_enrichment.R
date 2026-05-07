#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════
# Leave-one-anchor-out enrichment robustness
#
# Frozen scoring rule unchanged. Per-iteration: drop one anchor,
# recompute top-N enrichment using remaining anchors. Reports:
#   - hit_count, fold, p   distribution across LOO iterations
#   - whether the result is dominated by any single anchor
#   - bootstrap (resample anchors with replacement) for completeness
# ═══════════════════════════════════════════════════════════════
suppressPackageStartupMessages({ library(data.table) })
PKG <- rprojroot::find_package_root_file()
set.seed(42)

# Frozen audit_score
score_dt <- fread(file.path(PKG, "data-raw",
                  "audit_score_atlas_legacy_7feature.csv"))
n_atlas <- nrow(score_dt)
score_dt <- score_dt[order(-audit_score)]
gene_ranked <- score_dt$gene_symbol

# Anchors
anc <- fread(file.path(PKG, "data-raw",
              "external_positives_anchors_v2.csv"))
primary_anc   <- intersect(anc[include_primary_eval == TRUE, gene],
                            gene_ranked)
secondary_anc <- intersect(anc[include_secondary_eval == TRUE, gene],
                            gene_ranked)

eval_set <- function(anchors, top_n_vec = c(50, 100, 200, 500)) {
  res <- list()
  for (top_n in top_n_vec) {
    top_set <- gene_ranked[1:top_n]
    n_in <- length(anchors)
    hits <- sum(top_set %in% anchors)
    expect <- top_n * n_in / n_atlas
    fold <- hits / max(expect, 1e-12)
    pval <- phyper(hits - 1, n_in, n_atlas - n_in, top_n,
                    lower.tail = FALSE)
    res[[as.character(top_n)]] <- data.table(top_n = top_n, n_anchors = n_in,
                                                hits = hits, fold = fold,
                                                pval = pval)
  }
  rbindlist(res)
}

# ── Reference (full set) ───────────────────────────────────────
cat("=== Reference: all anchors ===\n")
ref_secondary <- eval_set(secondary_anc)
print(ref_secondary)

# ── LOO across secondary set ───────────────────────────────────
cat("\n=== LOO secondary tier (n =", length(secondary_anc), ") ===\n")
loo_res <- list()
for (drop_gene in secondary_anc) {
  remaining <- setdiff(secondary_anc, drop_gene)
  e <- eval_set(remaining)
  e[, dropped := drop_gene]
  loo_res[[drop_gene]] <- e
}
loo_dt <- rbindlist(loo_res)

cat("\n--- Top 100 LOO distribution ---\n")
top100_loo <- loo_dt[top_n == 100]
cat(sprintf("Reference fold: %.1fx\n",
    ref_secondary[top_n == 100, fold]))
cat(sprintf("LOO fold:        median=%.1fx, IQR=[%.1f, %.1f], range=[%.1f, %.1f]\n",
    median(top100_loo$fold), quantile(top100_loo$fold, 0.25),
    quantile(top100_loo$fold, 0.75),
    min(top100_loo$fold), max(top100_loo$fold)))
cat(sprintf("LOO hits:        median=%d, range=[%d, %d]\n",
    median(top100_loo$hits), min(top100_loo$hits), max(top100_loo$hits)))
cat(sprintf("LOO p-values:    median=%.2e, max=%.2e\n",
    median(top100_loo$pval), max(top100_loo$pval)))

cat("\n--- Most influential anchors (drops causing largest fold change) ---\n")
top100_loo[, fold_change := fold - ref_secondary[top_n == 100, fold]]
print(top100_loo[order(fold_change)][, .(dropped, hits, fold, pval, fold_change)])

# ── Bootstrap (resample anchors with replacement, k = original size) ──
cat("\n=== Bootstrap (1000 resamples, k = secondary size) ===\n")
B <- 1000
boot_fold <- numeric(B)
boot_hits <- integer(B)
for (b in 1:B) {
  resampled <- sample(secondary_anc, length(secondary_anc), replace = TRUE)
  unique_resampled <- unique(resampled)
  e <- eval_set(unique_resampled, top_n_vec = 100)
  boot_fold[b] <- e$fold
  boot_hits[b] <- e$hits
}
cat(sprintf("Bootstrap top-100 fold: median=%.1fx, 95%% CI=[%.1f, %.1f]\n",
    median(boot_fold),
    quantile(boot_fold, 0.025), quantile(boot_fold, 0.975)))
cat(sprintf("Bootstrap top-100 hits: median=%d, 95%% CI=[%d, %d]\n",
    median(boot_hits),
    quantile(boot_hits, 0.025), quantile(boot_hits, 0.975)))

# ── Save ──────────────────────────────────────────────────────
fwrite(loo_dt, file.path(PKG, "data-raw",
        "loo_anchor_enrichment_results.csv"))
fwrite(data.table(boot_iter = 1:B, top100_fold = boot_fold,
                   top100_hits = boot_hits),
        file.path(PKG, "data-raw",
                  "bootstrap_anchor_enrichment_results.csv"))

cat("\nSaved: data-raw/loo_anchor_enrichment_results.csv\n")
cat("Saved: data-raw/bootstrap_anchor_enrichment_results.csv\n")
