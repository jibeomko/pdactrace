#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════
# 3-tier anchor enrichment evaluation
#
# Scoring rule: FROZEN at v1 freeze (prototype/v1 weights).
# Evaluation universe: 3 tiers (provenance-audited).
#
#   Primary  (include_primary_eval==TRUE):    T1_validated direct mapping
#   Secondary (include_secondary_eval==TRUE):  T1 + T2_literature_db
#   Exploratory (all):                         T1 + T2 + cyst-fluid/multi-omic
#
# This script does NOT modify weights. It only adds breadth.
# ═══════════════════════════════════════════════════════════════
suppressPackageStartupMessages({ library(data.table) })
PKG <- rprojroot::find_package_root_file()

# ── Frozen audit_score (prototype) ──────────────────────
score_dt <- fread(file.path(PKG, "data-raw", "audit_score_atlas_legacy_7feature.csv"))
n_atlas  <- nrow(score_dt)
score_dt <- score_dt[order(-audit_score)]
score_dt[, rank_overall := seq_len(.N)]
# DO NOT setkey — preserves audit_score descending order (bug fix)
score_lookup <- setNames(score_dt$audit_score, score_dt$gene_symbol)
rank_lookup  <- setNames(score_dt$rank_overall, score_dt$gene_symbol)

# ── Anchor provenance v2 ──────────────────────────────────────
anc <- fread(file.path(PKG, "data-raw", "external_positives_anchors_v2.csv"))

primary_anc <- anc[include_primary_eval == TRUE, gene]
secondary_anc <- anc[include_secondary_eval == TRUE, gene]
exploratory_anc <- anc$gene

cat("================================================\n")
cat("=== Anchor v2 provenance breakdown ===\n")
cat("================================================\n")
print(anc[, .N, by = .(evidence_tier)])
cat(sprintf("Primary eval:    %d anchors\n", length(primary_anc)))
cat(sprintf("Secondary eval:  %d anchors\n", length(secondary_anc)))
cat(sprintf("Exploratory:     %d anchors\n", length(exploratory_anc)))

# ── Evaluation function ───────────────────────────────────────
eval_tier <- function(anchors, label) {
  in_atlas <- intersect(anchors, score_dt$gene_symbol)
  n_in <- length(in_atlas)
  cat(sprintf("\n--- %s (n=%d declared, %d in atlas) ---\n",
              label, length(anchors), n_in))
  results <- list()
  for (top_n in c(50, 100, 200, 500, 1000)) {
    top_set <- score_dt[1:top_n, gene_symbol]
    hits   <- sum(top_set %in% in_atlas)
    expect <- top_n * n_in / n_atlas
    fold   <- hits / max(expect, 1e-12)
    pval   <- phyper(hits - 1, n_in, n_atlas - n_in, top_n, lower.tail = FALSE)
    cat(sprintf("  Top %4d : %d hits  (expect %.2f, fold=%.1fx, p=%.2e)\n",
                top_n, hits, expect, fold, pval))
    results[[as.character(top_n)]] <- data.table(
      tier = label, top_n = top_n, n_anchor_in_atlas = n_in,
      hits = hits, expect = round(expect, 2),
      fold = round(fold, 1), pval = pval)
  }
  rbindlist(results)
}

cat("\n================================================\n")
cat("=== 3-tier enrichment evaluation (frozen rule) ===\n")
cat("================================================\n")
res_primary  <- eval_tier(primary_anc, "PRIMARY (T1_validated, direct)")
res_secondary <- eval_tier(secondary_anc, "SECONDARY (T1+T2)")
res_exploratory <- eval_tier(exploratory_anc, "EXPLORATORY (all)")

# ── Per-anchor rank table ─────────────────────────────────────
cat("\n================================================\n")
cat("=== Per-anchor rank (sorted by audit_score) ===\n")
cat("================================================\n")
anc_atlas <- anc[gene %in% score_dt$gene_symbol]
anc_atlas[, audit_score := score_lookup[anc_atlas$gene]]
anc_atlas[, rank := rank_lookup[anc_atlas$gene]]
print(anc_atlas[order(rank),
                .(gene, evidence_tier, assay_context, early_detection,
                  audit_score = round(audit_score, 3), rank,
                  primary = include_primary_eval, sec = include_secondary_eval)])

# ── Summary table ─────────────────────────────────────────────
all_res <- rbind(res_primary, res_secondary, res_exploratory)
cat("\n================================================\n")
cat("=== Summary table (compare 3 tiers) ===\n")
cat("================================================\n")
print(all_res)

fwrite(all_res, file.path(PKG, "data-raw", "anchor_enrichment_3tier_results.csv"))

# ── Per-anchor table for manuscript supplement ────────────────
out_path <- file.path(PKG, "data-raw", "anchor_per_gene_audit_results.csv")
fwrite(anc_atlas[, .(gene, source, source_type, evidence_tier,
                      assay_context, early_detection, gene_level_mapping,
                      include_primary_eval, include_secondary_eval,
                      audit_score, rank, in_atlas = TRUE)],
       out_path)

# Anchors not in atlas
not_in <- anc[!gene %in% score_dt$gene_symbol]
if (nrow(not_in) > 0) {
  cat(sprintf("\n=== Anchors NOT in atlas (declared but not measurable) ===\n"))
  print(not_in[, .(gene, evidence_tier, assay_context, note)])
  fwrite(rbind(
    fread(out_path),
    not_in[, .(gene, source, source_type, evidence_tier, assay_context,
               early_detection, gene_level_mapping,
               include_primary_eval, include_secondary_eval,
               audit_score = NA_real_, rank = NA_integer_, in_atlas = FALSE)]),
    out_path)
}

cat(sprintf("\nSaved: %s\n", file.path(PKG, "data-raw", "anchor_enrichment_3tier_results.csv")))
cat(sprintf("Saved: %s\n", out_path))
