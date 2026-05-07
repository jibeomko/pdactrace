#!/usr/bin/env Rscript
# data-raw/build_audit_scores_v2.R   (pdactrace v0.3.0, 3+2 framework)
#
# Attaches the locked 3-axis + 2-gate + 4-class output of
# compute_audit_score() to data/pdactrace_reference.rda. Replaces the
# previous 7-feature audit_score / audit_confidence_class columns
# in-place; mapping detail kept for Supplement Table S1.
#
# Run from package root:
#   Rscript data-raw/build_audit_scores_v2.R
suppressPackageStartupMessages({ library(data.table) })

PKG <- rprojroot::find_package_root_file()
devtools::load_all(PKG, quiet = TRUE)

ref_path <- file.path(PKG, "data", "pdactrace_reference.rda")
load(ref_path)
ref <- as.data.table(pdactrace_reference)

# ── 3+2 atlas score ─────────────────────────────────────────────
score_v3 <- compute_audit_score(NULL)
setnames(score_v3,
  old = c("evidence_strength", "biological_coherence",
          "translational_relevance", "leakage_gate", "heterogeneity_gate",
          "positive_score", "audit_score", "audit_class"),
  new = c("audit_evidence_strength", "audit_biological_coherence",
          "audit_translational_relevance", "audit_leakage_gate",
          "audit_heterogeneity_gate", "audit_positive_score_v3",
          "audit_score_v3", "audit_class"))

# ── Join + replace canonical audit_score with v3 value ─────────
new_cols <- setdiff(names(score_v3), "gene_symbol")
for (col in new_cols) if (col %in% names(ref)) ref[, (col) := NULL]
ref <- score_v3[ref, on = "gene_symbol"]

# Promote v3 to the canonical audit_score column (preserve column name
# downstream code already references) and drop the *_v3 alias.
ref[, audit_score := audit_score_v3]
ref[, audit_positive_score := audit_positive_score_v3]
ref[, c("audit_score_v3", "audit_positive_score_v3") := NULL]

# ── Refresh MC quantile columns from audit_score_mc_v1.csv ─────
mc_path <- file.path(PKG, "data-raw", "audit_score_mc_v1.csv")
if (file.exists(mc_path)) {
  mc <- fread(mc_path)
  mc_keep <- intersect(
    c("gene_symbol", "audit_score_median", "audit_score_lo95",
      "audit_score_hi95", "uncertainty_width", "rank_median",
      "rank_lo95", "rank_hi95", "confidence_class"),
    names(mc))
  mc <- mc[, ..mc_keep]
  setnames(mc,
    old = setdiff(names(mc), "gene_symbol"),
    new = c("audit_score_median", "audit_score_lo95",
            "audit_score_hi95", "audit_uncertainty_width",
            "audit_rank_median", "audit_rank_lo95",
            "audit_rank_hi95", "audit_confidence_class")[
              seq_along(setdiff(names(mc), "gene_symbol"))])
  for (col in setdiff(names(mc), "gene_symbol")) {
    if (col %in% names(ref)) ref[, (col) := NULL]
  }
  ref <- mc[ref, on = "gene_symbol"]
}

# ── Order columns: identifiers first, audit block grouped ──────
audit_block <- c(
  # 3-axis
  "audit_evidence_strength", "audit_biological_coherence",
  "audit_translational_relevance",
  # 2-gate
  "audit_leakage_gate", "audit_heterogeneity_gate",
  # final
  "audit_positive_score", "audit_score", "audit_class",
  # MC
  "audit_score_median", "audit_score_lo95", "audit_score_hi95",
  "audit_uncertainty_width", "audit_rank_median", "audit_rank_lo95",
  "audit_rank_hi95", "audit_confidence_class")
audit_block <- intersect(audit_block, names(ref))
other_cols <- setdiff(names(ref), audit_block)
setcolorder(ref, c(other_cols, audit_block))

# ── Write back ─────────────────────────────────────────────────
pdactrace_reference <- as.data.frame(ref, stringsAsFactors = FALSE)
save(pdactrace_reference, file = ref_path, compress = "xz")

# ── Verification ───────────────────────────────────────────────
cat("\n── 3+2 attach summary ──\n")
cat(sprintf("Rows: %d   Cols: %d\n", nrow(ref), ncol(ref)))
cat("Class distribution:\n")
print(table(ref$audit_class, useNA = "ifany"))
cat("\nCase verification:\n")
cases <- ref[gene_symbol %in% c("LTBP1", "GAPDH", "ALB", "THBS2", "LGALS3BP"),
             .(gene_symbol, audit_score, audit_class,
               audit_leakage_gate, audit_heterogeneity_gate)]
print(cases)
message("\nAttached 3+2 framework columns to ", ref_path)
