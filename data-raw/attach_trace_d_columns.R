#!/usr/bin/env Rscript
# data-raw/attach_trace_d_columns.R  (pdactrace v0.99.19, Algorithm 1:
# TRACE-D tissue-to-serum directional evidence concordance)
#
# Attaches the deterministic TRACE-D class + confidence + pancreatitis
# specificity columns produced by `compute_trace_d()` to
# `data/pdactrace_reference.rda`. The frozen 3+2 audit_score formula
# AND the legacy `translation_class` column are NOT modified — TRACE-D
# coexists with `translation_class` for back-compat (concordance table
# is written for the manuscript suppl).
#
# Run from package root (after attach_pareto_columns.R has populated
# pareto_*):
#   Rscript data-raw/attach_trace_d_columns.R
suppressPackageStartupMessages({ library(data.table) })

PKG <- rprojroot::find_package_root_file()
devtools::load_all(PKG, quiet = TRUE)

ref_path <- file.path(PKG, "data", "pdactrace_reference.rda")
load(ref_path)
ref <- as.data.table(pdactrace_reference)

cache_tracd       <- file.path(PKG, "data-raw",
                                "tracd_columns_v0.99.19.csv")
cache_concordance <- file.path(PKG, "data-raw",
                                "tracd_concordance_v0.99.19.csv")

# ── 1. Compute TRACE-D ────────────────────────────────────────
t0 <- Sys.time()
out <- compute_trace_d()
fwrite(out, cache_tracd)
cat(sprintf("Wrote %s (elapsed %.2fs, n=%d)\n",
            basename(cache_tracd),
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            nrow(out)))

# ── 2. Concordance with legacy translation_class ─────────────
legacy <- ref[!is.na(translation_class),
              .(gene_symbol, translation_class)]
conc <- merge(legacy, out[, .(gene_symbol, tracd_class)],
              by = "gene_symbol", all.x = TRUE)
conc[, agree := tracd_class == translation_class]
fwrite(conc, cache_concordance)
n_legacy <- nrow(conc)
n_agree <- sum(conc$agree, na.rm = TRUE)
cat(sprintf("Concordance with legacy translation_class: %d/%d = %.1f%%\n",
            n_agree, n_legacy,
            ifelse(n_legacy > 0, 100 * n_agree / n_legacy, NA)))

# ── 3. Join into atlas ────────────────────────────────────────
tracd_cols <- setdiff(names(out), "gene_symbol")
for (col in tracd_cols) if (col %in% names(ref)) ref[, (col) := NULL]
ref <- out[ref, on = "gene_symbol"]

# ── 4. Order columns: keep tracd block after pareto block ─────
tracd_block <- c(
  "tracd_tissue_dir", "tracd_serum_dir",
  "tracd_class", "tracd_confidence",
  "tracd_pancreatitis_overlap_score",
  "tracd_pancreatitis_specificity",
  "tracd_tissue_weight", "tracd_decision_path")
tracd_block <- intersect(tracd_block, names(ref))
pareto_anchor <- "pareto_top10_pct_stability"
if (pareto_anchor %in% names(ref)) {
  ord <- setdiff(names(ref), tracd_block)
  ins <- which(ord == pareto_anchor)
  ord <- append(ord, tracd_block, after = ins)
  setcolorder(ref, ord)
}

# ── 5. Write back ──────────────────────────────────────────────
pdactrace_reference <- as.data.frame(ref, stringsAsFactors = FALSE)
save(pdactrace_reference, file = ref_path, compress = "xz")

# ── 6. Verification ────────────────────────────────────────────
cat("\n── TRACE-D attach summary ──\n")
cat(sprintf("Rows: %d   Cols: %d\n", nrow(ref), ncol(ref)))
cat("\ntracd_class distribution:\n")
print(table(ref$tracd_class, useNA = "ifany"))
cat("\ntracd_pancreatitis_specificity distribution:\n")
print(table(ref$tracd_pancreatitis_specificity, useNA = "ifany"))

# Domain sanity check
sanity_genes <- c("LGALS3BP", "LTBP1", "GAPDH", "THBS2", "ALB")
cat("\nDomain sanity check:\n")
print(ref[gene_symbol %in% sanity_genes,
          .(gene_symbol, translation_class, tracd_class,
            tracd_confidence, tracd_pancreatitis_specificity,
            tracd_decision_path)])

message("\nAttached TRACE-D columns to ", ref_path)
