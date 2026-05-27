#!/usr/bin/env Rscript
# data-raw/attach_pareto_columns.R  (pdactrace v0.99.19, weight-free
# Pareto layer)
#
# Attaches the deterministic Pareto layer + MC stability columns
# produced by `compute_pareto_layers()` and
# `evaluate_pareto_stability()` to `data/pdactrace_reference.rda`.
# The frozen 3+2 audit_score / audit_class formula is NOT modified —
# this script appends 9 new columns next to the existing
# `audit_*` block.
#
# Run from package root:
#   Rscript data-raw/attach_pareto_columns.R
suppressPackageStartupMessages({ library(data.table) })

PKG <- rprojroot::find_package_root_file()
devtools::load_all(PKG, quiet = TRUE)

ref_path <- file.path(PKG, "data", "pdactrace_reference.rda")
load(ref_path)
ref <- as.data.table(pdactrace_reference)

cache_layers <- file.path(PKG, "data-raw", "pareto_layers_v0.99.19.csv")
cache_stab   <- file.path(PKG, "data-raw", "pareto_stability_v0.99.19.csv")

# ── 1. Deterministic Pareto layers (top-2000 by audit_score) ─────
det <- compute_pareto_layers(top_n = 2000L)
fwrite(det, cache_layers)
cat(sprintf("Wrote %s (n=%d, layer-1 count=%d)\n",
            basename(cache_layers), nrow(det),
            sum(det$pareto_layer == 1L, na.rm = TRUE)))

# ── 2. Monte Carlo Pareto stability (1000 draws) ────────────────
if (file.exists(cache_stab)) {
  cat("Reusing cached MC stability from ", basename(cache_stab), "\n")
  mc <- fread(cache_stab)
} else {
  cat("Computing MC Pareto stability (n_draws = 1000) ...\n")
  t0 <- Sys.time()
  mc <- evaluate_pareto_stability(n_draws = 1000L, top_n = 2000L,
                                   seed = 20250506L)
  fwrite(mc, cache_stab)
  cat(sprintf("Wrote %s (elapsed %.1fs)\n",
              basename(cache_stab),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# ── 3. Join into atlas ────────────────────────────────────────
join_cols <- c("pareto_layer", "pareto_rank", "crowding_distance",
               "pareto_excluded_by_gate")
for (col in join_cols) if (col %in% names(ref)) ref[, (col) := NULL]
ref <- det[ref, on = "gene_symbol"]

mc_cols <- setdiff(names(mc), "gene_symbol")
for (col in mc_cols) if (col %in% names(ref)) ref[, (col) := NULL]
ref <- mc[ref, on = "gene_symbol"]

# ── 4. Order columns: keep audit block contiguous, append pareto ──
pareto_block <- c(
  "pareto_layer", "pareto_rank", "crowding_distance",
  "pareto_excluded_by_gate",
  "pareto_layer_median", "pareto_stability_top1",
  "pareto_layer_lo95", "pareto_layer_hi95",
  "pareto_top10_pct_stability")
pareto_block <- intersect(pareto_block, names(ref))
audit_anchor <- "audit_confidence_class"
if (audit_anchor %in% names(ref)) {
  ord <- names(ref)
  ord <- setdiff(ord, pareto_block)
  ins <- which(ord == audit_anchor)
  ord <- append(ord, pareto_block, after = ins)
  setcolorder(ref, ord)
}

# ── 5. Write back ──────────────────────────────────────────────
pdactrace_reference <- as.data.frame(ref, stringsAsFactors = FALSE)
save(pdactrace_reference, file = ref_path, compress = "xz")

# ── 6. Verification ────────────────────────────────────────────
cat("\n── Pareto attach summary ──\n")
cat(sprintf("Rows: %d   Cols: %d\n", nrow(ref), ncol(ref)))
cat("Pareto layer distribution (top-2000 pool):\n")
print(table(ref$pareto_layer, useNA = "ifany"))
cat("\nGate-excluded count: ",
    sum(ref$pareto_excluded_by_gate, na.rm = TRUE), "\n")
cat("\nCase verification (top-10 by pareto_rank):\n")
top_pareto <- ref[!is.na(pareto_rank)][order(pareto_rank)][1:10,
  .(gene_symbol, audit_score, pareto_layer, pareto_rank,
    pareto_stability_top1)]
print(top_pareto)
message("\nAttached Pareto columns to ", ref_path)
