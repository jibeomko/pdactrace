#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════
# data-raw/build_meta_analysis.R  (pdactrace v0.2.0 T1 stat-track)
#
# Random-effects meta-analysis of per-cohort DESeq2 stage
# coefficients, producing per-gene β_meta, SE, τ², I² for two
# primary contrasts:
#   * "Mid_vs_Early"  — k=4 cohorts (TCGA, CPTAC, GSE224564, GSE79668)
#   * "Late_vs_Early" — k=4 cohorts
# Plus auxiliary 2-cohort meta:
#   * "Early_vs_Normal"  — k=2 (TCGA, CPTAC only)
#   * "Mid_vs_Normal"    — k=2 (TCGA, CPTAC)
#   * "Late_vs_Normal"   — k=2 (TCGA, CPTAC)
#
# Method: metafor::rma() with REML estimator (modern default).
# Output: data/meta_analysis.rda inside the R package.
# ═══════════════════════════════════════════════════════════════
suppressPackageStartupMessages({
  library(metafor)
  library(data.table)
})

ROOT     <- "/home/kjb9412/PDAC_biomarker"
RNA_DIR  <- file.path(ROOT, "analysis/transcriptomics/results/figure1")
PKG      <- rprojroot::find_package_root_file()

cat("=== Loading per-cohort betas ===\n")
b <- fread(file.path(RNA_DIR, "per_cohort_betas_long.csv"))
cat(sprintf("  rows: %d\n", nrow(b)))
cat(sprintf("  unique genes: %d\n", uniqueN(b$gene_symbol)))

run_meta <- function(beta_vec, se_vec) {
  ok <- !is.na(beta_vec) & !is.na(se_vec) & se_vec > 0
  if (sum(ok) < 2) {
    return(list(beta = NA_real_, se = NA_real_, pval = NA_real_,
                tau2 = NA_real_, I2 = NA_real_, k = sum(ok)))
  }
  fit <- tryCatch(
    rma(yi = beta_vec[ok], sei = se_vec[ok], method = "REML",
        control = list(maxiter = 200), verbose = FALSE),
    error = function(e) NULL,
    warning = function(w) {
      tryCatch(
        rma(yi = beta_vec[ok], sei = se_vec[ok], method = "DL"),
        error = function(e) NULL)
    })
  if (is.null(fit)) {
    return(list(beta = NA_real_, se = NA_real_, pval = NA_real_,
                tau2 = NA_real_, I2 = NA_real_, k = sum(ok)))
  }
  list(beta = as.numeric(fit$beta),
       se   = as.numeric(fit$se),
       pval = as.numeric(fit$pval),
       tau2 = fit$tau2,
       I2   = fit$I2,
       k    = sum(ok))
}

contrasts_to_meta <- list(
  Mid_vs_Early    = "Mid_vs_Early",     # k=4 (all cohorts)
  Late_vs_Early   = "Late_vs_Early",    # k=4 (all cohorts)
  Normal_vs_Early = "Normal_vs_Early")  # k=2 (TCGA + CPTAC only)

cat("\n=== Running meta-analysis per gene ===\n")
genes <- unique(b$gene_symbol)
genes <- genes[!is.na(genes) & genes != ""]
out <- data.table(gene_symbol = genes)
for (ct_name in names(contrasts_to_meta)) {
  ct_target <- contrasts_to_meta[[ct_name]]
  cat(sprintf("  contrast: %s ...\n", ct_name))
  bsub <- b[contrast == ct_target]
  meta_res <- bsub[, run_meta(beta, se), by = gene_symbol]
  data.table::setnames(meta_res,
    c("beta","se","pval","tau2","I2","k"),
    paste0("meta_", ct_name, c("_beta","_se","_pval","_tau2","_I2","_k")))
  out <- merge(out, meta_res, by = "gene_symbol", all.x = TRUE)
}

# Add BH-padj for primary contrasts
for (ct_name in names(contrasts_to_meta)) {
  pcol <- paste0("meta_", ct_name, "_pval")
  qcol <- paste0("meta_", ct_name, "_padj")
  out[, (qcol) := p.adjust(get(pcol), method = "BH")]
}

# Cohort-divergence flag: high I² in either Mid_vs_Early or Late_vs_Early
out[, meta_cohort_divergent := pmax(meta_Mid_vs_Early_I2,
                                       meta_Late_vs_Early_I2,
                                       na.rm = TRUE) >= 50]

cat("\n=== Sanity checks (LTBP1 / SERPINA1 / ACTB / GAPDH) ===\n")
for (g in c("LTBP1", "SERPINA1", "ACTB", "GAPDH", "CDH13", "SPARC")) {
  cat("\n--- ", g, " ---\n", sep = "")
  r <- out[gene_symbol == g]
  if (nrow(r) == 0) { cat("  not found\n"); next }
  for (ct_name in names(contrasts_to_meta)) {
    cat(sprintf("  %s: β=%.2f, SE=%.2f, p=%.3g, I²=%.0f%%, k=%d\n",
        ct_name,
        r[[paste0("meta_", ct_name, "_beta")]],
        r[[paste0("meta_", ct_name, "_se")]],
        r[[paste0("meta_", ct_name, "_pval")]],
        r[[paste0("meta_", ct_name, "_I2")]],
        r[[paste0("meta_", ct_name, "_k")]]))
  }
}

cat("\n=== Atlas-wide summary ===\n")
cat(sprintf("  Total genes with meta: %d\n", nrow(out)))
cat(sprintf("  cohort_divergent (I² >= 50%%): %d (%.0f%%)\n",
            sum(out$meta_cohort_divergent, na.rm = TRUE),
            100 * sum(out$meta_cohort_divergent, na.rm = TRUE) / nrow(out)))
cat(sprintf("  Mid_vs_Early p < 0.05: %d (%.0f%%)\n",
            sum(out$meta_Mid_vs_Early_pval < 0.05, na.rm = TRUE),
            100 * sum(out$meta_Mid_vs_Early_pval < 0.05, na.rm = TRUE)/nrow(out)))
cat(sprintf("  Mid_vs_Early padj < 0.05: %d (%.0f%%)\n",
            sum(out$meta_Mid_vs_Early_padj < 0.05, na.rm = TRUE),
            100 * sum(out$meta_Mid_vs_Early_padj < 0.05, na.rm = TRUE)/nrow(out)))

meta_analysis <- out
out_path <- file.path(PKG, "data", "meta_analysis.rda")
save(meta_analysis, file = out_path, compress = "xz")
cat(sprintf("\nSaved: %s\n", out_path))
cat(sprintf("Cols: %d (incl. 5 contrasts × 6 stats + cohort_divergent flag)\n",
            ncol(meta_analysis)))
