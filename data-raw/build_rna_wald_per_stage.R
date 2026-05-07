#!/usr/bin/env Rscript
# ===========================================================
# data-raw/build_rna_wald_per_stage.R
#
# Refit dds_lrt_4group.rds with nbinomWaldTest() to extract
# genuine per-contrast Wald padj (Early/Mid/Late vs Normal).
#
# Why: under DESeq2 LRT mode, results(name = "stage_group_*")
# returns the same omnibus LRT padj for every contrast. The
# atlas columns rna_padj_E/M/L are therefore LRT padj copies,
# not per-contrast Wald padj. v0.4.0 adds rna_wald_padj_E/M/L
# so users can inspect contrast-specific significance without
# refitting.
#
# Run from package root:
#   Rscript data-raw/build_rna_wald_per_stage.R
# ===========================================================
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
})

PROJ <- "/home/kjb9412/PDAC_biomarker"
RNA  <- file.path(PROJ, "analysis/transcriptomics/results/figure1")
OUT  <- file.path(PROJ, "data-raw", "rna_wald_per_stage_padj.csv")

cat("[1/4] Loading dds_lrt_4group.rds (Normal/Early/Mid/Late LRT) ...\n")
dds_lrt <- readRDS(file.path(RNA, "dds_lrt_4group.rds"))

cat("[2/4] Refitting same design (~ dataset + stage_group) with Wald ...\n")
# nbinomWaldTest() reuses dispersion + size factors already estimated
# during DESeq(test="LRT"); only the per-coefficient Wald p-values are
# (re)computed. This is the recommended way to add Wald output to a
# pre-fit LRT model without redoing dispersion estimation.
dds_wald <- nbinomWaldTest(dds_lrt, quiet = TRUE)

cat("[3/4] Extracting per-contrast Wald padj (vs Normal) ...\n")
resE <- results(dds_wald, name = "stage_group_Early_vs_Normal",
                test = "Wald")
resM <- results(dds_wald, name = "stage_group_Mid_vs_Normal",
                test = "Wald")
resL <- results(dds_wald, name = "stage_group_Late_vs_Normal",
                test = "Wald")

wald_dt <- data.table(
  gene_id          = rownames(dds_wald),
  rna_wald_padj_E  = resE$padj,
  rna_wald_padj_M  = resM$padj,
  rna_wald_padj_L  = resL$padj,
  rna_wald_pval_E  = resE$pvalue,
  rna_wald_pval_M  = resM$pvalue,
  rna_wald_pval_L  = resL$pvalue
)

cat(sprintf("  rows: %d genes\n", nrow(wald_dt)))
cat(sprintf("  Early Wald padj < 0.05: %d\n",
            sum(wald_dt$rna_wald_padj_E < 0.05, na.rm = TRUE)))
cat(sprintf("  Mid   Wald padj < 0.05: %d\n",
            sum(wald_dt$rna_wald_padj_M < 0.05, na.rm = TRUE)))
cat(sprintf("  Late  Wald padj < 0.05: %d\n",
            sum(wald_dt$rna_wald_padj_L < 0.05, na.rm = TRUE)))

cat("[4/4] Writing data-raw/rna_wald_per_stage_padj.csv ...\n")
fwrite(wald_dt, OUT)
cat(sprintf("  wrote: %s\n", OUT))
cat("Done.\n")
