#!/usr/bin/env Rscript
# ===========================================================
# data-raw/attach_wald_padj_columns.R
#
# Inject rna_wald_padj_E/M/L into pdactrace_reference.rda
# without rebuilding the whole atlas (avoids losing the
# meta_analysis + tier + audit columns that are populated by
# downstream scripts after build_reference.R).
#
# Inputs:
#   data-raw/rna_wald_per_stage_padj.csv (from
#     build_rna_wald_per_stage.R)
#   data-raw/anchor_per_gene_audit_results.csv  (gene_id <-> gene_symbol)
#   data/pdactrace_reference.rda                (existing atlas)
#
# Output: overwrites data/pdactrace_reference.rda with 3 added cols.
# ===========================================================
suppressPackageStartupMessages({
  library(data.table)
})

PROJ <- "/home/kjb9412/PDAC_biomarker"
PKG  <- rprojroot::find_package_root_file()

cat("[1/4] Loading current atlas ...\n")
load(file.path(PKG, "data", "pdactrace_reference.rda"))
ref <- as.data.table(pdactrace_reference)
cat(sprintf("  ref rows=%d  cols=%d\n", nrow(ref), ncol(ref)))

cat("[2/4] Loading per-contrast Wald padj table ...\n")
wald_dt <- fread(file.path(PROJ, "data-raw",
                              "rna_wald_per_stage_padj.csv"))
wald_dt[, ens_join := sub("\\..*$", "", gene_id)]
cat(sprintf("  wald rows=%d\n", nrow(wald_dt)))

cat("[3/4] Joining Wald padj on ensembl_id ...\n")
ref[, ens_join := sub("\\..*$", "", ensembl_id)]
ref <- merge(ref, wald_dt[, .(ens_join, rna_wald_padj_E,
                                 rna_wald_padj_M, rna_wald_padj_L)],
              by = "ens_join", all.x = TRUE, sort = FALSE)
ref[, ens_join := NULL]
cat(sprintf("  matched: %d / %d\n",
            sum(!is.na(ref$rna_wald_padj_E)), nrow(ref)))

# Reorder so the new columns sit next to rna_padj_E/M/L
ord <- names(ref)
i_lfcSE_E <- which(ord == "rna_lfcSE_E")
new_cols <- c("rna_wald_padj_E", "rna_wald_padj_M", "rna_wald_padj_L")
ord_no_new <- setdiff(ord, new_cols)
i_lfcSE_E <- which(ord_no_new == "rna_lfcSE_E")
ord_final <- c(ord_no_new[seq_len(i_lfcSE_E - 1)],
                new_cols,
                ord_no_new[i_lfcSE_E:length(ord_no_new)])
setcolorder(ref, ord_final)
setkey(ref, gene_symbol)

cat("[4/4] Saving updated atlas ...\n")
pdactrace_reference <- ref
save(pdactrace_reference,
     file = file.path(PKG, "data", "pdactrace_reference.rda"),
     compress = "xz")
cat(sprintf("  saved with %d cols (was 111).\n", ncol(ref)))
cat("Done.\n")
