# bundle_phase_csvs.R ----------------------------------------------
# Copies + trims + xz-compresses the manuscript-monorepo phase CSVs
# that build_reference.R / build_protein_betas.R consume, so the
# bundled atlas can be regenerated from the package alone (without
# the companion manuscript repo).
#
# Source: /home/kjb9412/PDAC_biomarker/analysis/manuscript/...
# Sink:   inst/extdata/phase*.csv.xz
#
# Re-run via:
#   PDAC_BASE_DIR=/home/kjb9412/PDAC_biomarker \
#     Rscript data-raw/bundle_phase_csvs.R

base_dir <- Sys.getenv("PDAC_BASE_DIR",
                        "/home/kjb9412/PDAC_biomarker")
ts_dir <- file.path(
  base_dir,
  "analysis/manuscript/tissue_to_serum_biomarker/results")
out_dir <- "inst/extdata"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_xz <- function(dt, name) {
  fp_csv <- file.path(out_dir, paste0(name, ".csv"))
  fp_xz  <- paste0(fp_csv, ".xz")
  data.table::fwrite(dt, fp_csv)
  if (file.exists(fp_xz)) unlink(fp_xz)
  status <- system2("xz", c("-9", "-e", "-f", shQuote(fp_csv)))
  if (status != 0L || !file.exists(fp_xz)) {
    stop("xz compression failed for ", fp_csv, call. = FALSE)
  }
  cat(sprintf("  %-55s %s\n", basename(fp_xz),
              format(structure(file.info(fp_xz)$size,
                               class = "object_size"),
                     units = "auto")))
}

cat("Trimming + bundling phase CSVs to inst/extdata/...\n")

# phase33 — RNA 12-template fit (largest, ~2.6 MB raw).  Keep only the
# columns build_reference.R selects.
p33 <- data.table::fread(
  file.path(ts_dir, "phase33_deseq2_coef_12template.csv"))
p33_keep <- c("gene_id", "gene_symbol", "lrt_padj",
               "beta_N", "beta_E", "beta_M", "beta_L",
               "padj_E", "padj_M", "padj_L",
               "lfcSE_E", "lfcSE_M", "lfcSE_L",
               "coef_pat", "coef_rho", "max_abs_beta")
write_xz(p33[, intersect(p33_keep, names(p33)), with = FALSE],
          "phase33_deseq2_coef_12template")

# phase34 — protein 12-template
p34 <- data.table::fread(
  file.path(ts_dir, "phase34_protein_pooled_12template.csv"))
write_xz(p34, "phase34_protein_pooled_12template")

# Smaller phase CSVs — copied verbatim
for (nm in c("phase2c_celltype_specificity",
             "phase42_pancreatitis_check",
             "phase77_strict_RNAprotConvergent_serum",
             "phase29_resectable_markers",
             "phase60_signalP_pipeline",
             "phase80_ltbp1_pancreatitis_predeclared_panels")) {
  src <- file.path(ts_dir, paste0(nm, ".csv"))
  if (file.exists(src)) {
    write_xz(data.table::fread(src), nm)
  } else {
    cat("  (skip — not found)", basename(src), "\n")
  }
}

cat("\nTotal inst/extdata size:\n")
print(structure(sum(file.info(list.files(out_dir, full.names = TRUE,
                                           pattern = "\\.csv\\.xz$"))$size),
                class = "object_size"), units = "auto")
