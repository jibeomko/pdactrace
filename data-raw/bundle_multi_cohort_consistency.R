# bundle_multi_cohort_consistency.R --------------------------------
# Bundles multi_cohort_consistency.csv (Stouffer + per-cohort
# trend/monotonic, ~2.84 MB raw) into inst/extdata so that
# data-raw/build_reference.R can rebuild data/pdactrace_reference.rda
# without depending on the manuscript-monorepo path
# analysis/transcriptomics/results/figure1/.
#
# Re-run via:
#   PDAC_BASE_DIR=/home/kjb9412/PDAC_biomarker \
#     Rscript data-raw/bundle_multi_cohort_consistency.R

base_dir <- Sys.getenv("PDAC_BASE_DIR",
                        "/home/kjb9412/PDAC_biomarker")
src <- file.path(base_dir,
                  "analysis/transcriptomics/results/figure1",
                  "multi_cohort_consistency.csv")
out_dir <- "inst/extdata"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(src)) {
  stop("Source file not found:\n  ", src,
       "\nThis bundler must be run on a checkout that contains the ",
       "manuscript-monorepo analysis/ tree.", call. = FALSE)
}

mcc <- data.table::fread(src)
cat(sprintf("source: %s  (%d rows x %d cols)\n", basename(src),
            nrow(mcc), ncol(mcc)))

fp_csv <- file.path(out_dir, "multi_cohort_consistency.csv")
fp_xz  <- paste0(fp_csv, ".xz")
data.table::fwrite(mcc, fp_csv)
if (file.exists(fp_xz)) unlink(fp_xz)
status <- system2("xz", c("-9", "-e", "-f", shQuote(fp_csv)))
if (status != 0L || !file.exists(fp_xz)) {
  stop("xz compression failed for ", fp_csv, call. = FALSE)
}

cat(sprintf("bundled: %s  (%s)\n", basename(fp_xz),
            format(structure(file.info(fp_xz)$size,
                              class = "object_size"),
                   units = "auto")))
