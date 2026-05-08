# build_protein_betas.R --------------------------------------------
# Bundles the v0.99.x per-stage tissue-protein effect-size table that
# plot_template_atlas("protein") and plot_gene_template(layer="protein")
# need. Source-of-truth: phase34_protein_pooled_12template.csv from the
# manuscript-monorepo (PDAC_biomarker), which is the same pipeline that
# fed the bundled v0.99.x atlas's prot_pattern column.
#
# Re-run via:
#   PDAC_BASE_DIR=/home/kjb9412/PDAC_biomarker \
#     Rscript data-raw/build_protein_betas.R

.find_phase_csv <- function(stem) {
  # 1. Bundled in inst/extdata (canonical, self-contained path)
  bundled <- file.path("inst", "extdata",
                        paste0(stem, ".csv.xz"))
  if (file.exists(bundled)) return(bundled)
  # 2. Companion manuscript-monorepo (developer fallback)
  base_dir <- Sys.getenv("PDAC_BASE_DIR",
                          "/home/kjb9412/PDAC_biomarker")
  ext <- file.path(
    base_dir,
    "analysis/manuscript/tissue_to_serum_biomarker/results",
    paste0(stem, ".csv"))
  if (file.exists(ext)) return(ext)
  stop("Cannot find ", stem,
       ".csv(.xz). Re-run data-raw/bundle_phase_csvs.R or set ",
       "PDAC_BASE_DIR.", call. = FALSE)
}
src_path <- .find_phase_csv("phase34_protein_pooled_12template")

src <- if (grepl("\\.xz$", src_path)) {
  data.table::fread(cmd = paste("xz -dc", shQuote(src_path)))
} else {
  data.table::fread(src_path)
}
stopifnot(all(c("gene", "beta_E", "beta_M", "beta_L",
                 "prot_pat", "prot_rho") %in% names(src)))

pdactrace_protein_betas <- data.table::data.table(
  gene_symbol      = src$gene,
  prot_beta_N      = 0,
  prot_beta_E      = src$beta_E,
  prot_beta_M      = src$beta_M,
  prot_beta_L      = src$beta_L,
  prot_pattern_12  = src$prot_pat,
  prot_pattern_rho = src$prot_rho,
  prot_lrt_padj    = src$F_padj)
data.table::setkey(pdactrace_protein_betas, gene_symbol)

usethis::use_data(pdactrace_protein_betas,
                   overwrite = TRUE, compress = "xz")

cat(sprintf("pdactrace_protein_betas: %d rows x %d cols\n",
            nrow(pdactrace_protein_betas),
            ncol(pdactrace_protein_betas)))
cat("Pattern distribution:\n")
print(sort(table(pdactrace_protein_betas$prot_pattern_12,
                  useNA = "ifany"), decreasing = TRUE))
