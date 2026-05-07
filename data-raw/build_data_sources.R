# build_data_sources.R ----------------------------------------------
# Builds the bundled `pdactrace_data_sources` table that backs
# `list_data_sources()`. One row per public dataset that contributed
# to the v0.99.x reference atlas. Columns are designed so a Bioc
# reviewer can answer "where does this data come from?"
# programmatically rather than by chasing prose in the manuscript.
#
# Source-of-truth: `analysis/manuscript/BiB_framework/00_overview_summary.md`
# and the per-figure manifest tsv files. If a new cohort is added,
# update this script and re-run.

pdactrace_data_sources <- data.table::data.table(
  accession = c(
    # ---- RNA-seq cohorts (Tier 1 + Tier 2 + Tier 3) ----
    "TCGA-PAAD",     "GSE62452",      "GSE71729",
    "GSE15471",      "GSE16515",      "CPTAC-PDA-RNA",
    # ---- Tissue proteomics ----
    "PXD006511",     "PXD026311",     "CPTAC-PDA-Prot",
    # ---- Single-cell RNA ----
    "GSE111672",     "GSE154778",     "GSE155698",
    "GSE194247",     "GSE205013",     "GSE212966",
    "PRJNA784866",   "GSE176127",     "PMID-32675368",
    "GSE196007",
    # ---- Serum proteomics ----
    "PXD025705",     "PXD003626",     "MSV000084438",
    # ---- Pancreatitis context ----
    "GSE83796",      "PXD012789",
    # ---- Validation ----
    "ICGC-PACA-CA",  "GTEx-pancreas"
  ),
  layer = c(
    rep("RNA",        6L),
    rep("Protein",    3L),
    rep("scRNA",     10L),
    rep("Serum",      3L),
    rep("Pancreatitis", 2L),
    rep("Validation", 2L)),
  source_type = c(
    rep("GEO",        2L), "GEO", "GEO", "GEO", "PDC",
    "PRIDE", "PRIDE", "PDC",
    rep("GEO",        6L), "SRA", "GEO", "Other", "GEO",
    "PRIDE", "PRIDE", "MassIVE",
    "GEO", "PRIDE",
    "ICGC", "GTEx"),
  used_for = c(
    "Discovery + 12-template fit",
    "Tier 2 validation",
    "Tier 2 validation",
    "Tier 1 validation",
    "Tier 1 validation",
    "RNA + Protein paired Tier 1 validation",
    "Tier 1 tissue proteomics",
    "Tier 2 tissue proteomics",
    "Tier 1 paired tissue proteomics",
    rep("Cell-origin atlas (scVI integration)", 10L),
    "Serum DIA Tier 1",
    "Serum LC-MS/MS Tier 2",
    "Serum panel orthogonal validation",
    "Pancreatitis tissue RNA",
    "Pancreatitis serum proteomics",
    "Survival validation",
    "Normal-pancreas baseline"),
  citation = NA_character_,
  url = c(
    "https://portal.gdc.cancer.gov/projects/TCGA-PAAD",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62452",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71729",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15471",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16515",
    "https://pdc.cancer.gov/pdc/study/PDC000270",
    "https://www.ebi.ac.uk/pride/archive/projects/PXD006511",
    "https://www.ebi.ac.uk/pride/archive/projects/PXD026311",
    "https://pdc.cancer.gov/pdc/study/PDC000270",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111672",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154778",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155698",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194247",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205013",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212966",
    "https://www.ebi.ac.uk/ena/browser/view/PRJNA784866",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176127",
    "https://pubmed.ncbi.nlm.nih.gov/32675368",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196007",
    "https://www.ebi.ac.uk/pride/archive/projects/PXD025705",
    "https://www.ebi.ac.uk/pride/archive/projects/PXD003626",
    "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000084438",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83796",
    "https://www.ebi.ac.uk/pride/archive/projects/PXD012789",
    "https://dcc.icgc.org/projects/PACA-CA",
    "https://www.gtexportal.org/home/tissue/Pancreas"))

usethis::use_data(pdactrace_data_sources, overwrite = TRUE,
                   compress = "xz")
cat(sprintf("pdactrace_data_sources: %d rows x %d cols\n",
            nrow(pdactrace_data_sources),
            ncol(pdactrace_data_sources)))
