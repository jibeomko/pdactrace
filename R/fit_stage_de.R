#' Fit a stage-aware DESeq2 LRT model on user RNA-seq data
#'
#' Wraps the `phase33` canonical pipeline: `dataset + stage_group`
#' fixed-effect model on raw counts, LRT against `~ dataset` reduced
#' model, return per-gene coefficients (beta_E, beta_M, beta_L vs Normal) +
#' BH-adjusted LRT padj.
#'
#' Two interfaces are supported via S4 method dispatch on `object`:
#'
#' * **Matrix interface** (the default): `fit_stage_de(counts, stage,
#'   cohort)`. `counts` is a numeric matrix or data.frame with genes on
#'   rows and samples on columns; `stage` and `cohort` are vectors of
#'   length `ncol(counts)`. This is the original v0.99.1 signature and
#'   is preserved unchanged for backwards compatibility.
#' * **SummarizedExperiment interface** (new in v0.99.2):
#'   `fit_stage_de(se, stage_col = "stage", cohort_col = "cohort")`.
#'   The count matrix is pulled from `assay(se, assay_name)` and the
#'   stage / cohort vectors are pulled from `colData(se)`. This is the
#'   Bioconductor-native entry point.
#'
#' @param object Either a count matrix/data.frame (genes by samples,
#'   rownames = gene symbols) or a `SummarizedExperiment` whose
#'   `assay()` is the count matrix and whose `colData()` carries the
#'   stage / cohort columns.
#' @param ... Method-specific arguments — see Details for the matrix
#'   and SummarizedExperiment interfaces.
#' @param stage *(matrix/data.frame interface)* Character or factor
#'   vector of stage labels per sample. Required levels: at minimum
#'   `"Normal"`, `"Early"`, `"Mid"`, `"Late"`. Other levels are dropped.
#' @param cohort *(matrix/data.frame interface)* Optional character /
#'   factor vector of dataset/cohort labels per sample. Used as a fixed
#'   effect to absorb cohort-level technical variance. If only one
#'   cohort, omit (default `NULL`).
#' @param stage_col *(SummarizedExperiment interface)* Name of the
#'   column in `colData(object)` that carries the stage labels.
#' @param cohort_col *(SummarizedExperiment interface)* Optional name
#'   of the column in `colData(object)` that carries the cohort labels.
#'   Default `NULL` (single-cohort).
#' @param assay_name *(SummarizedExperiment interface)* Name of the
#'   assay slot to pull as the count matrix. Default `"counts"`.
#' @param min_count Numeric. Pre-filter: genes with row sum below
#'   `min_count` are dropped before LRT. Default 10.
#' @param padj_cutoff Numeric. BH-adjusted LRT padj cutoff for
#'   "stage-progressive" gene set. Default 0.05.
#' @return A `data.table` with one row per gene and columns:
#'   `gene_symbol, beta_N, beta_E, beta_M, beta_L, lrt_padj`,
#'   `lrt_significant` (logical, padj < cutoff), and lfcSE_*.
#' @examples
#' data(toy_counts)
#' data(toy_coldata)
#' if (requireNamespace("DESeq2", quietly = TRUE)) {
#'   fit <- fit_stage_de(toy_counts, toy_coldata$stage,
#'                       toy_coldata$cohort)
#'   head(fit)
#' }
#'
#' \donttest{
#'   # Matrix interface (original signature):
#'   fit <- fit_stage_de(my_counts, my_stage, my_cohort)
#'
#'   # SummarizedExperiment interface (Bioconductor-native):
#'   library(SummarizedExperiment)
#'   se <- SummarizedExperiment(
#'     assays = list(counts = my_counts),
#'     colData = DataFrame(stage = my_stage, cohort = my_cohort))
#'   fit <- fit_stage_de(se, stage_col = "stage", cohort_col = "cohort")
#'   pat <- classify_trajectory(fit)
#' }
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
setGeneric("fit_stage_de",
           function(object, ...) standardGeneric("fit_stage_de"))

#' @rdname fit_stage_de
#' @export
setMethod("fit_stage_de", "ANY",
          function(object, stage, cohort = NULL,
                   min_count = 10, padj_cutoff = 0.05) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required for fit_stage_de(). ",
         "Install via BiocManager::install('DESeq2').")
  }
  counts <- object
  # ---- friendly stage-label check (v0.99.9) -----------------------
  # The 12-template trajectory framework is stage-aware; binary
  # Normal vs Tumor input does NOT apply. Surface the mismatch
  # explicitly instead of silent factor() coercion + downstream
  # "Need at least 8 samples" which obscures the root cause.
  stage_in <- as.character(stage)
  required <- c("Normal", "Early", "Mid", "Late")
  unmatched <- setdiff(unique(stage_in), required)
  matched   <- intersect(unique(stage_in), required)
  if (length(matched) < 2L) {
    msg <- sprintf(
      paste0("Stage labels you provided: %s\n",
             "Required levels:           %s\n",
             "Matched: %d / %d unique stage values.\n\n",
             "pdactrace tissue layer needs stage-aware data ",
             "(Normal/Early/Mid/Late, >=2 levels with samples).\n",
             "  - If you only have Normal vs Tumor: this framework ",
             "does not apply at the tissue layer.\n",
             "  - For serum proteomics (binary group contrasts), ",
             "use project_user_serum_cohort() instead.\n",
             "  - To map clinical TNM/AJCC stages to ",
             "Normal/Early/Mid/Late, see vignette('user_cohort_extension')."),
      paste(sort(unique(stage_in)), collapse = ", "),
      paste(required, collapse = ", "),
      length(matched), length(unique(stage_in)))
    stop(msg, call. = FALSE)
  }

  stage <- factor(stage_in, levels = required)
  if (any(is.na(stage))) {
    keep <- !is.na(stage)
    counts <- counts[, keep, drop = FALSE]
    stage  <- stage[keep]
    if (!is.null(cohort)) cohort <- cohort[keep]
    message(sprintf(
      "Dropped %d samples with stage labels not in Normal/Early/Mid/Late: %s",
      sum(!keep), paste(unmatched, collapse = ", ")))
  }
  if (length(stage) < 8L) stop(
    "Need at least 8 samples across stages for DESeq2 LRT. ",
    "After mapping to Normal/Early/Mid/Late you have ",
    length(stage), " samples.", call. = FALSE)

  cd <- data.frame(stage_group = stage)
  if (!is.null(cohort)) cd$cohort <- factor(cohort)
  rownames(cd) <- colnames(counts)

  full_formula <- if (!is.null(cohort)) ~ cohort + stage_group
                   else ~ stage_group
  reduced_formula <- if (!is.null(cohort)) ~ cohort else ~ 1

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData   = cd,
    design    = full_formula)
  keep_g <- rowSums(DESeq2::counts(dds)) >= min_count
  dds <- dds[keep_g, ]
  dds <- DESeq2::DESeq(dds, test = "LRT", reduced = reduced_formula,
                        quiet = TRUE)

  resE <- DESeq2::results(dds, name = "stage_group_Early_vs_Normal")
  resM <- DESeq2::results(dds, name = "stage_group_Mid_vs_Normal")
  resL <- DESeq2::results(dds, name = "stage_group_Late_vs_Normal")

  out <- data.table::data.table(
    gene_symbol = rownames(dds),
    beta_N      = 0,
    beta_E      = resE$log2FoldChange,
    beta_M      = resM$log2FoldChange,
    beta_L      = resL$log2FoldChange,
    lfcSE_E     = resE$lfcSE,
    lfcSE_M     = resM$lfcSE,
    lfcSE_L     = resL$lfcSE,
    lrt_padj    = stats::p.adjust(
      DESeq2::results(dds)$pvalue, method = "BH"))
  out[, lrt_significant := !is.na(lrt_padj) & lrt_padj < padj_cutoff]
  attr(out, "padj_cutoff") <- padj_cutoff
  attr(out, "min_count")   <- min_count
  attr(out, "n_samples")   <- length(stage)
  attr(out, "n_cohorts")   <- if (is.null(cohort)) 1L else
                                length(unique(cohort))
  out
})

#' @rdname fit_stage_de
#' @export
setMethod("fit_stage_de", "SummarizedExperiment",
          function(object, stage_col, cohort_col = NULL,
                   assay_name = "counts",
                   min_count = 10, padj_cutoff = 0.05) {
  cd <- SummarizedExperiment::colData(object)
  if (!stage_col %in% colnames(cd)) {
    stop(sprintf(
      "colData(object) does not have a column named '%s'. ",
      stage_col),
      "Available columns: ",
      paste(colnames(cd), collapse = ", "), call. = FALSE)
  }
  if (!is.null(cohort_col) && !cohort_col %in% colnames(cd)) {
    stop(sprintf(
      "colData(object) does not have a column named '%s'. ",
      cohort_col),
      "Set cohort_col = NULL for single-cohort input.", call. = FALSE)
  }
  assay_names <- SummarizedExperiment::assayNames(object)
  if (!assay_name %in% assay_names) {
    stop(sprintf(
      "Assay '%s' not found in SummarizedExperiment. ", assay_name),
      "Available assays: ",
      paste(assay_names, collapse = ", "), call. = FALSE)
  }
  counts  <- SummarizedExperiment::assay(object, assay_name)
  stage   <- as.character(cd[[stage_col]])
  cohort  <- if (!is.null(cohort_col)) as.character(cd[[cohort_col]])
             else NULL
  fit_stage_de(counts, stage = stage, cohort = cohort,
                min_count = min_count, padj_cutoff = padj_cutoff)
})
