#' Fit a stage-aware DESeq2 LRT model on user RNA-seq data
#'
#' Wraps the `phase33` canonical pipeline: `dataset + stage_group`
#' fixed-effect model on raw counts, LRT against `~ dataset` reduced
#' model, return per-gene coefficients (beta_E, beta_M, beta_L vs Normal) +
#' BH-adjusted LRT padj.
#'
#' @param counts Integer matrix or data.frame: genes (rows) x samples
#'   (cols). Row names are gene symbols.
#' @param stage Character/factor vector of stage labels per sample.
#'   Required levels: at minimum `"Normal"`, `"Early"`, `"Mid"`,
#'   `"Late"`. Other levels are dropped.
#' @param cohort Character/factor vector of dataset/cohort labels
#'   per sample. Used as a fixed effect to absorb cohort-level
#'   technical variance. If only one cohort, omit (default `NULL`).
#' @param min_count Numeric. Pre-filter: genes with row sum below
#'   `min_count` are dropped before LRT. Default 10.
#' @param padj_cutoff Numeric. BH-adjusted LRT padj cutoff for
#'   "stage-progressive" gene set. Default 0.05.
#' @return A `data.table` with one row per gene and columns:
#'   `gene_symbol, beta_N, beta_E, beta_M, beta_L, lrt_padj`,
#'   `lrt_significant` (logical, padj < cutoff), and lfcSE_*.
#' @examples
#' \dontrun{
#'   fit <- fit_stage_de(my_counts, my_stage, my_cohort)
#'   pat <- classify_trajectory(fit)
#' }
#' @export
fit_stage_de <- function(counts, stage, cohort = NULL,
                          min_count = 10,
                          padj_cutoff = 0.05) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required for fit_stage_de(). ",
         "Install via BiocManager::install('DESeq2').")
  }
  stage <- factor(stage, levels = c("Normal", "Early", "Mid", "Late"))
  if (any(is.na(stage))) {
    keep <- !is.na(stage)
    counts <- counts[, keep, drop = FALSE]
    stage  <- stage[keep]
    if (!is.null(cohort)) cohort <- cohort[keep]
    message(sprintf(
      "Dropped %d samples with stage outside Normal/Early/Mid/Late.",
      sum(!keep)))
  }
  if (length(stage) < 8L) stop(
    "Need at least 8 samples across stages for DESeq2 LRT.")

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
}
