#' Fit a stage-aware limma model on user tissue-protein intensity data
#'
#' Wraps the `phase34` canonical pipeline: `cohort + stage_group`
#' fixed-effect model on log-transformed protein intensities, returning
#' per-gene contrasts (beta_E, beta_M, beta_L vs Normal) plus
#' BH-adjusted F-test padj on the stage_group main effect.
#'
#' Designed to be the protein-side parallel of [fit_stage_de()], so
#' downstream `classify_trajectory()` / `assemble_user_evidence()` can
#' consume both layers symmetrically.
#'
#' @param intensity Numeric matrix or data.frame: genes (rows) x
#'   samples (cols). Row names are gene symbols. Values are assumed
#'   log2-transformed and median-normalised; if your input is raw
#'   intensities, log-transform first.
#' @param stage Character/factor vector of stage labels per sample.
#'   Required levels: at minimum `"Normal"`, `"Early"`, `"Mid"`,
#'   `"Late"`. Other levels are dropped.
#' @param cohort Character/factor vector of dataset/cohort labels per
#'   sample. Used as a fixed effect to absorb cohort-level technical
#'   variance. If only one cohort, omit (default `NULL`).
#' @param min_nonNA Numeric. Pre-filter: rows with fewer than this
#'   many non-NA samples are dropped. Default 10.
#' @param padj_cutoff Numeric. BH-adjusted F-test padj cutoff for
#'   "stage-progressive" gene set. Default 0.05.
#' @return A `data.table` with the same column schema as
#'   [fit_stage_de()]: `gene_symbol, beta_N, beta_E, beta_M, beta_L,
#'   lfcSE_E, lfcSE_M, lfcSE_L, lrt_padj, lrt_significant`. The
#'   `lrt_padj` column reports the BH-adjusted F-test p-value on the
#'   stage_group main effect (`limma::topTable(..., coef =
#'   c("Early", "Mid", "Late"))`).
#' @examples
#' \dontrun{
#'   # User log2-intensity matrix and matched stage labels
#'   prot_fit <- fit_stage_de_protein(my_intensity, my_stage, my_cohort)
#' }
#' @export
fit_stage_de_protein <- function(intensity, stage, cohort = NULL,
                                  min_nonNA = 10,
                                  padj_cutoff = 0.05) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("limma is required for fit_stage_de_protein(). ",
         "Install via BiocManager::install('limma').")
  }
  stage <- factor(stage, levels = c("Normal", "Early", "Mid", "Late"))
  if (any(is.na(stage))) {
    keep <- !is.na(stage)
    intensity <- intensity[, keep, drop = FALSE]
    stage <- stage[keep]
    if (!is.null(cohort)) cohort <- cohort[keep]
    message(sprintf(
      "Dropped %d samples with stage outside Normal/Early/Mid/Late.",
      sum(!keep)))
  }
  if (length(stage) < 8L) stop(
    "Need at least 8 samples across stages for limma.")

  mat <- as.matrix(intensity)
  keep_g <- rowSums(!is.na(mat)) >= min_nonNA
  mat <- mat[keep_g, , drop = FALSE]

  if (!is.null(cohort)) {
    cohort <- factor(cohort)
    design <- stats::model.matrix(~ cohort + stage)
  } else {
    design <- stats::model.matrix(~ stage)
  }

  fit <- limma::lmFit(mat, design)
  fit <- limma::eBayes(fit)

  stage_coefs <- grep("^stage", colnames(design), value = TRUE)
  expected <- c("stageEarly", "stageMid", "stageLate")
  miss <- setdiff(expected, stage_coefs)
  if (length(miss) > 0L) stop(
    sprintf("Missing stage levels in design: %s. ",
            paste(miss, collapse = ", ")),
    "Each of Normal/Early/Mid/Late must have >=1 sample.")

  res_list <- lapply(expected, function(co) {
    tt <- limma::topTable(fit, coef = co, number = Inf, sort.by = "none")
    tt[rownames(mat), , drop = FALSE]
  })
  names(res_list) <- expected

  f_tt <- limma::topTable(fit, coef = expected, number = Inf,
                           sort.by = "none")
  f_tt <- f_tt[rownames(mat), , drop = FALSE]

  out <- data.table::data.table(
    gene_symbol = rownames(mat),
    beta_N = 0,
    beta_E = res_list[["stageEarly"]]$logFC,
    beta_M = res_list[["stageMid"]]$logFC,
    beta_L = res_list[["stageLate"]]$logFC,
    lfcSE_E = sqrt(fit$s2.post) * fit$stdev.unscaled[, "stageEarly"],
    lfcSE_M = sqrt(fit$s2.post) * fit$stdev.unscaled[, "stageMid"],
    lfcSE_L = sqrt(fit$s2.post) * fit$stdev.unscaled[, "stageLate"],
    lrt_padj = stats::p.adjust(f_tt$P.Value, method = "BH"))
  out[, lrt_significant := !is.na(lrt_padj) & lrt_padj < padj_cutoff]
  attr(out, "padj_cutoff") <- padj_cutoff
  attr(out, "min_nonNA") <- min_nonNA
  attr(out, "n_samples") <- length(stage)
  attr(out, "n_cohorts") <- if (is.null(cohort)) 1L else
                             length(unique(cohort))
  out
}
