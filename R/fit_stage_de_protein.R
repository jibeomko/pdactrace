#' Fit a stage-aware limma model on user tissue-protein intensity data
#'
#' Wraps the `phase34` canonical pipeline: `cohort + stage_group`
#' fixed-effect model on log-transformed protein intensities, returning
#' per-gene contrasts (beta_E, beta_M, beta_L vs Normal) plus
#' BH-adjusted F-test padj on the stage_group main effect.
#'
#' Designed to be the protein-side parallel of [fit_stage_de()], so
#' downstream `classify_trajectory()` / `assemble_user_evidence()` can
#' consume both layers symmetrically. Two interfaces are supported via
#' S4 dispatch: a numeric matrix / data.frame (the original v0.99.1
#' signature) or a `SummarizedExperiment` (new in v0.99.2; the
#' Bioconductor-native entry point).
#'
#' @param object Either a numeric intensity matrix/data.frame
#'   (genes by samples; rownames = gene symbols; values are log2
#'   intensities) or a `SummarizedExperiment` whose `assay()` is the
#'   intensity matrix and whose `colData()` carries stage / cohort.
#' @param ... Method-specific arguments — see the matrix and
#'   SummarizedExperiment interfaces below.
#' @param stage *(matrix/data.frame interface)* Character or factor
#'   vector of stage labels per sample. Required levels: at minimum
#'   `"Normal"`, `"Early"`, `"Mid"`, `"Late"`. Other levels are dropped.
#' @param cohort *(matrix/data.frame interface)* Optional character /
#'   factor vector of dataset/cohort labels per sample. `NULL` for
#'   single-cohort.
#' @param stage_col *(SummarizedExperiment interface)* Name of the
#'   column in `colData(object)` that carries the stage labels.
#' @param cohort_col *(SummarizedExperiment interface)* Optional name
#'   of the column in `colData(object)` that carries cohort labels.
#'   Default `NULL`.
#' @param assay_name *(SummarizedExperiment interface)* Name of the
#'   assay slot to pull as the intensity matrix. Default `"intensity"`.
#' @param min_nonNA Numeric. Pre-filter: rows with fewer than this many
#'   non-NA samples are dropped. Default 10.
#' @param padj_cutoff Numeric. BH-adjusted F-test padj cutoff for
#'   "stage-progressive" gene set. Default 0.05.
#' @return A `data.table` with the same column schema as
#'   [fit_stage_de()]: `gene_symbol, beta_N, beta_E, beta_M, beta_L,
#'   lfcSE_E, lfcSE_M, lfcSE_L, lrt_padj, lrt_significant`.
#' @examples
#' data(toy_protein)
#' data(toy_coldata)
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   prot_fit <- fit_stage_de_protein(toy_protein,
#'                                    toy_coldata$stage,
#'                                    toy_coldata$cohort)
#'   head(prot_fit)
#' }
#'
#' \donttest{
#'   # Matrix interface:
#'   prot_fit <- fit_stage_de_protein(my_intensity, my_stage, my_cohort)
#'
#'   # SummarizedExperiment interface:
#'   library(SummarizedExperiment)
#'   se <- SummarizedExperiment(
#'     assays = list(intensity = my_intensity),
#'     colData = DataFrame(stage = my_stage, cohort = my_cohort))
#'   prot_fit <- fit_stage_de_protein(
#'     se, stage_col = "stage", cohort_col = "cohort",
#'     assay_name = "intensity")
#' }
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
setGeneric("fit_stage_de_protein",
           function(object, ...) standardGeneric("fit_stage_de_protein"))

#' @rdname fit_stage_de_protein
#' @export
setMethod("fit_stage_de_protein", "ANY",
          function(object, stage, cohort = NULL,
                   min_nonNA = 10, padj_cutoff = 0.05) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("limma is required for fit_stage_de_protein(). ",
         "Install via BiocManager::install('limma').")
  }
  intensity <- object
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
})

#' @rdname fit_stage_de_protein
#' @export
setMethod("fit_stage_de_protein", "SummarizedExperiment",
          function(object, stage_col, cohort_col = NULL,
                   assay_name = "intensity",
                   min_nonNA = 10, padj_cutoff = 0.05) {
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
  intensity <- SummarizedExperiment::assay(object, assay_name)
  stage  <- as.character(cd[[stage_col]])
  cohort <- if (!is.null(cohort_col)) as.character(cd[[cohort_col]])
            else NULL
  fit_stage_de_protein(intensity, stage = stage, cohort = cohort,
                        min_nonNA = min_nonNA,
                        padj_cutoff = padj_cutoff)
})
