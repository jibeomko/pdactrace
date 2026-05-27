#' Convert pdactrace_reference into a SummarizedExperiment
#'
#' Returns a `SummarizedExperiment` view of the bundled v0.4.0 atlas
#' in which the four stage-effect coefficient columns
#' (`rna_beta_N`, `rna_beta_E`, `rna_beta_M`, `rna_beta_L`) become a
#' 4-column `assay`, the standard errors `rna_lfcSE_E/M/L` (with `NA`
#' for the `Normal` reference column) become a parallel `assay`, and
#' all of the per-gene metadata (audit components, RNA / protein
#' trajectory pattern, translation class, scRNA cell origin, serum
#' direction, meta-analysis I2, etc.) goes into `rowData()`. The
#' `colData()` carries one row per stage with a `reference_level`
#' flag for the Normal column.
#'
#' Intended use: Bioconductor-native consumption of the atlas with
#' `assay()`, `rowData()`, `subsetByOverlaps()` (when paired with
#' a separate ranges object), and downstream tooling that expects
#' the SE container — without forcing a one-time mass conversion of
#' the bundled object (which remains a `data.table` for fast
#' query-based use cases).
#'
#' @param reference Optional override of the bundled
#'   `pdactrace_reference`. Default `NULL` uses the bundled atlas.
#'   Useful for unit tests on a synthetic atlas-shaped table.
#' @return A `SummarizedExperiment` with:
#'   * 2 assays — `rna_beta` (10,113 genes by 4 stages) and
#'     `rna_lfcSE` (same shape; `NA` in the `Normal` column).
#'   * `colData` — 4 rows: `stage = c("Normal", "Early", "Mid", "Late")`,
#'     `reference_level = c(TRUE, FALSE, FALSE, FALSE)`.
#'   * `rowData` — every non-stage-axis column of the bundled atlas
#'     (~109 columns: identifiers, audit components, trajectory
#'     pattern calls, meta-analysis summaries, etc.).
#'   * `metadata` — `list(atlas_version, build_date, n_cohorts, ...)`.
#' @examples
#' se <- as_summarized_experiment()
#' dim(se)
#' SummarizedExperiment::assayNames(se)
#'
#' \donttest{
#'   if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'     se <- as_summarized_experiment()
#'     se
#'     SummarizedExperiment::assay(se, "rna_beta")[1:3, ]
#'     head(SummarizedExperiment::rowData(se))
#'   }
#' }
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
as_summarized_experiment <- function(reference = NULL) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment is required. Install via ",
         "BiocManager::install('SummarizedExperiment').", call. = FALSE)
  }
  ref <- .get_reference(reference)

  # Stage-axis assay (genes x 4 stages) -----------------------------
  beta_cols <- c("rna_beta_N", "rna_beta_E", "rna_beta_M", "rna_beta_L")
  for (cn in beta_cols) {
    if (!cn %in% names(ref)) {
      stop(sprintf(
        "Reference is missing required stage-effect column '%s'.",
        cn), call. = FALSE)
    }
  }
  rna_beta <- as.matrix(ref[, ..beta_cols])
  rownames(rna_beta) <- ref$gene_symbol
  colnames(rna_beta) <- c("Normal", "Early", "Mid", "Late")

  # SE assay; pad NA into the Normal column (no SE for reference) ---
  se_cols <- c("rna_lfcSE_E", "rna_lfcSE_M", "rna_lfcSE_L")
  if (all(se_cols %in% names(ref))) {
    rna_lfcSE <- cbind(
      Normal = NA_real_,
      Early  = ref$rna_lfcSE_E,
      Mid    = ref$rna_lfcSE_M,
      Late   = ref$rna_lfcSE_L)
    rownames(rna_lfcSE) <- ref$gene_symbol
  } else {
    rna_lfcSE <- NULL
  }

  # colData ---------------------------------------------------------
  cd <- S4Vectors::DataFrame(
    stage = c("Normal", "Early", "Mid", "Late"),
    reference_level = c(TRUE, FALSE, FALSE, FALSE))
  rownames(cd) <- cd$stage

  # rowData = everything except stage-axis numeric assays -----------
  drop_cols <- c(beta_cols, se_cols)
  keep_cols <- setdiff(names(ref), drop_cols)
  rd <- S4Vectors::DataFrame(ref[, ..keep_cols])
  rownames(rd) <- ref$gene_symbol

  # Build SE --------------------------------------------------------
  assays_list <- if (!is.null(rna_lfcSE)) {
    list(rna_beta = rna_beta, rna_lfcSE = rna_lfcSE)
  } else {
    list(rna_beta = rna_beta)
  }

  meta <- tryCatch({
    am <- list_atlas_metadata()
    list(atlas_version = am$version,
          build_date    = am$build_date,
          n_cohorts     = am$n_cohorts,
          source_repo   = "https://github.com/jibeomko/pdactrace")
  }, error = function(e) list(
    atlas_version = NA_character_,
    source_repo   = "https://github.com/jibeomko/pdactrace"))

  SummarizedExperiment::SummarizedExperiment(
    assays  = assays_list,
    colData = cd,
    rowData = rd,
    metadata = meta)
}
