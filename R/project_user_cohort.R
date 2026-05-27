#' End-to-end user-cohort projection through the audit framework
#'
#' One-call wrapper that takes a user RNA-seq cohort (and optionally a
#' tissue-protein cohort) and walks it end-to-end through the four
#' v0.99 framework steps: stage-aware DE fit, 12-template trajectory
#' matching, evidence-table assembly, and frozen audit scoring.
#' Returns a structured list with each intermediate object plus the
#' final audit table, so the user can inspect at any step.
#'
#' This is a convenience around the chain
#' [fit_stage_de()] -> [classify_trajectory()] ->
#' [assemble_user_evidence()] -> [compute_audit_score()].
#' Calling those four explicitly continues to work; this wrapper
#' reduces a typical user-cohort workflow to a single line and
#' provides a consistent return shape for reporting and review.
#'
#' Bioconductor-native input is supported via the same S4 method
#' dispatch as `fit_stage_de()`: the `rna` argument can be either a
#' count matrix / data.frame plus `coldata`, or a
#' `SummarizedExperiment` whose `assay()` is the count matrix and
#' whose `colData()` carries the stage / cohort columns.
#'
#' @param rna Either an integer count matrix / data.frame (genes by
#'   samples) or a `SummarizedExperiment`. Required.
#' @param coldata Optional `data.frame` / `DataFrame` of sample-level
#'   metadata. Required when `rna` is a matrix; ignored when `rna` is
#'   a `SummarizedExperiment` (its `colData()` is used instead).
#' @param stage_col Name of the column in `coldata` (or
#'   `colData(rna)`) that carries the stage labels.
#' @param cohort_col Optional name of the cohort column. Default
#'   `NULL` (single-cohort).
#' @param protein Optional second `SummarizedExperiment` or numeric
#'   intensity matrix carrying log2 protein intensities for the same
#'   stage / cohort design. When supplied, the protein side runs
#'   through [fit_stage_de_protein()] +
#'   [classify_protein_trajectory()] and is layered into
#'   `assemble_user_evidence(rna_fit, prot_fit)`. Default `NULL`.
#' @param protein_assay_name Assay slot to use when `protein` is a
#'   `SummarizedExperiment`. Default `"intensity"`.
#' @param signal_peptide Optional character vector of gene symbols
#'   carrying a signal peptide; passed straight to
#'   [assemble_user_evidence()].
#' @param sig_only Passed to [classify_protein_trajectory()] (default
#'   `FALSE` here so small / noisy cohorts still surface
#'   best-template calls).
#' @return A `list` of class `pdactrace_user_projection` with
#'   elements:
#'   * `rna_fit`     -- output of `fit_stage_de()`.
#'   * `rna_pattern` -- output of `classify_trajectory()`.
#'   * `prot_fit`    -- output of `fit_stage_de_protein()` (or `NULL`).
#'   * `prot_pattern`-- output of `classify_protein_trajectory()`
#'                      (or `NULL`).
#'   * `evidence`    -- output of `assemble_user_evidence()`.
#'   * `audit`       -- output of `compute_audit_score(evidence = ...)`.
#'   * `summary`     -- one-row `data.table` with input shape +
#'                      audit-class counts.
#' @examples
#' data(toy_counts)
#' data(toy_coldata)
#' if (requireNamespace("DESeq2", quietly = TRUE)) {
#'   res <- project_user_cohort(
#'     rna = toy_counts, coldata = toy_coldata,
#'     stage_col = "stage", cohort_col = "cohort")
#'   res$summary
#' }
#'
#' \donttest{
#'   if (requireNamespace("DESeq2", quietly = TRUE)) {
#'     res <- project_user_cohort(
#'       rna        = toy_counts,
#'       coldata    = toy_coldata,
#'       stage_col  = "stage",
#'       cohort_col = "cohort",
#'       protein    = toy_protein)
#'     res$summary
#'     head(res$audit[order(-audit_score)])
#'   }
#' }
#' @export
project_user_cohort <- function(rna,
                                  coldata = NULL,
                                  stage_col,
                                  cohort_col = NULL,
                                  protein    = NULL,
                                  protein_assay_name = "intensity",
                                  signal_peptide = NULL,
                                  sig_only = FALSE) {
  # ── 1. RNA fit + trajectory ────────────────────────────────
  if (is(rna, "SummarizedExperiment")) {
    rna_fit <- fit_stage_de(rna,
                              stage_col = stage_col,
                              cohort_col = cohort_col)
  } else {
    if (is.null(coldata)) {
      stop("`coldata` is required when `rna` is a matrix or ",
           "data.frame. Either pass coldata or a SummarizedExperiment.",
           call. = FALSE)
    }
    if (!stage_col %in% colnames(coldata)) {
      stop(sprintf("coldata has no column '%s'.", stage_col),
           call. = FALSE)
    }
    cohort_vec <- if (!is.null(cohort_col)) coldata[[cohort_col]]
                  else NULL
    rna_fit <- fit_stage_de(rna,
                              stage  = coldata[[stage_col]],
                              cohort = cohort_vec)
  }
  rna_pattern <- classify_trajectory(rna_fit)

  # ── 2. Optional protein fit + trajectory ───────────────────
  prot_fit <- NULL
  prot_pattern <- NULL
  if (!is.null(protein)) {
    if (is(protein, "SummarizedExperiment")) {
      prot_fit <- fit_stage_de_protein(
        protein, stage_col = stage_col, cohort_col = cohort_col,
        assay_name = protein_assay_name)
    } else {
      prot_fit <- fit_stage_de_protein(
        protein,
        stage  = coldata[[stage_col]],
        cohort = if (!is.null(cohort_col)) coldata[[cohort_col]]
                else NULL)
    }
    prot_pattern <- classify_protein_trajectory(
      prot_fit, sig_only = sig_only)
  }

  # ── 3. Assemble evidence + audit score ─────────────────────
  evidence <- assemble_user_evidence(
    rna_fit  = rna_pattern,
    prot_fit = prot_pattern,
    signal_peptide = signal_peptide)
  audit <- compute_audit_score(evidence = evidence)

  # ── 4. Summary ─────────────────────────────────────────────
  cls_counts <- if ("audit_class" %in% names(audit)) {
    table(audit$audit_class, useNA = "ifany")
  } else integer(0)
  summary <- data.table::data.table(
    n_genes_input  = nrow(rna_fit),
    n_genes_audit  = nrow(audit),
    has_protein    = !is.null(prot_fit),
    n_high_conf    = unname(cls_counts["high_confidence"] %||% 0L),
    n_supp_uncert  = unname(cls_counts["supported_uncertain"] %||% 0L),
    n_penalized    = unname(cls_counts["penalized"] %||% 0L),
    n_excluded     = unname(cls_counts["excluded"] %||% 0L),
    n_low          = unname(cls_counts["low"] %||% 0L))

  out <- list(
    rna_fit      = rna_fit,
    rna_pattern  = rna_pattern,
    prot_fit     = prot_fit,
    prot_pattern = prot_pattern,
    evidence     = evidence,
    audit        = audit,
    summary      = summary)
  class(out) <- c("pdactrace_user_projection", "list")
  out
}

#' @export
print.pdactrace_user_projection <- function(x, ...) {
  cat("<pdactrace_user_projection>\n")
  cat(sprintf("  Input:   %d genes; protein layer: %s\n",
              x$summary$n_genes_input,
              if (x$summary$has_protein) "yes" else "no"))
  cat(sprintf("  Audit:   %d genes scored\n",
              x$summary$n_genes_audit))
  cat(sprintf("  Classes: high=%d  supp_unc=%d  penalized=%d  excluded=%d  low=%d\n",
              x$summary$n_high_conf,  x$summary$n_supp_uncert,
              x$summary$n_penalized,  x$summary$n_excluded,
              x$summary$n_low))
  invisible(x)
}

# Small NULL-coalesce helper, kept private.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
