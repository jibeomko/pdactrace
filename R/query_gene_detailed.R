#' Detailed gene-level evidence - 5 granularity tables
#'
#' Returns a named list of `data.table`s exposing the *hidden* per-stage,
#' per-cohort, per-celltype, per-filter-step, and per-serum-cohort detail
#' that `query_gene()` summarizes only at the headline level.
#'
#' Use this for reviewer-grade audit (e.g., *"is the LTBP1 upregulation
#' consistent across all 4 RNA cohorts?"*) or for downstream computation.
#'
#' @param gene_symbol HGNC standard symbol (case-sensitive).
#' @return Named list with 5 slots (each a `data.table`) plus
#'   `$summary` (1-line text) and `$provenance`. Returns `NULL`
#'   invisibly with a message if `gene_symbol` is not in the atlas.
#' @section Slots:
#' * **per_stage** - 4 rows (Normal, Early, Mid, Late) with
#'   `log2FC`, `lfcSE`, `ci_lo`, `ci_hi`, `lrt_padj`, `wald_padj`,
#'   `significant`. `lrt_padj` is the omnibus DESeq2 LRT padj
#'   (identical across stages by construction); `wald_padj` is the
#'   per-contrast Wald padj from a `nbinomWaldTest()` refit of the
#'   same model (v0.4.0). The `significant` column uses Wald padj.
#' * **per_cohort** - 4 rows (TCGA, CPTAC, GSE224564, GSE79668) with
#'   `trend` (-1/0/+1), `monotonic` (lgl), plus Stouffer summary as
#'   table attribute.
#' * **per_celltype** - 11 rows with `mean_expression`, `pct_of_total`,
#'   ranked descending.
#' * **filter_diag** - 7 rows (phase60 steps) with `pass` (lgl),
#'   `underlying_metric`, `cutoff`, `note`.
#' * **serum_per_cohort** - 3 rows (PDAC, Pancreatitis, HC) with
#'   `mean`, `n_unknown` (cohort sample n not always available),
#'   plus PDAC-vs-Pan t-test p-value as attribute.
#' @examples
#'   d <- query_gene_detailed("LTBP1")
#'   d$per_stage      # log2FC + 95% CI per stage
#'   d$per_cohort     # cohort x trend x monotonic
#'   d$filter_diag    # 7-step audit with underlying numbers
#' @export
query_gene_detailed <- function(gene_symbol) {
  ref <- .get_reference()
  idx <- which(ref$gene_symbol == gene_symbol)
  row <- ref[idx]
  if (nrow(row) == 0L) {
    message(sprintf("No evidence for '%s' in pdactrace atlas (v%s).",
                     gene_symbol, .pkg_version()))
    return(invisible(NULL))
  }

  out <- list()

  # -- 1. per_stage: forest-ready table -------------------------
  # lrt_padj    = DESeq2 LRT padj (omnibus stage-progression test);
  #               identical across stages by construction.
  # wald_padj   = nbinomWaldTest padj for that contrast vs Normal (v0.4.0).
  beta  <- c(row$rna_beta_N, row$rna_beta_E, row$rna_beta_M, row$rna_beta_L)
  se    <- c(0, row$rna_lfcSE_E, row$rna_lfcSE_M, row$rna_lfcSE_L)
  lrt_padj <- c(NA_real_, row$rna_padj_E, row$rna_padj_M, row$rna_padj_L)
  wald_padj <- if (all(c("rna_wald_padj_E", "rna_wald_padj_M",
                           "rna_wald_padj_L") %in% names(row))) {
    c(NA_real_, row$rna_wald_padj_E, row$rna_wald_padj_M,
      row$rna_wald_padj_L)
  } else {
    rep(NA_real_, 4)
  }
  out$per_stage <- data.table::data.table(
    stage        = factor(c("Normal", "Early", "Mid", "Late"),
                            levels = c("Normal", "Early", "Mid", "Late")),
    log2FC       = beta,
    lfcSE        = se,
    ci_lo        = beta - 1.96 * se,
    ci_hi        = beta + 1.96 * se,
    lrt_padj     = lrt_padj,
    wald_padj    = wald_padj,
    significant  = !is.na(wald_padj) & wald_padj < 0.05)
  attr(out$per_stage, "lrt_padj_kind")  <- "LRT (omnibus stage progression)"
  attr(out$per_stage, "wald_padj_kind") <- "Wald (per-contrast vs Normal)"

  # -- 2. per_cohort: vote + monotonic --------------------------
  trend_lst <- row$rna_per_cohort_trend[[1]]
  mono_lst  <- row$rna_per_cohort_monotonic[[1]]
  if (!is.null(trend_lst) && length(trend_lst) > 0) {
    out$per_cohort <- data.table::data.table(
      cohort    = names(trend_lst),
      trend     = unlist(trend_lst, use.names = FALSE),
      monotonic = unlist(mono_lst,  use.names = FALSE))
    attr(out$per_cohort, "stouffer_z")    <- row$rna_stouffer_z
    attr(out$per_cohort, "stouffer_p")    <- row$rna_stouffer_p
    attr(out$per_cohort, "stouffer_padj") <- row$rna_stouffer_padj
    attr(out$per_cohort, "agreement_pct") <- row$rna_cohort_agreement
  } else {
    out$per_cohort <- data.table::data.table()
  }

  # -- 3. per_celltype: full 11-celltype expression -------------
  d <- row$cell_origin_distrib[[1]]
  if (!is.null(d) && length(d) > 0) {
    pct <- 100 * d / sum(d)
    out$per_celltype <- data.table::data.table(
      celltype        = names(d),
      mean_expression = as.numeric(d),
      pct_of_total    = as.numeric(pct))
    data.table::setorder(out$per_celltype, -mean_expression)
    attr(out$per_celltype, "specificity_tau") <- row$cell_specificity_tau
  } else {
    out$per_celltype <- data.table::data.table()
  }

  # -- 4. filter_diag: 7-step audit with underlying metric ------
  flt_steps <- c("signal_peptide", "serum_measurable",
                  "serum_significant", "pancreatitis_pdac",
                  "pancreatitis_hc", "direction_match", "final")
  flt_pass <- c(row$flt_signal_peptide, row$flt_serum_measurable,
                 row$flt_serum_significant, row$flt_pancreatitis_pdac,
                 row$flt_pancreatitis_hc, row$flt_direction_match,
                 row$flt_final)
  flt_metric <- c(
    sprintf("UniProt-SignalP: %s",
            ifelse(is.na(row$flt_signal_peptide), "NA",
                    ifelse(row$flt_signal_peptide, "yes", "no"))),
    sprintf("Detected in serum cohort union: %s",
            ifelse(is.na(row$flt_serum_measurable), "NA",
                    ifelse(row$flt_serum_measurable, "yes", "no"))),
    sprintf("pool padj = %s",
            if (is.na(row$ann_pool_padj)) "NA" else
              formatC(row$ann_pool_padj, format = "g", digits = 2)),
    sprintf("PDAC vs Pancreatitis t-pval = %s",
            if (is.na(row$ann_pdac_vs_pan_pval)) "NA" else
              formatC(row$ann_pdac_vs_pan_pval, format = "g", digits = 2)),
    sprintf("pan-vs-HC pval = %s",
            if (is.na(row$ann_pan_vs_hc_pval)) "NA" else
              formatC(row$ann_pan_vs_hc_pval, format = "g", digits = 2)),
    sprintf("Direction agreement: %s",
            ifelse(is.na(row$flt_direction_match), "NA",
                    ifelse(row$flt_direction_match, "yes", "no"))),
    "All upstream steps must pass")
  flt_cutoff <- c("UniProt SignalP positive",
                   ">= 1 cohort detection",
                   "padj < 0.05",
                   "t-pval < 0.05",
                   "pan-vs-HC ns or HC-in-middle",
                   "matching direction",
                   "all 6 upstream pass")
  out$filter_diag <- data.table::data.table(
    step              = flt_steps,
    pass              = flt_pass,
    underlying_metric = flt_metric,
    cutoff            = flt_cutoff)

  # -- 5. serum_per_cohort: PDAC / Pan / HC means ---------------
  out$serum_per_cohort <- data.table::data.table(
    arm  = c("PDAC", "Pancreatitis", "HC"),
    mean = c(row$ann_pdac_mean, row$ann_pan_mean, row$ann_hc_mean))
  attr(out$serum_per_cohort, "pdac_vs_pan_pval") <-
    row$ann_pdac_vs_pan_pval
  attr(out$serum_per_cohort, "log2fc_PDAC_vs_HC") <-
    row$serum_log2fc_PDAC_vs_HC
  attr(out$serum_per_cohort, "log2fc_Pan_vs_HC") <-
    row$serum_log2fc_Pan_vs_HC

  # -- 6. meta_analysis: random-effects meta + tier (v0.2.0) ----
  out$meta <- data.table::data.table(
    contrast       = c("Normal_vs_Early", "Mid_vs_Early", "Late_vs_Early"),
    beta_meta      = c(row$meta_NvE_beta, row$meta_MvE_beta,
                        row$meta_LvE_beta),
    pval           = c(row$meta_NvE_pval, row$meta_MvE_pval,
                        row$meta_LvE_pval),
    padj           = c(row$meta_NvE_padj, row$meta_MvE_padj,
                        row$meta_LvE_padj),
    I_squared      = c(row$meta_NvE_I2, row$meta_MvE_I2,
                        row$meta_LvE_I2),
    n_cohorts      = c(row$meta_NvE_k, row$meta_MvE_k, row$meta_LvE_k))
  attr(out$meta, "cohort_divergent")     <- row$meta_cohort_divergent
  attr(out$meta, "max_abs_beta_meta")    <- row$max_abs_beta_meta
  attr(out$meta, "max_I2_meta")          <- row$max_I2_meta
  attr(out$meta, "audit_class")          <- row$audit_class
  attr(out$meta, "audit_score")          <- row$audit_score

  out$summary    <- .one_line_summary(row)
  out$provenance <- row$provenance

  class(out) <- c("pdactrace_gene_detailed", "list")
  out
}

#' @export
print.pdactrace_gene_detailed <- function(x, ...) {
  cat(x$summary, "\n")
  cat("Provenance:", x$provenance, "\n\n")
  cat("Slots: $per_stage / $per_cohort / $per_celltype /",
      "$filter_diag / $serum_per_cohort\n")
  cat("(Each slot is a data.table; see ?query_gene_detailed.)\n")
  invisible(x)
}
