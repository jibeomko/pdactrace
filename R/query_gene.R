#' Query the pdactrace reference atlas for a single gene
#'
#' Returns a layered list with all evidence available for `gene_symbol`
#' across the requested layers. If the gene is not in the LRT-significant
#' atlas universe, returns `NULL` with an informative message.
#'
#' @param gene_symbol HGNC standard symbol (case-sensitive).
#' @param layers Character vector subset of
#'   `c("rna", "protein", "scrna", "serum", "clinical", "filter_status",
#'   "annotation")`. Default: all seven.
#' @return Named list with one slot per requested layer plus `$summary`
#'   (1-line text) and `$provenance` (comma-separated source phase
#'   scripts). Returns `NULL` invisibly if `gene_symbol` is not in the
#'   atlas.
#' @examples
#'   q <- query_gene("LTBP1")
#'   q$rna$pattern        # "Early_Burst_Up"
#'   q$serum$translation_class  # "B"
#' @export
query_gene <- function(gene_symbol,
                        layers = c("rna", "protein", "scrna", "serum",
                                    "clinical", "filter_status",
                                    "annotation")) {
  layers <- match.arg(layers, several.ok = TRUE)
  ref <- .get_reference()
  idx <- which(ref$gene_symbol == gene_symbol)
  row <- ref[idx]

  if (nrow(row) == 0L) {
    message(sprintf(
      "No evidence for '%s' in pdactrace atlas (v%s). The atlas covers ",
      gene_symbol, .pkg_version()),
      "LRT-significant tissue-expressed genes only.")
    return(invisible(NULL))
  }

  out <- list()

  if ("rna" %in% layers) {
    out$rna <- data.table::data.table(
      pattern         = row$rna_pattern,
      pattern_rho     = row$rna_pattern_rho,
      runner_up_rho   = row$rna_pattern_rho_runner_up,
      lrt_padj        = row$rna_lrt_padj,
      beta_N          = row$rna_beta_N,
      beta_E          = row$rna_beta_E,
      beta_M          = row$rna_beta_M,
      beta_L          = row$rna_beta_L,
      stouffer_z      = row$rna_stouffer_z,
      cohort_agreement = row$rna_cohort_agreement,
      excluded_mid_pattern = row$excluded_mid_pattern)
  }

  if ("protein" %in% layers) {
    out$protein <- data.table::data.table(
      pattern    = row$prot_pattern,
      tier       = row$prot_tier,
      concordant = row$rnaprot_concordant)
  }

  if ("scrna" %in% layers) {
    out$scrna <- data.table::data.table(
      top_celltype = row$cell_origin_top,
      enrichment_padj = row$cell_origin_padj)
    out$scrna$distribution <- list(row$cell_origin_distrib[[1]])
  }

  if ("serum" %in% layers) {
    out$serum <- data.table::data.table(
      detected           = row$serum_detected,
      n_cohorts_detected = row$serum_n_cohorts_detected,
      log2fc_PDAC_vs_HC  = row$serum_log2fc_PDAC_vs_HC,
      log2fc_Pan_vs_HC   = row$serum_log2fc_Pan_vs_HC,
      translation_class  = row$translation_class,
      phase77_strict     = row$phase77_strict)
  }

  if ("clinical" %in% layers) {
    out$clinical <- data.table::data.table(
      resectable_marker          = row$resectable_marker,
      resectable_pattern_phase29 = row$resectable_pattern_phase29,
      panel_member               = row$panel_member,
      evidence_scope             = row$evidence_scope)
  }

  if ("filter_status" %in% layers) {
    out$filter_status <- data.table::data.table(
      signal_peptide      = row$flt_signal_peptide,
      serum_measurable    = row$flt_serum_measurable,
      serum_significant   = row$flt_serum_significant,
      pancreatitis_pdac   = row$flt_pancreatitis_pdac,
      pancreatitis_hc     = row$flt_pancreatitis_hc,
      direction_match     = row$flt_direction_match,
      final_pass          = row$flt_final)
    if (all(is.na(unlist(out$filter_status)))) {
      attr(out$filter_status, "note") <-
        "Gene not evaluated by phase60 7-step pipeline (outside funnel)."
    }
  }

  if ("annotation" %in% layers) {
    out$annotation <- data.table::data.table(
      signalp_score = row$ann_signalp_score,
      pool_logfc    = row$ann_pool_logfc,
      pool_padj     = row$ann_pool_padj)
  }

  out$summary    <- .one_line_summary(row)
  out$provenance <- row$provenance

  class(out) <- c("pdactrace_gene_evidence", "list")
  out
}

#' @export
print.pdactrace_gene_evidence <- function(x, ...) {
  cat(x$summary, "\n")
  cat("Evidence:  ", format_provenance(x$provenance, "compact"),
      "\n", sep = "")
  cat("Technical: ", format_provenance(x$provenance, "raw"),
      "\n", sep = "")
  cat("\nLayers loaded:", paste(setdiff(names(x),
                                          c("summary", "provenance")),
                                  collapse = ", "), "\n")
  cat("Use $rna / $protein / $scrna / $serum / $clinical /\n")
  cat("$filter_status / $annotation for full evidence.\n")
  invisible(x)
}

# -- Internal helpers -----------------------------------------
.get_reference <- function(reference = NULL) {
  if (!is.null(reference)) {
    if (!data.table::is.data.table(reference)) {
      reference <- data.table::as.data.table(reference)
    }
    return(reference)
  }
  e <- new.env()
  data("pdactrace_reference", package = "pdactrace", envir = e)
  ref <- e$pdactrace_reference
  if (!data.table::is.data.table(ref)) {
    ref <- data.table::as.data.table(ref)
  }
  ref
}

.pkg_version <- function() {
  tryCatch(as.character(utils::packageVersion("pdactrace")),
            error = function(e) "0.1.0")
}

.one_line_summary <- function(row) {
  pat <- if (is.na(row$rna_pattern)) "no Early-onset pattern" else
    sprintf("%s (rho=%.2f, padj=%.1e)",
            row$rna_pattern, row$rna_pattern_rho, row$rna_lrt_padj)
  tier <- if (is.na(row$prot_tier)) "no protein tier" else row$prot_tier
  serum_part <- if (isTRUE(row$serum_detected)) {
    cls <- if (is.na(row$translation_class)) "" else
      sprintf(", Class %s", row$translation_class)
    sprintf("serum-detected (n=%d cohort%s)",
            row$serum_n_cohorts_detected, cls)
  } else {
    "no serum evidence"
  }
  sprintf("%s: %s | %s | %s.",
          row$gene_symbol, pat, tier, serum_part)
}
