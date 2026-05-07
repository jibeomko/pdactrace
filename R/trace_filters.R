#' Trace which selection filters a gene passed or failed
#'
#' Returns the gene's pass/fail status for each step of pdactrace's
#' two parallel narrowing funnels:
#'
#' * **`phase60` 7-step pipeline** - strict cumulative funnel
#'   (signal peptide -> serum measurable -> pool significance -> vs
#'   pancreatitis -> HC-in-middle -> direction match -> final).
#' * **`phase77` strict candidate set** - independent narrow funnel
#'   producing the 22 strict tissue+serum convergent candidates
#'   that anchor manual Class B curation.
#'
#' Most genes (including LTBP1) pass *one* of the two funnels but not
#' both - this is the audit-trail message of pdactrace.
#'
#' @param gene Character vector of HGNC symbols.
#' @return A `data.table` with one row per requested gene and one
#'   column per filter step (logical). Includes `n_phase60_pass`
#'   (0-7) and `class_route` summary column.
#' @examples
#'   trace_filters("LTBP1")
#'   trace_filters(c("LTBP1", "SPARC", "CDH13", "SERPINA1"))
#' @export
trace_filters <- function(gene) {
  ref <- .get_reference()
  hit <- ref[ref$gene_symbol %in% gene, ]
  miss <- setdiff(gene, hit$gene_symbol)
  if (length(miss) > 0L) {
    message(sprintf("Genes not in atlas (skipped): %s",
                     paste(miss, collapse = ", ")))
  }
  if (nrow(hit) == 0L) return(data.table::data.table())

  out <- hit[, list(
    gene_symbol            = gene_symbol,
    # phase60 7-step funnel
    flt_signal_peptide     = flt_signal_peptide,
    flt_serum_measurable   = flt_serum_measurable,
    flt_serum_significant  = flt_serum_significant,
    flt_pancreatitis_pdac  = flt_pancreatitis_pdac,
    flt_pancreatitis_hc    = flt_pancreatitis_hc,
    flt_direction_match    = flt_direction_match,
    flt_final              = flt_final,
    # phase77 + manual paths
    phase77_strict         = phase77_strict,
    translation_class      = translation_class,
    panel_member           = panel_member,
    resectable_marker      = resectable_marker)]

  # Summary count + route classification
  flt_cols <- c("flt_signal_peptide", "flt_serum_measurable",
                 "flt_serum_significant", "flt_pancreatitis_pdac",
                 "flt_pancreatitis_hc", "flt_direction_match",
                 "flt_final")
  out[, n_phase60_pass := rowSums(.SD == TRUE, na.rm = TRUE),
      .SDcols = flt_cols]

  out[, class_route := data.table::fcase(
    isTRUE_vec(flt_final),                                   "phase60_final",
    isTRUE_vec(phase77_strict) & translation_class == "B",   "phase77_classB",
    isTRUE_vec(phase77_strict) & translation_class == "A",   "phase77_classA",
    isTRUE_vec(panel_member),                                "panel_only",
    isTRUE_vec(resectable_marker),                           "resectable_only",
    n_phase60_pass > 0,                                      "partial_phase60",
    default                                                  = "no_pass")]

  # Reorder by pass intensity
  data.table::setorder(out, -n_phase60_pass, na.last = TRUE)
  out[]
}
