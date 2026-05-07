#' Compare multiple candidate genes side-by-side
#'
#' Returns a single `data.table` row-per-gene with the columns most
#' useful for panel-design discussion: audit class + score, RNA and
#' protein trajectory pattern, translation class, dominant scRNA cell
#' origin, serum-detectability flags, and a redundancy hint based on
#' shared cell origin and shared RNA pattern. The result is sorted by
#' `audit_score` (descending) so the strongest candidates surface
#' first.
#'
#' Input genes that are not in the bundled atlas are returned with
#' `audit_class = NA`, so the call never fails on stray symbols — the
#' returned table makes the gap explicit.
#'
#' @param gene_symbols Character vector of HGNC gene symbols
#'   (length >= 2; length-1 input is allowed but the comparison is
#'   trivial).
#' @return A `data.table` (one row per gene) with columns:
#'   * `gene_symbol`, `audit_class`, `audit_score`
#'   * `rna_pattern`, `prot_pattern`, `translation_class`
#'   * `cell_origin_top`, `serum_detected`,
#'     `serum_log2fc_PDAC_vs_HC`
#'   * `max_I2_meta`, `redundancy_with` — the latter lists other genes
#'     in the input set that share `rna_pattern` *and* `cell_origin_top`
#'     (semicolon-separated, "" if none).
#' @examples
#' compare_candidates(c("LGALS3BP", "LTBP1", "SERPINA1", "ALB", "GAPDH"))
#' @export
compare_candidates <- function(gene_symbols) {
  if (!is.character(gene_symbols) || length(gene_symbols) == 0L) {
    stop("`gene_symbols` must be a non-empty character vector.",
         call. = FALSE)
  }
  ref <- .get_reference()
  cols <- c("gene_symbol", "audit_class", "audit_score",
            "rna_pattern", "prot_pattern", "translation_class",
            "cell_origin_top", "serum_detected",
            "serum_log2fc_PDAC_vs_HC", "max_I2_meta")
  cols <- intersect(cols, names(ref))
  out <- ref[gene_symbol %in% gene_symbols,
             ..cols]

  # Add input genes that are not in the atlas as NA-rows
  missing_genes <- setdiff(gene_symbols, out$gene_symbol)
  if (length(missing_genes) > 0L) {
    pad <- data.table::data.table(gene_symbol = missing_genes)
    out <- data.table::rbindlist(list(out, pad), fill = TRUE)
  }

  # Redundancy hint: same rna_pattern AND same cell_origin_top -----
  out[, redundancy_with := vapply(seq_len(.N), function(i) {
    same <- which(
      seq_len(.N) != i &
      !is.na(rna_pattern) &
      !is.na(rna_pattern[i]) &
      rna_pattern == rna_pattern[i] &
      !is.na(cell_origin_top) &
      !is.na(cell_origin_top[i]) &
      cell_origin_top == cell_origin_top[i])
    if (length(same) == 0L) "" else
      paste(gene_symbol[same], collapse = "; ")
  }, character(1))]

  data.table::setorderv(out, "audit_score", order = -1L, na.last = TRUE)
  out[]
}
