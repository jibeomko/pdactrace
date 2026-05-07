#' Query the pdactrace reference atlas for multiple genes
#'
#' Convenience join of [query_gene()] outputs across multiple genes.
#'
#' @param genes Character vector of HGNC symbols.
#' @param joined Logical. If `TRUE` (default), returns a wide
#'   `data.table` with one row per gene and key columns from each
#'   layer. If `FALSE`, returns a long `data.table` with
#'   `gene_symbol x layer x variable x value` rows.
#' @param layers Character vector subset of available layers (passed
#'   through to [query_gene()]).
#' @return `data.table` (wide or long).
#' @examples
#'   query_panel(c("LTBP1", "SERPINA1", "CDH13"))
#' @export
query_panel <- function(genes,
                          joined = TRUE,
                          layers = c("rna", "protein", "scrna", "serum",
                                      "clinical", "filter_status",
                                      "annotation")) {
  if (length(genes) == 0L) return(data.table::data.table())
  genes <- unique(genes)

  ref <- .get_reference()
  hit <- ref[ref$gene_symbol %in% genes, ]
  miss <- setdiff(genes, hit$gene_symbol)
  if (length(miss) > 0L) {
    message(sprintf(
      "Genes not in atlas (skipped): %s",
      paste(miss, collapse = ", ")))
  }
  if (nrow(hit) == 0L) return(data.table::data.table())

  if (joined) {
    cols <- character()
    if ("rna"        %in% layers) cols <- c(cols,
      "rna_pattern", "rna_pattern_rho", "rna_lrt_padj", "rna_stouffer_z")
    if ("protein"    %in% layers) cols <- c(cols,
      "prot_pattern", "prot_tier", "rnaprot_concordant")
    if ("scrna"      %in% layers) cols <- c(cols,
      "cell_origin_top")
    if ("serum"      %in% layers) cols <- c(cols,
      "serum_detected", "translation_class", "serum_log2fc_PDAC_vs_HC")
    if ("clinical"   %in% layers) cols <- c(cols,
      "resectable_marker", "panel_member", "evidence_scope")
    if ("filter_status" %in% layers) cols <- c(cols,
      "flt_signal_peptide", "flt_serum_measurable",
      "flt_serum_significant", "flt_pancreatitis_pdac",
      "flt_pancreatitis_hc", "flt_direction_match", "flt_final")
    if ("annotation" %in% layers) cols <- c(cols,
      "ann_pool_logfc", "ann_pool_padj",
      "ann_pan_vs_hc_logfc", "ann_pan_vs_hc_pval",
      "ann_pdac_mean", "ann_pan_mean", "ann_hc_mean",
      "ann_pdac_vs_pan_pval")
    sel <- c("gene_symbol", cols)
    out <- hit[, sel, with = FALSE]
    # Reorder by user-supplied gene order (skip missing). setkey on
    # parent table makes integer i fall through to character lookup,
    # so we use base subset by integer position via .SD trick.
    idx <- match(genes, out$gene_symbol)
    idx <- idx[!is.na(idx)]
    return(out[idx, ])
  }

  # Long form
  rows <- vector("list", nrow(hit))
  for (i in seq_len(nrow(hit))) {
    g <- hit$gene_symbol[i]
    qg <- query_gene(g, layers = layers)
    if (is.null(qg)) next
    parts <- list()
    for (lyr in setdiff(names(qg), c("summary", "provenance"))) {
      df <- qg[[lyr]]
      if (is.null(df)) next
      ldf <- data.table::melt(df, id.vars = character(),
                                variable.name = "variable",
                                value.name = "value",
                                variable.factor = FALSE)
      ldf[, gene_symbol := g]
      ldf[, layer := lyr]
      parts[[lyr]] <- ldf[, .(gene_symbol, layer, variable, value)]
    }
    rows[[i]] <- data.table::rbindlist(parts, fill = TRUE)
  }
  data.table::rbindlist(rows, fill = TRUE)
}
