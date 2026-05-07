#' Compare a gene panel across evidence layers as a heatmap
#'
#' Builds a multi-row heatmap (genes x evidence columns) showing key
#' RNA / protein / serum / class / route evidence side-by-side. Useful
#' for paper Figure 7 (Class A vs Class B comparison).
#'
#' @param genes Character vector of HGNC symbols.
#' @param layers Character vector of evidence groups to include.
#'   Default: `c("rna", "protein", "serum", "class")`.
#' @return A `ggplot2` object.
#' @examples
#'   plot_panel_heatmap(c("LTBP1", "SERPINA1", "CDH13", "SPARC", "CP", "FGB"))
#' @export
plot_panel_heatmap <- function(genes,
                                 layers = c("rna", "protein", "serum",
                                             "class")) {
  qp <- query_panel(genes,
                     layers = c("rna", "protein", "scrna",
                                  "serum", "clinical", "filter_status"))
  if (nrow(qp) == 0L) return(invisible(NULL))

  # Build a long-format evidence table - each row = (gene, variable, value)
  rows <- list()

  if ("rna" %in% layers) {
    rna_long <- data.table::data.table(
      gene_symbol = qp$gene_symbol,
      variable    = "RNA pattern",
      label       = qp$rna_pattern,
      value       = qp$rna_pattern_rho)
    rows$rna <- rna_long
  }
  if ("protein" %in% layers) {
    rows$prot <- data.table::data.table(
      gene_symbol = qp$gene_symbol,
      variable    = "Tissue protein",
      label       = qp$prot_pattern,
      value       = ifelse(qp$rnaprot_concordant, 1, 0))
  }
  if ("serum" %in% layers) {
    rows$serum <- data.table::data.table(
      gene_symbol = qp$gene_symbol,
      variable    = "Serum log2FC PDAC vs HC",
      label       = ifelse(is.na(qp$serum_log2fc_PDAC_vs_HC), "-",
                            sprintf("%.2f", qp$serum_log2fc_PDAC_vs_HC)),
      value       = qp$serum_log2fc_PDAC_vs_HC)
  }
  if ("class" %in% layers) {
    rows$class <- data.table::data.table(
      gene_symbol = qp$gene_symbol,
      variable    = "Translation class",
      label       = qp$translation_class,
      value       = ifelse(is.na(qp$translation_class), 0,
                            ifelse(qp$translation_class == "A", 1,
                              ifelse(qp$translation_class == "B", -1,
                                ifelse(qp$translation_class == "C", 0.5, 0)))))
  }

  long <- data.table::rbindlist(rows, fill = TRUE)
  long[, gene_symbol := factor(gene_symbol, levels = rev(genes))]
  long[, variable := factor(variable,
    levels = c("RNA pattern", "Tissue protein",
                "Serum log2FC PDAC vs HC", "Translation class"))]

  ggplot2::ggplot(
    long,
    ggplot2::aes(x = variable, y = gene_symbol)) +
    ggplot2::geom_tile(ggplot2::aes(fill = value),
                         color = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = label),
                         size = 1.9, color = "grey15",
                         fontface = "bold") +
    ggplot2::scale_fill_gradient2(
      low = "#1565C0", mid = "white", high = "#C62828",
      midpoint = 0, na.value = "grey92",
      name = "value") +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(title = "Panel evidence heatmap",
                    x = NULL, y = NULL) +
    pdactrace_axes_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 0,
                                              size = 5.5),
      axis.text.y = ggplot2::element_text(size = 6, face = "bold"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right")
}
