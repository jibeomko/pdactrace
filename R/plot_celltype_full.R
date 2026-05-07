#' Full 11-celltype expression bar
#'
#' Replaces the v0.1.0 `plot_gene_evidence()` top-3 cell-type panel
#' with the *full* 11-celltype distribution from the scRNA atlas.
#' Useful when the dominant cell type is borderline or when the user
#' wants to see the long tail (e.g., mast cells, plasma cells).
#'
#' @param gene_symbol HGNC symbol.
#' @param highlight Optional character vector of cell types to draw
#'   in red (default `c("myCAF", "iCAF", "Acinar", "Ductal")` -
#'   the four PDAC-relevant origins).
#' @return A `ggplot2` object.
#' @examples
#'   plot_celltype_full("LTBP1")
#' @export
plot_celltype_full <- function(gene_symbol,
                                 highlight = c("myCAF", "iCAF",
                                                 "Acinar", "Ductal")) {
  d <- query_gene_detailed(gene_symbol)
  if (is.null(d)) return(invisible(NULL))
  ct <- d$per_celltype
  if (nrow(ct) == 0L) return(invisible(NULL))

  ct[, celltype := factor(celltype, levels = rev(celltype))]
  ct[, is_highlight := celltype %in% highlight]
  tau <- attr(ct, "specificity_tau")
  subt <- if (is.null(tau) || is.na(tau)) "" else
    sprintf("Specificity tau = %.2f", tau)

  ggplot2::ggplot(ct, ggplot2::aes(x = pct_of_total, y = celltype,
                                       fill = is_highlight)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", pct_of_total)),
      hjust = -0.15, size = 1.9, color = "grey20") +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = "#C62828", `FALSE` = "grey55"),
      guide = "none") +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.20))) +
    ggplot2::labs(
      title = sprintf("%s - scRNA cell origin (full 11 types)",
                          gene_symbol),
      subtitle = subt,
      x = "% of expression", y = NULL) +
    pdactrace_axes_theme() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}
