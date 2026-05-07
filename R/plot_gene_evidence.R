#' Multi-panel ggplot composite of gene evidence
#'
#' Builds a publication-ready composite (`patchwork`) showing:
#'   * 4-point trajectory (RNA + Protein overlay)
#'   * Cell origin distribution bar
#'   * Serum log2FC (PDAC, Pancreatitis vs HC)
#'   * Provenance + summary footer
#'
#' Width/height: NCS spec single column (NCS_W_SINGLE = 3.46 in) by ~3 in.
#'
#' @param gene_symbol HGNC standard symbol.
#' @param layers Subset of panels to draw. Default all four
#'   (`c("trajectory", "cell_origin", "serum", "summary")`).
#' @return A `patchwork` object. Returns `NULL` invisibly with message
#'   if gene is not in atlas.
#' @examples
#'   p <- plot_gene_evidence("LTBP1")
#'   print(p)
#' @export
plot_gene_evidence <- function(gene_symbol,
                                layers = c("trajectory", "cell_origin",
                                            "serum", "summary")) {
  q <- query_gene(gene_symbol)
  if (is.null(q)) return(invisible(NULL))

  panels <- list()

  # -- Panel 1: Trajectory (RNA beta + Protein beta) ----------------
  if ("trajectory" %in% layers) {
    rna_traj <- data.table::data.table(
      stage = factor(c("Normal", "Early", "Mid", "Late"),
                      levels = c("Normal", "Early", "Mid", "Late")),
      value = c(q$rna$beta_N, q$rna$beta_E, q$rna$beta_M, q$rna$beta_L),
      layer = "RNA")
    panels$trajectory <- ggplot2::ggplot(
      rna_traj,
      ggplot2::aes(x = stage, y = value, group = layer)) +
      ggplot2::geom_line(linewidth = 0.5, color = "#0D47A1") +
      ggplot2::geom_point(size = 2.0, color = "#0D47A1") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                           color = "grey60", linewidth = 0.3) +
      ggplot2::labs(
        title = sprintf("%s : %s",
                          gene_symbol,
                          q$rna$pattern %||% "Unclassified"),
        subtitle = sprintf("rho = %.2f, LRT padj = %s",
                              q$rna$pattern_rho,
                              fmt_p(q$rna$lrt_padj)),
        x = NULL, y = "log2FC vs Normal") +
      pdactrace_axes_theme()
  }

  # -- Panel 2: Cell origin distribution ----------------------
  if ("cell_origin" %in% layers && !is.na(q$scrna$top_celltype)) {
    distrib <- q$scrna$distribution[[1]]
    pct <- 100 * distrib / sum(distrib)
    pct <- sort(pct, decreasing = TRUE)[1:min(6, length(pct))]
    co_dt <- data.table::data.table(
      celltype = factor(names(pct), levels = rev(names(pct))),
      pct      = as.numeric(pct))
    panels$cell_origin <- ggplot2::ggplot(
      co_dt,
      ggplot2::aes(x = celltype, y = pct)) +
      ggplot2::geom_col(fill = "#37474F", width = 0.7) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.0f%%", pct)),
        hjust = -0.2, size = 2.0) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
      ggplot2::labs(title = "Cell origin (scRNA)",
                     x = NULL, y = "% expression") +
      pdactrace_axes_theme() +
      ggplot2::coord_flip()
  }

  # -- Panel 3: Serum log2FC vs HC ----------------------------
  if ("serum" %in% layers && isTRUE(q$serum$detected) &&
      (!is.na(q$serum$log2fc_PDAC_vs_HC) ||
        !is.na(q$serum$log2fc_Pan_vs_HC))) {
    sr_dt <- data.table::data.table(
      group = factor(c("PDAC vs HC", "Pancreatitis vs HC"),
                      levels = c("PDAC vs HC", "Pancreatitis vs HC")),
      value = c(q$serum$log2fc_PDAC_vs_HC, q$serum$log2fc_Pan_vs_HC))
    sr_dt <- sr_dt[!is.na(sr_dt$value)]
    panels$serum <- ggplot2::ggplot(
      sr_dt,
      ggplot2::aes(x = group, y = value, fill = group)) +
      ggplot2::geom_col(width = 0.65) +
      ggplot2::geom_hline(yintercept = 0, color = "grey60",
                           linewidth = 0.3) +
      ggplot2::scale_fill_manual(
        values = c("PDAC vs HC" = "#E57373",
                    "Pancreatitis vs HC" = "#FFB74D"),
        guide = "none") +
      ggplot2::labs(
        title = sprintf("Serum log2FC (Class %s)",
                          q$serum$translation_class %||% "?"),
        x = NULL, y = "log2FC") +
      pdactrace_axes_theme()
  }

  if (length(panels) == 0L) {
    message(sprintf(
      "No data layers available for plot of '%s'.", gene_symbol))
    return(invisible(NULL))
  }

  # -- Compose with patchwork ---------------------------------
  composite <- patchwork::wrap_plots(panels, ncol = length(panels)) +
    patchwork::plot_annotation(
      title = q$summary,
      caption = sprintf("Provenance: %s | pdactrace v%s",
                          q$provenance, .pkg_version()),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 7, face = "bold",
                                              color = "grey10"),
        plot.caption = ggplot2::element_text(size = 5,
                                                color = "grey45")))
  composite
}
