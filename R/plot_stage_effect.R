#' Forest plot of per-stage effect size with 95% CI
#'
#' Draws log2FC +/- 1.96 x lfcSE for each of the 4 stages
#' (Normal, Early, Mid, Late), with significance asterisks at
#' padj < 0.05.
#'
#' @param gene_symbol HGNC symbol.
#' @return A `ggplot2` object.
#' @examples
#'   plot_stage_effect("LTBP1")
#' @export
plot_stage_effect <- function(gene_symbol) {
  d <- query_gene_detailed(gene_symbol)
  if (is.null(d)) return(invisible(NULL))
  ps <- d$per_stage

  ps$sig_label <- ifelse(ps$significant, "*", "")

  ggplot2::ggplot(ps, ggplot2::aes(x = stage, y = log2FC)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60",
                          linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
      width = 0.18, linewidth = 0.45, color = "#37474F") +
    ggplot2::geom_point(size = 2.5, color = "#0D47A1") +
    ggplot2::geom_text(ggplot2::aes(y = ci_hi + 0.05,
                                        label = sig_label),
                         color = "#C62828", size = 4,
                         fontface = "bold", vjust = 0) +
    ggplot2::labs(
      title = sprintf("%s - per-stage log2FC (mean +/- 95%% CI)",
                          gene_symbol),
      subtitle = "* Wald padj < 0.05 (per-contrast vs Normal, cohort-adjusted)",
      x = NULL, y = "log2FC vs Normal") +
    pdactrace_axes_theme()
}
