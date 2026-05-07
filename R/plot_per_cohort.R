#' Per-cohort sign vote bar with Stouffer summary
#'
#' For each of the 4 RNA cohorts (TCGA, CPTAC, GSE224564, GSE79668),
#' shows the gene's trend direction (Increasing / Decreasing / Non-monotonic)
#' and whether the trend is monotonic. Headline annotation displays
#' Stouffer aggregate p-value and cohort agreement %.
#'
#' This plot reveals cohort divergence that the headline pattern
#' label may hide.
#'
#' @param gene_symbol HGNC symbol.
#' @return A `ggplot2` object.
#' @examples
#'   plot_per_cohort("LTBP1")  # reveals TCGAup but CPTACdown divergence
#' @export
plot_per_cohort <- function(gene_symbol) {
  d <- query_gene_detailed(gene_symbol)
  if (is.null(d)) return(invisible(NULL))
  pc <- d$per_cohort
  if (nrow(pc) == 0L) return(invisible(NULL))

  pc[, trend_factor := factor(trend,
    levels = c("Increasing", "Non-monotonic", "Decreasing"))]
  pc[, cohort := factor(cohort, levels = rev(cohort))]

  pal <- c(Increasing       = "#C62828",
            "Non-monotonic"  = "grey75",
            Decreasing       = "#1565C0")

  stouffer_p <- attr(pc, "stouffer_p")
  agreement  <- attr(pc, "agreement_pct")
  subt <- sprintf("Stouffer p = %s . %s%% cohort agreement",
                    if (is.na(stouffer_p)) "NA" else
                      formatC(stouffer_p, format = "g", digits = 2),
                    round(100 * agreement))

  ggplot2::ggplot(pc, ggplot2::aes(x = trend_factor, y = cohort,
                                       fill = trend_factor,
                                       shape = monotonic)) +
    ggplot2::geom_point(size = 4.5, stroke = 0.6, color = "grey20") +
    ggplot2::scale_fill_manual(values = pal, name = "Trend") +
    ggplot2::scale_shape_manual(
      values = c(`TRUE` = 21, `FALSE` = 4),
      labels = c(`TRUE` = "monotonic", `FALSE` = "non-monotonic"),
      name = NULL) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::labs(
      title = sprintf("%s - per-cohort RNA trend", gene_symbol),
      subtitle = subt,
      x = NULL, y = NULL) +
    pdactrace_axes_theme() +
    ggplot2::theme(legend.position = "bottom",
                     legend.box = "vertical")
}
