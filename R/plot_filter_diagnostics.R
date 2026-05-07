#' Filter diagnostic step bar with underlying metric annotation
#'
#' For each of the 7 phase60 filter steps, shows pass/fail status as
#' a horizontal bar and annotates each row with the underlying metric
#' (e.g., `pool padj = 0.32`) so the user can see *why* a step failed.
#'
#' Distinct from `plot_filter_trace()` which compares multiple genes
#' across PASS/FAIL only. `plot_filter_diagnostics()` focuses on a
#' single gene with full numeric audit.
#'
#' @param gene_symbol HGNC symbol.
#' @return A `ggplot2` object.
#' @examples
#'   plot_filter_diagnostics("LTBP1")
#' @export
plot_filter_diagnostics <- function(gene_symbol) {
  d <- query_gene_detailed(gene_symbol)
  if (is.null(d)) return(invisible(NULL))
  fd <- d$filter_diag
  if (nrow(fd) == 0L) return(invisible(NULL))

  fd[, step := factor(step, levels = rev(step))]
  fd[, status := data.table::fcase(
    is.na(pass),  "not evaluated",
    pass == TRUE, "PASS",
    pass == FALSE, "FAIL")]
  fd[, status := factor(status, levels = c("PASS", "FAIL", "not evaluated"))]
  fd[, x_pos := 1]

  pal <- c(PASS = "#2E7D32", FAIL = "#C62828",
            "not evaluated" = "grey75")

  ggplot2::ggplot(fd, ggplot2::aes(x = x_pos, y = step,
                                       fill = status)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4,
                         width = 1.0, height = 0.7) +
    ggplot2::geom_text(ggplot2::aes(x = 1.6,
                                        label = underlying_metric),
                         hjust = 0, size = 1.9, color = "grey15") +
    ggplot2::geom_text(ggplot2::aes(x = 1, label = status),
                         color = "white", size = 1.9,
                         fontface = "bold") +
    ggplot2::scale_fill_manual(values = pal, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0.4, 6),
                                  expand = c(0, 0)) +
    ggplot2::labs(
      title = sprintf("%s - phase60 filter diagnostics", gene_symbol),
      subtitle = "Pass/fail with underlying metric (cutoff in tooltip)",
      x = NULL, y = NULL) +
    pdactrace_axes_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     legend.position = "bottom")
}
