#' Visualise filter passage as a step bar chart
#'
#' For each gene, draws a 7-step horizontal bar showing which phase60
#' filters passed (filled) vs failed (open) vs not evaluated (grey),
#' alongside auxiliary route flags (phase77 strict, panel member,
#' resectable marker). Designed for direct comparison between
#' `phase60_final` exemplars (e.g. SERPINA1) and `phase77_classB`
#' exemplars (e.g. LTBP1, CDH13, SPARC).
#'
#' @param gene Character vector of HGNC symbols.
#' @param show_routes Logical. If `TRUE` (default), append the
#'   right-side panel showing phase77 strict / panel member /
#'   resectable marker flags.
#' @return A `ggplot2` object (or patchwork composite if
#'   `show_routes = TRUE`).
#' @examples
#'   plot_filter_trace(c("SERPINA1", "LTBP1", "SPARC", "CDH13"))
#' @export
plot_filter_trace <- function(gene, show_routes = TRUE) {
  tf <- trace_filters(gene)
  if (nrow(tf) == 0L) return(invisible(NULL))

  flt_cols <- c("flt_signal_peptide", "flt_serum_measurable",
                 "flt_serum_significant", "flt_pancreatitis_pdac",
                 "flt_pancreatitis_hc", "flt_direction_match",
                 "flt_final")
  flt_labels <- c("SignalP", "serum-measurable", "serum-significant",
                   "pan vs PDAC", "HC-in-middle", "direction-match",
                   "FINAL")

  long <- data.table::melt(tf, id.vars = "gene_symbol",
                              measure.vars = flt_cols,
                              variable.name = "step",
                              value.name = "pass",
                              variable.factor = FALSE)
  long[, step := factor(step, levels = flt_cols, labels = flt_labels)]
  long[, status := data.table::fcase(
    is.na(pass),     "not evaluated",
    pass == TRUE,    "PASS",
    pass == FALSE,   "FAIL")]
  long[, status := factor(status,
    levels = c("PASS", "FAIL", "not evaluated"))]
  # Preserve gene order from trace_filters (rank by n_phase60_pass desc)
  long[, gene_symbol := factor(gene_symbol,
                                  levels = tf$gene_symbol)]

  pal <- c(PASS = "#2E7D32", FAIL = "#FFFFFF",
            "not evaluated" = "grey92")

  p_left <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = step, y = gene_symbol, fill = status)) +
    ggplot2::geom_tile(color = "grey50", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = pal, name = NULL) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(
      title = "phase60 7-step filter trail",
      subtitle = "Filled = pass, open = fail, grey = not evaluated",
      x = NULL, y = NULL) +
    pdactrace_axes_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 0,
                                              size = 5.5),
      axis.text.y = ggplot2::element_text(size = 6, face = "bold"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom")

  if (!show_routes) return(p_left)

  # -- Right panel: route flags -------------------------------
  route_long <- data.table::melt(
    tf, id.vars = "gene_symbol",
    measure.vars = c("phase77_strict", "panel_member", "resectable_marker"),
    variable.name = "flag", value.name = "value",
    variable.factor = FALSE)
  route_long[, flag := factor(flag,
    levels = c("phase77_strict", "panel_member", "resectable_marker"),
    labels = c("phase77 strict", "panel member", "resectable marker"))]
  route_long[, status := data.table::fcase(
    is.na(value),    "not evaluated",
    value == TRUE,   "FLAG",
    value == FALSE,  "-")]
  route_long[, status := factor(status,
    levels = c("FLAG", "-", "not evaluated"))]
  route_long[, gene_symbol := factor(gene_symbol,
                                        levels = tf$gene_symbol)]

  pal_r <- c(FLAG = "#0D47A1", "-" = "#FFFFFF",
              "not evaluated" = "grey92")

  p_right <- ggplot2::ggplot(
    route_long,
    ggplot2::aes(x = flag, y = gene_symbol, fill = status)) +
    ggplot2::geom_tile(color = "grey50", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = pal_r, name = NULL) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(
      title = "Auxiliary routes",
      subtitle = NULL,
      x = NULL, y = NULL) +
    pdactrace_axes_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 0,
                                              size = 5.5),
      axis.text.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom")

  patchwork::wrap_plots(p_left, p_right) +
    patchwork::plot_layout(widths = c(2.5, 1))
}
