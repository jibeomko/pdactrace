#' Visualise tissue-to-serum filter passage with serum direction context
#'
#' For each gene, draws the 7-step phase60 filter trail (filled =
#' pass, open = fail, grey = not evaluated) with two extra layers
#' that the v0.99.6 baseline lacked:
#'
#' \itemize{
#'   \item A **serum direction strip** above the filter grid showing
#'     the per-gene serum log2FC vs healthy control and vs
#'     pancreatitis (the values that actually drive the
#'     `direction-match` filter step), coloured by translation
#'     class (Class A blue / Class B red / Class C grey / NA
#'     light-grey).
#'   \item A **per-gene title prefix** naming the translation
#'     class so the reader can read the trail in context.
#' }
#'
#' Multi-gene calls show one row per gene, ordered by phase60
#' pass count descending; this lets you compare `phase60_final`
#' exemplars (e.g. SERPINA1) against `phase77` Class B exemplars
#' (e.g. LTBP1, CDH13, SPARC) at a glance.
#'
#' @param gene Character vector of HGNC symbols.
#' @param show_routes Logical. If `TRUE` (default), append the
#'   right-side panel showing `phase77_strict` / `panel_member` /
#'   `resectable_marker` route flags.
#' @param show_serum Logical. If `TRUE` (default), append the
#'   serum-direction strip above the filter grid.
#' @return A `ggplot2` object (or `patchwork` composite if
#'   `show_routes = TRUE` or `show_serum = TRUE`).
#' @examples
#'   plot_filter_trace(c("SERPINA1", "LTBP1", "SPARC", "CDH13"))
#'   plot_filter_trace("LTBP1", show_routes = FALSE)
#' @export
plot_filter_trace <- function(gene,
                               show_routes = TRUE,
                               show_serum = TRUE) {
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

  filter_block <- patchwork::wrap_plots(p_left, p_right) +
    patchwork::plot_layout(widths = c(2.5, 1))

  if (!isTRUE(show_serum)) return(filter_block)

  # -- Serum direction strip (shows the values that drive the
  # 'direction-match' filter step + per-gene Class A/B/C) ----
  serum_block <- .ft_serum_strip(unique(as.character(tf$gene_symbol)))
  if (is.null(serum_block)) return(filter_block)

  patchwork::wrap_plots(serum_block, filter_block, ncol = 1L) +
    patchwork::plot_layout(heights = c(1.2, 2.5))
}

# Build the small per-gene serum direction strip that sits above
# the filter grid in plot_filter_trace().
.ft_serum_strip <- function(genes) {
  ref <- .get_reference()
  target <- genes
  rows <- ref[gene_symbol %in% target,
              .(gene_symbol,
                serum_log2fc_PDAC_vs_HC,
                serum_log2fc_Pan_vs_HC,
                translation_class)]
  if (nrow(rows) == 0L) return(NULL)
  long <- data.table::melt(
    rows, id.vars = c("gene_symbol", "translation_class"),
    measure.vars = c("serum_log2fc_PDAC_vs_HC",
                      "serum_log2fc_Pan_vs_HC"),
    variable.name = "contrast", value.name = "log2fc",
    variable.factor = FALSE)
  long[, contrast := factor(contrast,
    levels = c("serum_log2fc_PDAC_vs_HC",
               "serum_log2fc_Pan_vs_HC"),
    labels = c("PDAC vs HC", "Pancreatitis vs HC"))]
  long[, gene_symbol := factor(gene_symbol, levels = genes)]
  long[, class_lbl := data.table::fcase(
    is.na(translation_class), "no serum data",
    translation_class == "A",  "Class A (same dir.)",
    translation_class == "B",  "Class B (opposite)",
    translation_class == "C",  "Class C (decoupled)")]
  long[, class_lbl := factor(class_lbl, levels = c(
    "Class A (same dir.)", "Class B (opposite)",
    "Class C (decoupled)", "no serum data"))]
  pal <- c("Class A (same dir.)" = "#1565C0",
           "Class B (opposite)"  = "#C62828",
           "Class C (decoupled)" = "#616161",
           "no serum data"       = "grey85")
  ggplot2::ggplot(long,
    ggplot2::aes(x = log2fc, y = gene_symbol, fill = class_lbl)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey60",
                         linewidth = 0.3, linetype = "dashed") +
    ggplot2::geom_col(width = 0.55) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.na(log2fc), "NA",
                                   sprintf("%+.2f", log2fc))),
      size = 2.5, hjust = -0.15, color = "grey25") +
    ggplot2::facet_wrap(~ contrast, ncol = 2L, scales = "free_x") +
    ggplot2::scale_fill_manual(values = pal, name = NULL,
                                drop = FALSE) +
    ggplot2::labs(
      title = "Serum direction context",
      subtitle = paste0(
        "Inputs to the `direction-match` filter step: ",
        "PDAC log2FC vs HC and Pancreatitis log2FC vs HC, ",
        "coloured by translation class."),
      x = "log2FC", y = NULL) +
    pdactrace_axes_theme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 6, face = "bold"),
      strip.text  = ggplot2::element_text(size = 6, face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom")
}
