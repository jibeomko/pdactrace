#' Hexagonal evidence radar for one or more genes
#'
#' Draws a compact 6-axis hexagonal radar chart summarising the
#' pdactrace audit-feature components for one or more genes. Inspired
#' by the multi-polygon radar style used in graph-of-graphs analyses
#' (e.g., GATE; Liu and Cheng 2025), but reuses the
#' six interpretable feature axes from `extract_graph_features()`
#' rather than learned graph embeddings.
#'
#' The six axes (clockwise from the top vertex) are:
#'
#' 1. **Multi-layer** — `score_layer` (fraction of RNA/Protein/scRNA/Serum
#'    layers with a measurable signal).
#' 2. **Direction** — `score_direction` (cross-cohort + cross-layer
#'    direction concordance).
#' 3. **Stage-onset** — `score_early` (Pearson rho against the four
#'    Early templates, weighted by LRT significance).
#' 4. **Serum bridge** — `score_serum` (serum detected + signal
#'    peptide + tissue-to-serum direction concordance).
#' 5. **Leakage safety** — `leakage_mult` (1 if clean, 0.5 if plasma
#'    high-abundance, 0 if housekeeping).
#' 6. **Cohort consistency** — `het_mult` (1 if I² < 70, 0.7 if 70-90,
#'    0.3 if >= 90).
#'
#' Each axis is plotted on a fixed 0--1 scale; concentric hexagons
#' mark the 0.25 / 0.50 / 0.75 / 1.00 contours.
#'
#' @param gene_symbol Character vector of one or more HGNC gene
#'   symbols. Multiple genes are overlaid as separate polygons.
#' @param comparison Optional. Either `NULL` (default), a single
#'   gene symbol, or one of `"high_confidence_mean"`,
#'   `"supported_uncertain_mean"`, `"penalized_mean"`,
#'   `"excluded_mean"`. Adds a translucent reference polygon for
#'   visual context.
#' @param show_score Logical. If `TRUE` (default), annotate the
#'   centre with each gene's `audit_score` and `audit_class`.
#' @param palette Optional named character vector of polygon colors
#'   keyed by gene symbol. Defaults to the pdactrace audit-class
#'   palette where possible, otherwise to a Set1-style palette.
#' @return A `ggplot2` object (NCS-grade compact, cairo-pdf safe).
#' @examples
#' plot_gene_hexagon("LTBP1")
#' plot_gene_hexagon(c("LGALS3BP", "LTBP1"))
#' plot_gene_hexagon("LTBP1", comparison = "high_confidence_mean")
#' @export
plot_gene_hexagon <- function(gene_symbol,
                                comparison = NULL,
                                show_score = TRUE,
                                palette = NULL) {
  if (!is.character(gene_symbol) || length(gene_symbol) == 0L) {
    stop("gene_symbol must be a non-empty character vector.")
  }

  axes <- c("Multi-layer", "Direction", "Stage-onset",
            "Serum bridge", "Leakage safety", "Cohort consistency")
  feat_cols <- c("score_layer", "score_direction", "score_early",
                  "score_serum", "leakage_mult", "het_mult")

  feats <- extract_graph_features(gene_symbol)
  if (is.null(feats) || nrow(feats) == 0L) {
    message("No evidence available for the supplied gene_symbol(s).")
    return(invisible(NULL))
  }

  audit <- compute_audit_score(gene_symbol)
  feats <- merge(feats[, c("gene_symbol", feat_cols), with = FALSE],
                 audit[, list(gene_symbol, audit_score, audit_class)],
                 by = "gene_symbol", all.x = TRUE)
  data.table::setkey(feats, NULL)

  # Comparison polygon ---------------------------------------------
  cmp_row <- NULL
  cmp_label <- NULL
  cmp_color <- "#9E9E9E"
  if (!is.null(comparison)) {
    class_means <- c("high_confidence_mean", "supported_uncertain_mean",
                     "penalized_mean", "excluded_mean")
    if (length(comparison) == 1L && comparison %in% class_means) {
      target_class <- sub("_mean$", "", comparison)
      ref <- .get_reference()
      ref_audit <- compute_audit_score(NULL)
      ref_feat <- extract_graph_features(ref$gene_symbol)
      ref_full <- merge(ref_feat[, c("gene_symbol", feat_cols),
                                  with = FALSE],
                        ref_audit[, list(gene_symbol, audit_class)],
                        by = "gene_symbol", all.x = TRUE)
      sub <- ref_full[audit_class == target_class]
      if (nrow(sub) > 0L) {
        cmp_row <- as.numeric(sapply(feat_cols, function(c)
          mean(sub[[c]], na.rm = TRUE)))
        cmp_label <- sprintf("%s mean (n=%d)",
                              target_class, nrow(sub))
        cmp_color <- c(high_confidence = "#1B5E20",
                        supported_uncertain = "#F57C00",
                        penalized = "#9E9E9E",
                        excluded = "#424242")[target_class]
      }
    } else if (is.character(comparison) && length(comparison) == 1L) {
      cmp_feat <- extract_graph_features(comparison)
      if (!is.null(cmp_feat) && nrow(cmp_feat) > 0L) {
        cmp_row <- as.numeric(cmp_feat[1L, feat_cols, with = FALSE])
        cmp_label <- comparison
        cmp_color <- "#0D47A1"
      }
    }
  }

  # Polygon coordinates --------------------------------------------
  n_ax <- length(axes)
  angles <- seq(pi / 2, pi / 2 - 2 * pi, length.out = n_ax + 1L)[1:n_ax]
  ring_levels <- c(0.25, 0.50, 0.75, 1.00)

  # Build concentric hexagonal rings
  ring_dt <- data.table::rbindlist(lapply(ring_levels, function(r) {
    data.table::data.table(
      ring = r,
      x = r * cos(c(angles, angles[1L])),
      y = r * sin(c(angles, angles[1L]))
    )
  }))

  axis_lines <- data.table::data.table(
    x_end = 1.05 * cos(angles),
    y_end = 1.05 * sin(angles),
    label_x = 1.18 * cos(angles),
    label_y = 1.18 * sin(angles),
    label = axes)

  # Per-gene polygon ------------------------------------------------
  AUDIT_COL <- c(high_confidence = "#1B5E20",
                 supported_uncertain = "#F57C00",
                 penalized = "#9E9E9E",
                 excluded = "#424242",
                 low = "#BDBDBD")
  DEFAULT_PAL <- c("#1B5E20", "#0D47A1", "#C62828", "#6A1B9A",
                   "#EF6C00", "#00838F")

  poly_list <- list()
  for (i in seq_len(nrow(feats))) {
    g <- feats$gene_symbol[i]
    vals <- pmax(0, pmin(1,
      as.numeric(feats[i, feat_cols, with = FALSE])))
    fill_col <- if (!is.null(palette) && g %in% names(palette)) {
      palette[[g]]
    } else if (!is.na(feats$audit_class[i])) {
      AUDIT_COL[feats$audit_class[i]]
    } else {
      DEFAULT_PAL[((i - 1) %% length(DEFAULT_PAL)) + 1]
    }
    poly_list[[g]] <- data.table::data.table(
      gene = g,
      x = c(vals * cos(angles), vals[1L] * cos(angles[1L])),
      y = c(vals * sin(angles), vals[1L] * sin(angles[1L])),
      color = fill_col)
  }
  poly_dt <- data.table::rbindlist(poly_list)

  cmp_poly <- NULL
  if (!is.null(cmp_row)) {
    cmp_row <- pmax(0, pmin(1, cmp_row))
    cmp_poly <- data.table::data.table(
      gene = cmp_label,
      x = c(cmp_row * cos(angles), cmp_row[1L] * cos(angles[1L])),
      y = c(cmp_row * sin(angles), cmp_row[1L] * sin(angles[1L])),
      color = cmp_color)
  }

  # Plot ------------------------------------------------------------
  p <- ggplot2::ggplot()
  # Concentric hexagonal rings
  for (r in ring_levels) {
    p <- p + ggplot2::geom_path(
      data = ring_dt[ring == r],
      ggplot2::aes(x = x, y = y),
      color = "#BDBDBD", linewidth = 0.25,
      linetype = if (r == 1.00) "solid" else "dashed")
  }
  # Axis spokes
  p <- p + ggplot2::geom_segment(
    data = axis_lines,
    ggplot2::aes(x = 0, y = 0, xend = x_end, yend = y_end),
    color = "#9E9E9E", linewidth = 0.3)
  # Axis labels (outside)
  p <- p + ggplot2::geom_text(
    data = axis_lines,
    ggplot2::aes(x = label_x, y = label_y, label = label),
    size = 2.2, fontface = "bold", color = "grey15")
  # Comparison polygon (translucent, behind gene polygons)
  if (!is.null(cmp_poly)) {
    p <- p + ggplot2::geom_polygon(
      data = cmp_poly,
      ggplot2::aes(x = x, y = y),
      fill = cmp_color, alpha = 0.12,
      color = cmp_color, linewidth = 0.4, linetype = "dotted")
  }
  # Gene polygons
  p <- p + ggplot2::geom_polygon(
    data = poly_dt,
    ggplot2::aes(x = x, y = y, group = gene, color = color),
    fill = NA, linewidth = 0.7) +
    ggplot2::geom_polygon(
      data = poly_dt,
      ggplot2::aes(x = x, y = y, group = gene, fill = color),
      alpha = 0.18, color = NA)
  # Center scoring annotation
  if (show_score) {
    if (!is.null(cmp_label)) {
      p <- p + ggplot2::annotate("text", x = 0, y = 1.42,
        label = sprintf("vs %s", cmp_label),
        size = 1.9, color = cmp_color, fontface = "italic")
    }
    # One annotation line per gene (color-coded), centered below hexagon
    n_g <- nrow(feats)
    line_y <- seq(from = -1.30, by = -0.12, length.out = n_g)
    line_color <- vapply(seq_len(n_g), function(i) {
      g <- feats$gene_symbol[i]
      if (!is.null(palette) && g %in% names(palette)) {
        palette[[g]]
      } else if (!is.na(feats$audit_class[i])) {
        AUDIT_COL[feats$audit_class[i]]
      } else {
        DEFAULT_PAL[((i - 1) %% length(DEFAULT_PAL)) + 1]
      }
    }, character(1))
    score_dt <- data.table::data.table(
      x = 0, y = line_y,
      label = sprintf("%s -- score %.2f (%s)",
                      feats$gene_symbol, feats$audit_score,
                      feats$audit_class),
      color = line_color)
    p <- p + ggplot2::geom_text(
      data = score_dt,
      ggplot2::aes(x = x, y = y, label = label, color = color),
      size = 2.1, fontface = "bold")
  }
  # Ring scale labels (0.25 / 0.50 / 0.75 / 1.00) on top spoke
  p <- p + ggplot2::geom_text(
    data = data.table::data.table(r = ring_levels,
                                    x = 0.04, y = ring_levels - 0.03,
                                    label = sprintf("%.2f", ring_levels)),
    ggplot2::aes(x = x, y = y, label = label),
    size = 1.7, color = "grey55", hjust = 0)

  # Pad ylim to fit multi-line annotation when many genes
  n_g_total <- nrow(feats)
  ymin <- if (n_g_total > 1L) -1.30 - 0.12 * n_g_total - 0.10 else -1.50
  p +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_fixed(xlim = c(-1.5, 1.5),
                          ylim = c(ymin, 1.5), expand = FALSE) +
    pdactrace_panel_theme()
}
