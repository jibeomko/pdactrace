#' Per-gene template-overlay plot
#'
#' For one gene/protein, draws the trajectory-template panel that
#' the gene was matched to (12-template Pearson-rho argmax) with the
#' gene's own z-scored Normal/Early/Mid/Late trajectory overlaid as
#' a highlight line. Companion to the atlas-wide
#' [plot_template_atlas()]: that one walks all 12 templates without
#' a gene; this one zeroes in on a single gene's matched panel.
#'
#' The matched template is recovered via the same Pearson-rho argmax
#' across the 12 bundled templates that backs [classify_trajectory()]
#' and [score_trajectory()]. The matched rho is annotated as a small
#' caption.
#'
#' @param gene_symbol HGNC gene symbol (length-1 character).
#' @param layer Either `"rna"` (default) or `"protein"`.
#' @param highlight_color Colour of the gene-highlight line. Default
#'   `"#C62828"` (NCS-grade muted red).
#' @param output_file Optional output PDF path. If `NULL` (default),
#'   the function returns a `ggplot` object; otherwise it writes one
#'   cairo PDF via [pdactrace_save()] and returns the path invisibly.
#' @param width,height PDF size in inches when `output_file` is set.
#'   Defaults to the compact 1.55 x 1.40 fig2C style.
#' @param reference Optional override of the bundled
#'   `pdactrace_reference` (RNA layer) for unit-test injection.
#' @return Either a `ggplot` object (when `output_file` is `NULL`) or
#'   the absolute path of the written PDF (invisibly).
#' @examples
#' p <- plot_gene_template("LGALS3BP", "rna")
#' class(p)
#' p$labels$title
#' @export
plot_gene_template <- function(gene_symbol,
                                 layer = c("rna", "protein"),
                                 highlight_color = "#C62828",
                                 output_file = NULL,
                                 width  = 1.55,
                                 height = 1.40,
                                 reference = NULL) {
  if (!is.character(gene_symbol) || length(gene_symbol) != 1L ||
      !nzchar(gene_symbol)) {
    stop("`gene_symbol` must be a length-1 non-empty character.",
         call. = FALSE)
  }
  layer <- match.arg(layer)

  target <- gene_symbol
  asg <- .assign_templates_12(layer, reference = reference)
  hit <- asg[gene_symbol == target]
  if (nrow(hit) == 0L) {
    stop(sprintf(
      "Gene '%s' is not in the %s atlas.", gene_symbol, layer),
      call. = FALSE)
  }
  template <- hit$template_argmax[1L]
  rho      <- hit$rho_argmax[1L]
  if (is.na(template)) {
    stop(sprintf(
      "Gene '%s' has no usable per-stage profile in the %s layer ",
      gene_symbol, layer),
      "(SD = 0 or NAs in beta_*).", call. = FALSE)
  }

  is_non_early <- !template %in% early_pattern_names()
  subtitle <- if (is_non_early) {
    "non-Early best-match (not surfaced in atlas rna_pattern)"
  } else {
    NULL
  }
  rho_caption <- sprintf(
    "%s | matched %s (rho = %.2f)",
    gene_symbol, template, rho)

  agg <- .template_aggregate(layer, template, reference = reference)
  p <- .plot_template_panel(
    agg,
    highlight_gene  = gene_symbol,
    highlight_color = highlight_color,
    template_label  = rho_caption,
    subtitle        = subtitle)

  if (!is.null(output_file)) {
    out_dir  <- dirname(output_file)
    out_stem <- sub("\\.pdf$", "", basename(output_file))
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    fp <- pdactrace_save(p, dir = out_dir, name = out_stem,
                          w = width, h = height)
    return(invisible(fp))
  }
  p
}
