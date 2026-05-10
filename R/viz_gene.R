#' Single-call visual evidence canvas for one gene
#'
#' Returns a 2x2 `patchwork` composite that puts the four most
#' load-bearing evidence layers for one gene on a single page:
#'
#' \enumerate{
#'   \item Top-left: per-stage trajectory forest (log2FC +/- 95%
#'     CI across Normal / Early / Mid / Late) -- the *shape*.
#'   \item Top-right: per-cohort sign-vote bar with Stouffer
#'     summary -- *cohort robustness*.
#'   \item Bottom-left: scRNA cell-of-origin distribution across
#'     all 11 cell types -- *biological coherence*.
#'   \item Bottom-right: 7-step tissue-to-serum filter trace
#'     coloured by Class A / B / C translation outcome --
#'     *translational relevance*.
#' }
#'
#' A title strip across the top names the gene, its matched
#' template, audit class, and translation class. This single
#' figure is the visual analog of [summarize_gene_evidence()] and
#' is the recommended **first call** for a gene a clinician or
#' biologist hands you -- one plot, the whole evidence story, no
#' pre-existing knowledge of the per-axis function names required.
#'
#' Internally composes [plot_stage_effect()], [plot_per_cohort()],
#' [plot_celltype_full()], and [plot_filter_trace()] via
#' [patchwork::wrap_plots()]. All four sub-panels remain available
#' as standalone functions for users who want a single layer.
#'
#' @param gene_symbol HGNC gene symbol (length-1 character).
#' @param title Optional title override. `NULL` (default) builds
#'   the title from the bundled atlas headline.
#' @param ncol Layout: `2` (default; 2x2) or `1` (4x1 vertical
#'   strip suitable for narrow embedding).
#' @return A `patchwork` object printable to any active graphics
#'   device or saveable via `ggsave()` / [pdactrace_save()].
#' @examples
#' viz_gene("LGALS3BP")
#' viz_gene("LTBP1", ncol = 1)   # vertical strip
#' @seealso [summarize_gene_evidence()] for the text counterpart;
#'   [plot_gene_evidence()] for the older single-row composite;
#'   [report_gene()] for an HTML report (requires pandoc).
#' @export
viz_gene <- function(gene_symbol, title = NULL, ncol = 2L) {
  if (!is.character(gene_symbol) || length(gene_symbol) != 1L) {
    stop("`gene_symbol` must be a length-1 character string.",
         call. = FALSE)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("`patchwork` is required for viz_gene(). ",
         "Install via install.packages(\"patchwork\").",
         call. = FALSE)
  }

  # Resolve a one-line headline from the atlas.
  ref <- .get_reference()
  target <- gene_symbol
  row <- ref[gene_symbol == target]
  if (nrow(row) == 0L) {
    stop(sprintf("Gene '%s' is not in the bundled atlas.",
                 gene_symbol), call. = FALSE)
  }
  if (is.null(title)) title <- .vg_headline(row)

  # Build the 4 panels. Each returns a ggplot.
  p_traj <- plot_stage_effect(gene_symbol)
  p_coh  <- plot_per_cohort(gene_symbol)
  p_cell <- suppressWarnings(plot_celltype_full(gene_symbol))
  p_filt <- plot_filter_trace(gene_symbol, show_routes = FALSE)

  composed <- patchwork::wrap_plots(p_traj, p_coh, p_cell, p_filt,
                                     ncol = ncol) +
              patchwork::plot_annotation(
                title = title,
                theme = ggplot2::theme(
                  plot.title = ggplot2::element_text(
                    face = "bold", size = 11, hjust = 0)))
  composed
}

# ---- internal --------------------------------------------------------

.vg_headline <- function(row) {
  pat <- if (is.na(row$rna_pattern)) "no Early-onset pattern" else
    sprintf("%s (rho=%.2f)", row$rna_pattern, row$rna_pattern_rho)
  cls <- if (is.na(row$audit_class)) "no audit class" else
    as.character(row$audit_class)
  trc <- if (is.na(row$translation_class)) "no serum data" else
    paste0("Class ", row$translation_class)
  sprintf("%s -- %s | %s | %s",
          row$gene_symbol, pat, cls, trc)
}
