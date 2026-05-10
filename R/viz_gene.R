#' Single-call visual evidence + scoring canvas for one gene
#'
#' Returns a 2x3 `patchwork` composite that puts **one panel per
#' evidence layer plus the final audit-axis radar** for a gene on
#' a single page (v0.99.8 expansion from the previous 4-panel
#' layout). The figure shows the entire prioritisation flow on one
#' page: four input layers (RNA, tissue protein, scRNA cell
#' origin, serum direction), the 7-step tissue-to-serum filter
#' trail that gates them, and the 6-axis hexagon that compresses
#' all of the above into the final audit-score components.
#'
#' \enumerate{
#'   \item Top-left -- **bulk RNA-seq**: per-stage log2FC forest
#'     (mean +/- 95% CI vs Normal).
#'   \item Top-middle -- **tissue proteomics**: per-stage protein
#'     log2FC trajectory (point estimates; per-stage SE not
#'     available from the bundled limma F-test fit).
#'   \item Top-right -- **scRNA cell origin**: 11-celltype
#'     expression distribution with the dominant lineage
#'     highlighted.
#'   \item Bottom-left -- **serum direction**: per-gene PDAC vs HC
#'     and Pancreatitis vs HC log2FC, coloured by translation
#'     class (Class A blue, Class B red, Class C grey).
#'   \item Bottom-middle -- **7-step filter trace**: the phase60
#'     tissue-to-serum audit trail with per-gene pass-count badge.
#'   \item Bottom-right -- **6-axis audit hexagon**: the final
#'     score's six component axes (multi-layer, direction,
#'     stage-onset, serum bridge, leakage safety, cohort
#'     consistency), the same axes [explain_score()] decomposes.
#' }
#'
#' A title strip across the top names the gene, its matched
#' template, audit class, and translation class. This single
#' figure is the visual analog of [summarize_gene_evidence()] +
#' [explain_score()] combined and is the recommended **first
#' call** for a gene a clinician or biologist hands you -- one
#' plot, the whole evidence-to-score chain, no per-axis function
#' names required.
#'
#' Internally composes [plot_stage_effect()] (RNA + protein),
#' [plot_celltype_full()], [plot_serum_direction()],
#' [plot_filter_trace()], and [plot_gene_hexagon()] via
#' [patchwork::wrap_plots()]. Every sub-panel remains available
#' standalone; the previous RNA per-cohort sign-vote view is at
#' [plot_per_cohort()].
#'
#' @param gene_symbol HGNC gene symbol (length-1 character).
#' @param title Optional title override. `NULL` (default) builds
#'   the title from the bundled atlas headline.
#' @param ncol Layout: `3` (default; 2x3), `2` (3x2), or `1`
#'   (6x1 vertical strip suitable for narrow embedding).
#' @return A `patchwork` object printable to any active graphics
#'   device or saveable via `ggsave()` / [pdactrace_save()].
#' @examples
#' viz_gene("LGALS3BP")
#' \donttest{
#'   # Alternative layouts (skipped in R CMD check examples for time;
#'   # all run normally in interactive use).
#'   viz_gene("LTBP1", ncol = 2)            # 3x2 grid
#'   viz_gene("LTBP1", ncol = 1)            # vertical strip
#'
#'   # Split mode: get each panel as its own full-size figure
#'   panels <- viz_gene("LTBP1", layout = "split")
#'   panels$rna       # bulk RNA per-stage forest
#'   panels$protein   # tissue protein per-stage trajectory
#'   panels$cell      # scRNA cell-of-origin distribution
#'   panels$serum     # serum direction strip
#'   panels$filter    # 7-step filter trail
#'   panels$hexagon   # 6-axis audit hexagon
#' }
#' \dontrun{
#'   # Or write each panel to a separate PDF in one call:
#'   viz_gene("LTBP1", layout = "split", output_dir = tempdir())
#' }
#' @param layout One of `"compact"` (default; the 2x3 patchwork
#'   composite) or `"split"` (returns a named list of six
#'   ggplot panels so the caller can render each as a full-size
#'   figure on its own page).
#' @param output_dir Optional directory path. Only used when
#'   `layout = "split"`. If supplied, each of the six panels is
#'   written to a separate cairo-PDF file via [pdactrace_save()]
#'   and the named list is returned invisibly with file paths
#'   attached as `attr(., "files")`.
#' @param width,height Per-panel width and height in inches when
#'   `layout = "split"` and `output_dir` is set. Defaults
#'   `6 x 4.5` (single-column, NCS-grade).
#' @seealso [summarize_gene_evidence()] for the text counterpart;
#'   [report_gene()] for an HTML report (requires pandoc).
#' @export
viz_gene <- function(gene_symbol,
                     title = NULL,
                     ncol = 3L,
                     layout = c("compact", "split"),
                     output_dir = NULL,
                     width = 6,
                     height = 4.5) {
  layout <- match.arg(layout)
  if (!is.character(gene_symbol) || length(gene_symbol) != 1L) {
    stop("`gene_symbol` must be a length-1 character string.",
         call. = FALSE)
  }
  if (layout == "compact" &&
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop("`patchwork` is required for viz_gene(layout=\"compact\"). ",
         "Install via install.packages(\"patchwork\"), or use ",
         "layout = \"split\".", call. = FALSE)
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

  # Build the 6 panels -- 4 input layers + filter trace + final
  # audit-axis radar. Every builder is wrapped so a NULL / error
  # return is replaced by a graceful "no data" placeholder; this
  # keeps the 6-panel grid intact for genes with sparse evidence.
  p_rna   <- .vg_safe(plot_stage_effect(gene_symbol, layer = "rna"),
                       gene_symbol, "RNA per-stage")
  p_prot  <- .vg_safe(plot_stage_effect(gene_symbol, layer = "protein"),
                       gene_symbol, "tissue protein per-stage")
  p_cell  <- .vg_safe(suppressWarnings(plot_celltype_full(gene_symbol)),
                       gene_symbol, "scRNA cell origin",
                       msg = "no scRNA cell-origin data")
  p_serum <- .vg_safe(suppressWarnings(plot_serum_direction(gene_symbol)),
                       gene_symbol, "serum direction",
                       msg = "no serum data")
  p_filt  <- .vg_safe(plot_filter_trace(gene_symbol,
                                          show_routes = FALSE,
                                          show_serum  = FALSE),
                       gene_symbol, "7-step filter trail",
                       msg = "no phase60 filter data")
  p_hex   <- .vg_safe(tryCatch(plot_gene_hexagon(gene_symbol),
                                  error = function(e) NULL),
                       gene_symbol, "audit hexagon",
                       msg = "audit components missing")

  if (layout == "split") {
    panels <- list(rna = p_rna, protein = p_prot, cell = p_cell,
                    serum = p_serum, filter = p_filt,
                    hexagon = p_hex)
    if (!is.null(output_dir)) {
      dir.create(output_dir, recursive = TRUE,
                  showWarnings = FALSE)
      files <- character(0L)
      for (nm in names(panels)) {
        fname <- sprintf("viz_gene_%s_%s", gene_symbol, nm)
        suppressMessages(pdactrace_save(panels[[nm]], output_dir,
                                         fname,
                                         w = width, h = height))
        files <- c(files, file.path(output_dir,
                                     paste0(fname, ".pdf")))
      }
      attr(panels, "files") <- files
      message(sprintf(
        "Wrote %d panels to %s", length(files), output_dir))
      return(invisible(panels))
    }
    return(panels)
  }

  # layout == "compact"
  composed <- patchwork::wrap_plots(p_rna, p_prot, p_cell,
                                     p_serum, p_filt, p_hex,
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

.vg_empty_panel <- function(title, msg) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 1, y = 1, label = msg,
                        color = "grey40", size = 3) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    pdactrace_axes_theme() +
    ggplot2::theme(axis.text  = ggplot2::element_blank(),
                    axis.ticks = ggplot2::element_blank())
}

# Wrap a panel-builder result so NULL / non-ggplot returns become
# a graceful "no data" placeholder; keeps the 6-panel grid intact
# for genes with sparse evidence layers.
.vg_safe <- function(p, gene_symbol, panel_label,
                      msg = "no data") {
  if (inherits(p, "ggplot")) return(p)
  .vg_empty_panel(sprintf("%s -- %s", gene_symbol, panel_label),
                   msg)
}
