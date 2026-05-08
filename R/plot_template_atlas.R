#' Per-template trajectory atlas plots
#'
#' Renders one panel per trajectory template in the bundled
#' 12-template catalogue. Each panel shows the cohort of
#' genes / proteins whose z-scored Normal / Early / Mid / Late
#' trajectory matches that template best (Pearson rho argmax across
#' the 12 templates; deterministic given the bundled atlas), as thin
#' translucent lines plus a mean line and a +/-1 SD ribbon. Reproduces
#' the figures in
#' `analysis/manuscript/PDAC_졸업v1/fig2/individual/fig2C_RNA_*.pdf`
#' inside the package and extends them to the four templates that
#' did not exist in the v0.3.0 era (Late_Burst_Up, Late_Loss_Down,
#' Monotonic_Up, Monotonic_Down).
#'
#' Cohort assignment is recomputed inside the package via the same
#' z-scored Pearson-rho argmax that backs [classify_trajectory()] and
#' [score_trajectory()]; the bundled atlas's `rna_pattern` column
#' surfaces only the four Early templates, so the eight non-Early
#' panels rely on this argmax recovery rather than on a directly
#' bundled label.
#'
#' Edge case: very small cohorts degrade gracefully — singletons drop
#' the ribbon, cohorts of 2-4 drop the ribbon but keep individual
#' lines, cohorts of >=5 show the full thin-lines + mean + ribbon.
#'
#' @param layer Either `"rna"` (default) or `"protein"`. The protein
#'   layer reads the bundled [pdactrace_protein_betas] (5,917
#'   measurable proteins).
#' @param templates Optional character vector of template names to
#'   restrict the output to. Default `NULL` returns all 12.
#' @param output_dir Optional directory. If non-NULL, one cairo PDF
#'   per template is written via [pdactrace_save()] and the file
#'   paths are reported in `attr(out, "files")`. Default `NULL`
#'   returns ggplot objects only.
#' @param filename_pattern File-name template containing
#'   `{layer}` and `{template}` placeholders. Default
#'   `"template_{layer}_{template}.pdf"`.
#' @param width,height Output PDF size in inches. Defaults match the
#'   compact `fig2C` style (1.55 x 1.40 in; ~ 4 cm square).
#' @param reference Optional override of the bundled
#'   `pdactrace_reference` (RNA layer) for unit-test injection. The
#'   protein layer always reads `pdactrace_protein_betas`.
#' @return A named `list` of `ggplot` objects, one per template, in
#'   `c(early_pattern_names(), mid_pattern_names_excluded(),
#'   "Late_Burst_Up", "Late_Loss_Down", "Monotonic_Up",
#'   "Monotonic_Down")` order. When `output_dir` is non-NULL, the
#'   list also carries the file paths in `attr(out, "files")`.
#' @examples
#' panels <- plot_template_atlas("rna",
#'                                 templates = "Early_Burst_Up")
#' length(panels)
#' names(panels)
#' @export
plot_template_atlas <- function(layer = c("rna", "protein"),
                                  templates = NULL,
                                  output_dir = NULL,
                                  filename_pattern =
                                    "template_{layer}_{template}.pdf",
                                  width = 1.55,
                                  height = 1.40,
                                  reference = NULL) {
  layer <- match.arg(layer)
  all_templates <- c(
    early_pattern_names(),
    mid_pattern_names_excluded(),
    "Late_Burst_Up", "Late_Loss_Down",
    "Monotonic_Up", "Monotonic_Down")
  if (is.null(templates)) {
    templates <- all_templates
  } else {
    bad <- setdiff(templates, all_templates)
    if (length(bad) > 0L) {
      stop("Unknown template name(s): ",
           paste(bad, collapse = ", "), ". Valid: ",
           paste(all_templates, collapse = ", "), call. = FALSE)
    }
  }

  panels <- list()
  for (tpl in templates) {
    agg <- .template_aggregate(layer, tpl, reference = reference)
    panels[[tpl]] <- .plot_template_panel(
      agg,
      highlight_gene  = NULL,
      template_label  = tpl,
      subtitle        = NULL)
  }

  files <- character(0)
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    for (tpl in templates) {
      stem <- gsub("\\.pdf$", "", basename(
        gsub("\\{layer\\}", layer,
             gsub("\\{template\\}", tpl, filename_pattern))))
      files <- c(files,
                 pdactrace_save(panels[[tpl]],
                                 dir  = output_dir,
                                 name = stem,
                                 w = width, h = height))
    }
  }

  attr(panels, "files") <- files
  panels
}
