#' Visualise the candidate landscape across tissue and serum
#'
#' Scatter plot of atlas-surfaced Early × 4 genes (12-template classifier surface) with measurable
#' tissue + serum evidence. X-axis: tissue protein effect (max |beta|
#' across stages, signed by direction). Y-axis: serum log2FC PDAC vs
#' HC. Color: translation_class (A / B / C / unclassified). Selected
#' genes (default LTBP1, SERPINA1, CDH13, SPARC) get text labels.
#'
#' @param highlight Character vector of gene symbols to label
#'   prominently.
#' @param min_serum_n Integer. Drop genes detected in fewer cohorts.
#' @return A `ggplot2` object.
#' @examples
#'   plot_candidate_landscape()
#'   plot_candidate_landscape(highlight = c("LTBP1", "SERPINA1",
#'                                           "CDH13", "SPARC"))
#' @export
plot_candidate_landscape <- function(
    highlight = c("LTBP1", "SERPINA1", "CDH13", "SPARC"),
    min_serum_n = 1L) {
  ref <- .get_reference()

  d <- ref[!is.na(rna_pattern) &
            !is.na(serum_log2fc_PDAC_vs_HC) |
            phase77_strict == TRUE]

  # Tissue effect: signed max beta
  d[, tissue_effect := pmax(abs(rna_beta_E),
                              abs(rna_beta_M),
                              abs(rna_beta_L), na.rm = TRUE)]
  d[, tissue_dir := data.table::fcase(
    rna_pattern %in% c("Early_Burst_Up", "Early_Peak"),  +1,
    rna_pattern %in% c("Early_Loss_Down", "Early_Trough"), -1,
    default = 0)]
  d[, tissue_signed := tissue_effect * tissue_dir]

  # Serum effect - keep raw logFC; for highlighted genes that lack
  # raw values (LTBP1 has NA log2fc_PDAC_vs_HC), substitute phase77
  # vs_serum direction -> +/-0.5 placeholder
  d[is.na(serum_log2fc_PDAC_vs_HC) & phase77_strict == TRUE,
    serum_log2fc_PDAC_vs_HC := data.table::fcase(
      translation_class == "A",  +0.3,
      translation_class == "B",  -0.3,
      default = 0)]
  d <- d[!is.na(serum_log2fc_PDAC_vs_HC) & !is.na(tissue_signed)]

  d[, class_label := data.table::fcase(
    translation_class == "A",            "Class A (concordant)",
    translation_class == "B",            "Class B (inverse stromal)",
    translation_class == "C",            "Class C (decoupled)",
    default                              = "Other")]

  pal <- pdactrace_pal_class
  pal["Other"] <- "grey80"

  hl <- d[gene_symbol %in% highlight]

  ggplot2::ggplot(d, ggplot2::aes(x = tissue_signed,
                                       y = serum_log2fc_PDAC_vs_HC,
                                       color = class_label)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60",
                          linewidth = 0.3, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, color = "grey60",
                          linewidth = 0.3, linetype = "dashed") +
    ggplot2::geom_point(alpha = 0.55, size = 1.2) +
    ggplot2::geom_point(data = hl, size = 2.4, stroke = 0.7,
                          color = "black", shape = 21,
                          ggplot2::aes(fill = class_label)) +
    ggplot2::geom_text(data = hl,
                         ggplot2::aes(label = gene_symbol),
                         size = 2.2, hjust = -0.2, vjust = -0.5,
                         color = "grey10", fontface = "bold") +
    ggplot2::scale_color_manual(values = pal,
                                  name = "Translation class") +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::labs(
      title = "Candidate landscape: tissue vs serum",
      subtitle = sprintf(
        "n=%d genes; LTBP1/CDH13/SPARC = inverse stromal exemplars",
        nrow(d)),
      x = "Tissue effect (signed max |beta|)",
      y = "Serum log2FC (PDAC vs HC)") +
    pdactrace_axes_theme() +
    ggplot2::theme(legend.position = "bottom")
}
