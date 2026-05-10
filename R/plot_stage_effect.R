#' Forest plot of per-stage effect size
#'
#' For RNA, draws log2FC +/- 1.96 x lfcSE per stage with
#' significance asterisks at per-stage Wald padj < 0.05. For
#' tissue protein, draws log2FC point estimates per stage; the
#' bundled per-stage SE is unavailable so no CI is drawn (the
#' overall F-test padj is reported in the subtitle instead).
#'
#' @param gene_symbol HGNC symbol.
#' @param layer One of `"rna"` (default) or `"protein"`. The
#'   protein layer reads from [pdactrace_protein_betas].
#' @return A `ggplot2` object.
#' @examples
#'   plot_stage_effect("LTBP1")
#'   plot_stage_effect("LTBP1", layer = "protein")
#' @export
plot_stage_effect <- function(gene_symbol,
                               layer = c("rna", "protein")) {
  layer <- match.arg(layer)
  if (layer == "rna") return(.plot_stage_effect_rna(gene_symbol))
  .plot_stage_effect_protein(gene_symbol)
}

.plot_stage_effect_rna <- function(gene_symbol) {
  d <- query_gene_detailed(gene_symbol)
  if (is.null(d)) return(invisible(NULL))
  ps <- d$per_stage

  ps$sig_label <- ifelse(ps$significant, "*", "")

  ggplot2::ggplot(ps, ggplot2::aes(x = stage, y = log2FC)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60",
                          linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
      width = 0.18, linewidth = 0.45, color = "#37474F") +
    ggplot2::geom_point(size = 2.5, color = "#0D47A1") +
    ggplot2::geom_text(ggplot2::aes(y = ci_hi + 0.05,
                                        label = sig_label),
                         color = "#C62828", size = 4,
                         fontface = "bold", vjust = 0) +
    ggplot2::labs(
      title = sprintf("%s - RNA per-stage log2FC (mean +/- 95%% CI)",
                          gene_symbol),
      subtitle = "* Wald padj < 0.05 (per-contrast vs Normal, cohort-adjusted)",
      x = NULL, y = "log2FC vs Normal") +
    pdactrace_axes_theme()
}

.plot_stage_effect_protein <- function(gene_symbol) {
  e <- new.env()
  utils::data("pdactrace_protein_betas", package = "pdactrace",
              envir = e)
  pb <- data.table::as.data.table(e$pdactrace_protein_betas)
  target <- gene_symbol
  row <- pb[gene_symbol == target]
  if (nrow(row) == 0L) {
    return(.plot_stage_effect_empty(
      sprintf("%s - tissue protein (not in atlas)", gene_symbol)))
  }
  ps <- data.table::data.table(
    stage  = factor(c("Normal", "Early", "Mid", "Late"),
                    levels = c("Normal", "Early", "Mid", "Late")),
    log2FC = c(as.numeric(row$prot_beta_N),
               as.numeric(row$prot_beta_E),
               as.numeric(row$prot_beta_M),
               as.numeric(row$prot_beta_L)))
  padj <- as.numeric(row$prot_lrt_padj)
  sub  <- if (is.na(padj))
    "tissue protein 12-template fit (no per-stage SE; F-test padj NA)"
    else
    sprintf("F-test padj = %s   (no per-stage SE: limma F-test; per-contrast Wald deferred to v0.99.8+)",
            formatC(padj, format = "g", digits = 3))
  ggplot2::ggplot(ps, ggplot2::aes(x = stage, y = log2FC)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60",
                          linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_line(ggplot2::aes(group = 1),
                        color = "#6A1B9A", linewidth = 0.4,
                        linetype = "dotted") +
    ggplot2::geom_point(size = 2.5, color = "#6A1B9A") +
    ggplot2::labs(
      title = sprintf("%s - tissue protein per-stage log2FC",
                          gene_symbol),
      subtitle = sub,
      x = NULL, y = "log2FC vs Normal") +
    pdactrace_axes_theme()
}

.plot_stage_effect_empty <- function(title) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 1, y = 1,
                        label = "no data",
                        color = "grey40", size = 3) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    pdactrace_axes_theme() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                    axis.ticks = ggplot2::element_blank())
}
