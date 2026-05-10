#' Per-gene serum direction context, coloured by translation class
#'
#' Bar plot of per-gene serum log2FC against healthy control (HC)
#' and against pancreatitis (Pan), faceted by contrast and
#' coloured by translation class (Class A / B / C). These two
#' values together drive the `direction-match` step of the phase60
#' tissue-to-serum filter audit; the panel is the visual analog of
#' the `serum_log2fc_PDAC_vs_HC` and `serum_log2fc_Pan_vs_HC`
#' columns from [pdactrace_reference], with the categorical
#' `translation_class` mapped to the bar fill.
#'
#' Used as a standalone diagnostic and as one of the six
#' [viz_gene()] panels.
#'
#' @param gene Character vector of HGNC gene symbols.
#' @return A `ggplot2` object, or `NULL` invisibly if none of the
#'   requested genes are in the bundled atlas.
#' @examples
#'   plot_serum_direction("LTBP1")
#'   plot_serum_direction(c("SERPINA1", "LTBP1", "LGALS3BP"))
#' @seealso [plot_filter_trace()] which embeds the same panel
#'   above its filter grid by default; [viz_gene()] for the
#'   6-panel composite.
#' @export
plot_serum_direction <- function(gene) {
  if (!is.character(gene) || length(gene) == 0L) {
    stop("`gene` must be a non-empty character vector.",
         call. = FALSE)
  }
  .ft_serum_strip(unique(as.character(gene)))
}
