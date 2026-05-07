#' Random-effects meta-analysis forest plot per gene
#'
#' Visualises the per-cohort log2FoldChange estimates with 95%
#' confidence intervals, plus the random-effects meta-summary diamond
#' for each of three contrasts: `Mid_vs_Early`, `Late_vs_Early`,
#' `Normal_vs_Early`. Heterogeneity (I2) is annotated.
#'
#' Reveals cohort divergence that the headline pattern label can hide.
#' E.g., for LTBP1 the forest plot shows TCGA mostly null, CPTAC
#' clearly elevated, and the meta-summary diamond crossing zero.
#'
#' @param gene_symbol HGNC standard symbol.
#' @param contrast Character, one of `"Mid_vs_Early"`,
#'   `"Late_vs_Early"`, `"Normal_vs_Early"`. Default `"Mid_vs_Early"`.
#' @param per_cohort_path Optional path to the per-cohort betas CSV
#'   (long-format with `gene_symbol`, `contrast`, `cohort`, `beta`,
#'   `se`, `n_ref`, `n_target`). Defaults to
#'   `${PDAC_BASE_DIR}/analysis/transcriptomics/results/figure1/per_cohort_betas_long.csv`,
#'   where `PDAC_BASE_DIR` is read from the environment. Returns
#'   `invisible(NULL)` with a message if the file is unavailable.
#' @return A `ggplot2` object.
#' @examples
#'   plot_meta_forest("LTBP1")
#'   plot_meta_forest("SERPINA1", contrast = "Normal_vs_Early")
#' @export
plot_meta_forest <- function(gene_symbol,
                              contrast = c("Mid_vs_Early",
                                            "Late_vs_Early",
                                            "Normal_vs_Early"),
                              per_cohort_path = NULL) {
  contrast <- match.arg(contrast)
  if (is.null(per_cohort_path)) {
    base_dir <- Sys.getenv("PDAC_BASE_DIR", unset = "")
    if (!nzchar(base_dir)) {
      message("plot_meta_forest(): set PDAC_BASE_DIR env var or pass ",
              "per_cohort_path; per-cohort betas not located.")
      return(invisible(NULL))
    }
    per_cohort_path <- file.path(
      base_dir, "analysis/transcriptomics",
      "results/figure1/per_cohort_betas_long.csv")
  }
  if (!file.exists(per_cohort_path)) {
    message("Per-cohort betas file not found at: ", per_cohort_path)
    return(invisible(NULL))
  }

  b <- data.table::fread(per_cohort_path)
  # rename function args to break NSE conflict with data.table columns
  .target_gene     <- gene_symbol
  .target_contrast <- contrast
  bsub <- b[b$gene_symbol == .target_gene & b$contrast == .target_contrast]
  if (nrow(bsub) == 0L) {
    message(sprintf("No %s data for %s.", contrast, gene_symbol))
    return(invisible(NULL))
  }

  ref <- .get_reference()
  meta_row <- ref[which(ref$gene_symbol == .target_gene)]
  if (nrow(meta_row) == 0L) {
    message("Gene not in atlas.")
    return(invisible(NULL))
  }
  meta_prefix <- switch(contrast,
    "Mid_vs_Early"   = "meta_MvE",
    "Late_vs_Early"  = "meta_LvE",
    "Normal_vs_Early" = "meta_NvE")
  meta_beta <- meta_row[[paste0(meta_prefix, "_beta")]]
  meta_se   <- meta_row[[paste0(meta_prefix, "_se")]]
  meta_I2   <- meta_row[[paste0(meta_prefix, "_I2")]]
  meta_k    <- meta_row[[paste0(meta_prefix, "_k")]]
  meta_pval <- meta_row[[paste0(meta_prefix, "_pval")]]
  meta_ci_lo <- meta_beta - 1.96 * meta_se
  meta_ci_hi <- meta_beta + 1.96 * meta_se

  bsub[, ci_lo := beta - 1.96 * se]
  bsub[, ci_hi := beta + 1.96 * se]
  bsub[, label := sprintf("%s (n=%d/%d)", cohort, n_ref, n_target)]

  meta_label <- sprintf("Meta (k=%d, I2=%.0f%%, p=%.3g)",
                          meta_k, meta_I2, meta_pval)

  rows <- data.table::data.table(
    label = c(bsub$label, meta_label),
    beta  = c(bsub$beta, meta_beta),
    ci_lo = c(bsub$ci_lo, meta_ci_lo),
    ci_hi = c(bsub$ci_hi, meta_ci_hi),
    is_meta = c(rep(FALSE, nrow(bsub)), TRUE))
  rows[, label := factor(label, levels = rev(rows$label))]

  ggplot2::ggplot(rows, ggplot2::aes(x = beta, y = label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                          color = "grey60", linewidth = 0.3) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = ci_lo, xmax = ci_hi, color = is_meta),
      height = 0.18, linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(color = is_meta, shape = is_meta),
                          size = 2.5, stroke = 0.6) +
    ggplot2::scale_color_manual(values = c("FALSE" = "#37474F",
                                              "TRUE" = "#C62828"),
                                   guide = "none") +
    ggplot2::scale_shape_manual(values = c("FALSE" = 19, "TRUE" = 23),
                                   guide = "none") +
    ggplot2::labs(
      title = sprintf("%s - random-effects meta forest (%s)",
                          gene_symbol, contrast),
      subtitle = sprintf(
        "diamond = REML meta-estimate (beta=%.2f, 95%% CI [%.2f, %.2f])",
        meta_beta, meta_ci_lo, meta_ci_hi),
      x = "log2FC (per-cohort) / meta beta",
      y = NULL) +
    pdactrace_axes_theme()
}
