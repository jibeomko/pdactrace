#' Render a pdactrace evidence report for one gene or a panel
#'
#' Renders a self-contained HTML report summarising the bundled
#' atlas's evidence for a single gene OR for a multi-gene panel. The
#' single-gene layout (length-1 input) carries the audit class + score,
#' the audit-component table, the six-axis evidence radar, the
#' cohort-adjusted stage trajectory, the per-cohort breakdown, and the
#' tissue-to-serum filter trace. The panel layout (length >= 2)
#' starts with a side-by-side comparison table from
#' [compare_candidates()] (audit class, RNA pattern, translation
#' class, cell origin, serum detectability, redundancy hint) and
#' follows with one condensed section per gene: explain_score()
#' rationale plus the evidence radar.
#'
#' The render is deterministic given (gene_symbol, atlas) because the
#' audit rule is frozen.
#'
#' @param gene_symbol HGNC gene symbol — character vector of length
#'   1 or more. Length 1 produces the single-gene template; length
#'   2+ produces the panel template.
#' @param output_file Optional output filename (without directory).
#'   Defaults to `paste0(<stem>, "_pdactrace_report.html")` where
#'   `<stem>` is the gene name for length-1 input or
#'   `paste0("panel_", length(gene_symbol), "gene")` otherwise.
#' @param output_dir Output directory. Defaults to `tempdir()`.
#' @param quiet Logical. If `TRUE` (default), suppress rmarkdown's
#'   chatty progress output.
#' @return Invisibly returns the absolute path to the rendered HTML
#'   file.
#' @examples
#' \donttest{
#'   if (requireNamespace("rmarkdown", quietly = TRUE) &&
#'       requireNamespace("knitr",     quietly = TRUE)) {
#'     # Single-gene report:
#'     fp1 <- report_gene("LTBP1", output_dir = tempdir())
#'
#'     # Panel report (3 genes):
#'     fp2 <- report_gene(c("LGALS3BP", "LTBP1", "GAPDH"),
#'                          output_dir = tempdir())
#'     file.exists(c(fp1, fp2))
#'   }
#' }
#' @export
report_gene <- function(gene_symbol,
                          output_file = NULL,
                          output_dir  = tempdir(),
                          quiet       = TRUE) {
  if (!is.character(gene_symbol) || length(gene_symbol) < 1L ||
      any(!nzchar(gene_symbol))) {
    stop("`gene_symbol` must be a non-empty character vector ",
         "(length >= 1).", call. = FALSE)
  }
  for (pkg in c("rmarkdown", "knitr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("`%s` is required for report_gene(). ", pkg),
           "Install via install.packages(\"", pkg, "\").",
           call. = FALSE)
    }
  }
  is_panel <- length(gene_symbol) > 1L
  template_name <- if (is_panel) "panel_report.Rmd" else "gene_report.Rmd"
  template <- system.file("rmd", template_name, package = "pdactrace")
  if (!nzchar(template)) {
    stop(sprintf(
      "Could not find inst/rmd/%s inside the installed pdactrace ",
      template_name),
      "package. Reinstall the package or check that inst/rmd/ is ",
      "present.", call. = FALSE)
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (is.null(output_file)) {
    output_file <- if (is_panel) {
      sprintf("panel_%dgene_pdactrace_report.html", length(gene_symbol))
    } else {
      paste0(gene_symbol, "_pdactrace_report.html")
    }
  }
  params <- if (is_panel) list(genes = gene_symbol)
            else list(gene = gene_symbol)
  out <- rmarkdown::render(
    input         = template,
    output_file   = output_file,
    output_dir    = output_dir,
    params        = params,
    envir         = new.env(parent = globalenv()),
    quiet         = quiet)
  invisible(normalizePath(out))
}
