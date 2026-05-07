#' Render a one-gene pdactrace evidence report
#'
#' Renders a self-contained HTML report summarising the bundled
#' atlas's evidence for a single gene: the audit class + score, key
#' panel fields, the six-axis evidence radar, the cohort-adjusted
#' stage trajectory, the per-cohort breakdown, and the
#' tissue-to-serum filter trace. The report uses the package's
#' NCS-grade ggplot theme throughout and is intended for sharing a
#' one-page evidence dossier with a collaborator.
#'
#' This is a thin wrapper around [rmarkdown::render()] on the bundled
#' template `inst/rmd/gene_report.Rmd`. The render is deterministic
#' given (gene_symbol, atlas) because the audit rule is frozen.
#'
#' @param gene_symbol HGNC gene symbol (length-1 character).
#' @param output_file Optional output filename (without directory).
#'   Defaults to `paste0(gene_symbol, "_pdactrace_report.html")`.
#' @param output_dir Output directory. Defaults to `tempdir()`.
#' @param quiet Logical. If `TRUE` (default), suppress rmarkdown's
#'   chatty progress output.
#' @return Invisibly returns the absolute path to the rendered HTML
#'   file.
#' @examples
#' \donttest{
#'   if (requireNamespace("rmarkdown", quietly = TRUE) &&
#'       requireNamespace("knitr",     quietly = TRUE)) {
#'     fp <- report_gene("LTBP1", output_dir = tempdir())
#'     file.exists(fp)
#'   }
#' }
#' @export
report_gene <- function(gene_symbol,
                          output_file = NULL,
                          output_dir  = tempdir(),
                          quiet       = TRUE) {
  if (!is.character(gene_symbol) || length(gene_symbol) != 1L ||
      !nzchar(gene_symbol)) {
    stop("`gene_symbol` must be a length-1, non-empty character string.")
  }
  for (pkg in c("rmarkdown", "knitr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("`%s` is required for report_gene(). ", pkg),
           "Install it via install.packages(\"", pkg, "\").",
           call. = FALSE)
    }
  }
  template <- system.file("rmd", "gene_report.Rmd", package = "pdactrace")
  if (!nzchar(template)) {
    stop("Could not find inst/rmd/gene_report.Rmd inside the installed ",
         "pdactrace package. Reinstall the package or check that ",
         "inst/rmd/ is present.", call. = FALSE)
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (is.null(output_file)) {
    output_file <- paste0(gene_symbol, "_pdactrace_report.html")
  }
  out <- rmarkdown::render(
    input         = template,
    output_file   = output_file,
    output_dir    = output_dir,
    params        = list(gene = gene_symbol),
    envir         = new.env(parent = globalenv()),
    quiet         = quiet)
  invisible(normalizePath(out))
}
