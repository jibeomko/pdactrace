#' Human-readable evidence labels for the bundled phase tags
#'
#' Maps the internal `phase33`, `phase34`, `stouffer_consistency`, ...
#' provenance tags carried in `pdactrace_reference$provenance` onto
#' plain-English evidence-source labels. The phase tags remain the
#' canonical IDs for reproducibility (they correspond to the numbered
#' analysis scripts in the companion manuscript-monorepo) but the
#' user-facing print path now leads with the labels and folds the raw
#' tag list into a `Technical:` footer.
#'
#' Intended for internal use by [query_gene()],
#' [query_gene_detailed()], and [summarize_gene_evidence()] —
#' exported so end users who consume the raw `provenance` column can
#' relabel it the same way.
#'
#' @param provenance Character of length 1 or longer. Either a
#'   comma-separated phase string (e.g. `"phase33,phase34,phase60"`)
#'   or a character vector of individual tags. Unknown tags fall
#'   through unchanged with a `(unmapped)` suffix in `verbose` mode.
#' @param style One of:
#'   * `"compact"` (default) — returns a one-line plus-separated
#'     summary suitable for inline print, e.g.
#'     `"RNA + Tissue protein + Multi-cohort RNA + ..."`.
#'   * `"verbose"` — returns a multi-line hyphen-prefixed list with
#'     "label: explanation" per evidence source.
#'   * `"raw"` — returns the input phase tags themselves (useful for
#'     the technical-provenance footer).
#' @return Character of length 1.
#' @examples
#' format_provenance("phase33,phase34,phase60", style = "compact")
#' format_provenance("phase33,phase34,phase60", style = "verbose")
#' format_provenance("phase33,phase34,phase60", style = "raw")
#' @export
format_provenance <- function(provenance,
                                style = c("compact", "verbose", "raw")) {
  style <- match.arg(style)
  if (is.null(provenance) || length(provenance) == 0L ||
      all(is.na(provenance))) {
    return("(no provenance recorded)")
  }
  tags <- if (length(provenance) == 1L) {
    strsplit(as.character(provenance), "[,;[:space:]]+")[[1L]]
  } else {
    as.character(provenance)
  }
  tags <- tags[nzchar(tags)]
  if (length(tags) == 0L) return("(no provenance recorded)")

  if (style == "raw") {
    return(paste(tags, collapse = ", "))
  }

  labels <- vapply(tags, function(t) {
    .PROVENANCE_LABEL_MAP[[t]] %||% paste0(t, " (unmapped)")
  }, character(1L))
  details <- vapply(tags, function(t) {
    .PROVENANCE_DETAIL_MAP[[t]] %||% NA_character_
  }, character(1L))

  if (style == "compact") {
    return(paste(labels, collapse = " + "))
  }

  # verbose
  out <- character(0L)
  for (i in seq_along(tags)) {
    if (is.na(details[i])) {
      out <- c(out, paste0("- ", labels[i]))
    } else {
      out <- c(out, paste0("- ", labels[i], ": ", details[i]))
    }
  }
  paste(out, collapse = "\n")
}

# Internal lookup tables ----------------------------------------------
.PROVENANCE_LABEL_MAP <- c(
  phase33              = "RNA trajectory",
  phase34              = "Tissue protein trajectory",
  stouffer_consistency = "Multi-cohort RNA consistency",
  phase2c              = "scRNA cell origin",
  phase42              = "Serum / pancreatitis comparison",
  phase60              = "7-step serum filter audit",
  phase77              = "Strict RNA-protein-serum bridge",
  phase29              = "Resectable-stage marker screen",
  phase80              = "Predeclared panel member")

.PROVENANCE_DETAIL_MAP <- c(
  phase33 = "matched in bulk RNA-seq stage model (12-template)",
  phase34 = "matched in pooled tissue proteomics (12-template)",
  stouffer_consistency =
    "supported across RNA cohorts by Stouffer meta-analysis",
  phase2c = "cell-type specificity available from scRNA atlas",
  phase42 = "evaluated for serum abundance / pancreatitis context",
  phase60 = "evaluated by the 7-step tissue-to-serum filter audit",
  phase77 = paste0(
    "strict RNA-protein-serum detected/concordance class available"),
  phase29 = "included in resectable-stage marker screen",
  phase80 = "member of a predeclared serum panel evaluation")

# null-coalesce
`%||%` <- function(a, b) if (is.null(a)) b else a
