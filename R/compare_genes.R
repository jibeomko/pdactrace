#' Side-by-side Evidence Math comparison across multiple genes
#'
#' Calls [evidence_math()] for each input gene and pivots the
#' per-axis lists into a tidy `data.table`. Useful for assembling
#' the "candidate vs. competitors" comparison tables that show up
#' in the manuscript and in `report_gene()` output.
#'
#' By default returns a long table with one row per
#' (gene x axis x metric) tuple — friendly to `data.table::dcast()`
#' if you need a custom wide layout. Set `wide = TRUE` for the
#' convenience pivot to one row per gene with metric columns
#' (`axis.metric`, e.g. `trajectory_fit.delta_rho`).
#'
#' Note: this function reports the *math evidence* per axis. For
#' the audit-score-driven candidate ranking with redundancy
#' grouping, see [compare_candidates()] (different layer, same
#' atlas).
#'
#' @param gene_symbols Character vector of HGNC gene symbols.
#' @param axes Optional character vector to filter the axis set
#'   reported. Defaults to all axes from [evidence_math()].
#'   Recognised values: `"trajectory_fit"`, `"effect_magnitude"`,
#'   `"cohort_consistency"`, `"rna_protein_coupling"`,
#'   `"serum_bridge"`, `"cell_specificity"`,
#'   `"filter_survival"`, `"clinical_role"`.
#' @param wide Logical. If `TRUE`, returns one row per gene with
#'   `axis.metric` columns. Default `FALSE` (long form).
#' @param reference Optional `data.table` to inject in place of
#'   the bundled atlas (used by tests).
#' @return A `data.table`. Long form columns: `gene`, `axis`,
#'   `metric`, `value` (character; numerics formatted to 4 sig
#'   figs to keep the table renderable). Wide form: `gene`
#'   followed by `axis.metric` columns.
#' @examples
#' compare_genes(c("LGALS3BP", "LTBP1", "TIMP1"))
#' compare_genes(c("LGALS3BP", "LTBP1"),
#'               axes = c("trajectory_fit", "rna_protein_coupling"),
#'               wide = TRUE)
#' @seealso [evidence_math()], [explain_gene()],
#'   [compare_candidates()].
#' @export
compare_genes <- function(gene_symbols,
                           axes = NULL,
                           wide = FALSE,
                           reference = NULL) {
  if (!is.character(gene_symbols) || length(gene_symbols) == 0L) {
    stop("`gene_symbols` must be a non-empty character vector.",
         call. = FALSE)
  }
  all_axes <- c("trajectory_fit", "effect_magnitude",
                "cohort_consistency", "rna_protein_coupling",
                "serum_bridge", "cell_specificity",
                "filter_survival", "clinical_role")
  if (is.null(axes)) {
    axes <- all_axes
  } else {
    bad <- setdiff(axes, all_axes)
    if (length(bad) > 0L) {
      stop("Unknown axis name(s): ", paste(bad, collapse = ", "),
           ".\n  Valid: ", paste(all_axes, collapse = ", "),
           call. = FALSE)
    }
  }

  rows <- list()
  for (g in gene_symbols) {
    m <- tryCatch(evidence_math(g, reference = reference),
                  error = function(e) NULL)
    if (is.null(m)) {
      rows[[g]] <- data.table::data.table(
        gene = g, axis = NA_character_,
        metric = "missing",
        value = "gene not in bundled atlas")
      next
    }
    for (a in axes) {
      ax_list <- m[[a]]
      if (is.null(ax_list)) next
      flat <- .cg_flatten_axis(ax_list)
      rows[[paste(g, a)]] <- data.table::data.table(
        gene = g, axis = a,
        metric = names(flat),
        value = flat)
    }
  }
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)

  if (isTRUE(wide)) {
    out <- out[!is.na(axis)]
    out[, key := paste(axis, metric, sep = ".")]
    wide_dt <- data.table::dcast(out, gene ~ key, value.var = "value")
    return(wide_dt[])
  }
  out[]
}

# ---- internal: flatten one axis list into named character vector ------

.cg_flatten_axis <- function(ax) {
  parts <- character(0L)
  for (nm in names(ax)) {
    v <- ax[[nm]]
    if (is.null(v)) {
      parts[[nm]] <- NA_character_
      next
    }
    if (length(v) > 1L) {
      # per-step booleans (filter_survival$per_step)
      for (sub in names(v)) {
        parts[[paste(nm, sub, sep = ".")]] <- .cg_format_scalar(v[[sub]])
      }
      next
    }
    parts[[nm]] <- .cg_format_scalar(v)
  }
  parts
}

.cg_format_scalar <- function(x) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) return(NA_character_)
  if (is.logical(x)) return(as.character(x))
  if (is.numeric(x)) return(formatC(x, format = "g", digits = 4))
  as.character(x)
}
