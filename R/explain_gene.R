#' Print a per-axis evidence summary for one gene
#'
#' User-facing text formatter that wraps the two evidence layers
#' already in the package: [format_provenance()] (plain-English
#' phase-tag rollup) and [evidence_math()] (per-axis mathematical
#' values). Mirrors [explain_score()] in pattern: prints to console
#' when `verbose = TRUE`, returns the underlying structured data
#' invisibly either way.
#'
#' Three views are supported:
#'
#' \describe{
#'   \item{`"evidence"`}{Plain-English provenance summary only —
#'     same content [format_provenance()] produces, plus the
#'     gene's `audit_class` headline.}
#'   \item{`"math"`}{Per-axis Evidence Math output — eight axes
#'     (trajectory_fit, effect_magnitude, cohort_consistency,
#'     rna_protein_coupling, serum_bridge, cell_specificity,
#'     filter_survival, clinical_role).}
#'   \item{`"both"`}{Both, separated by a horizontal rule.}
#' }
#'
#' Per-axis printing is the load-bearing design choice: it keeps
#' the output interpretable rather than hiding the evidence behind
#' a single composite score.
#'
#' @param gene_symbol HGNC gene symbol (length-1 character).
#' @param view One of `"evidence"` (default), `"math"`, or `"both"`.
#' @param verbose Logical. If `TRUE` (default), prints the formatted
#'   sections to the console. The structured list is returned
#'   invisibly either way.
#' @param reference Optional `data.table` to inject in place of the
#'   bundled atlas (for tests / downstream pipelines).
#' @return Invisibly, a list with elements `gene`, `view`,
#'   `provenance` (length-1 character or NA), and `math` (list-of-
#'   axes from [evidence_math()] or NULL).
#' @examples
#' explain_gene("LGALS3BP", view = "math", verbose = FALSE)
#' explain_gene("LTBP1", view = "evidence", verbose = FALSE)
#' explain_gene("LTBP1", view = "both", verbose = FALSE)
#' @seealso [evidence_math()], [format_provenance()],
#'   [explain_score()].
#' @export
explain_gene <- function(gene_symbol,
                          view = c("evidence", "math", "both"),
                          verbose = TRUE,
                          reference = NULL) {
  view <- match.arg(view)
  if (!is.character(gene_symbol) || length(gene_symbol) != 1L) {
    stop("`gene_symbol` must be a length-1 character string.",
         call. = FALSE)
  }
  ref <- .get_reference(reference)
  target <- gene_symbol
  row <- ref[gene_symbol == target]
  if (nrow(row) == 0L) {
    stop(sprintf("Gene '%s' is not in the bundled atlas.",
                 gene_symbol), call. = FALSE)
  }

  prov <- as.character(row$provenance)
  audit_class <- as.character(row$audit_class)

  prov_text <- if (length(prov) == 0L || is.na(prov))
    "(no provenance recorded)" else format_provenance(prov, "verbose")
  prov_compact <- if (length(prov) == 0L || is.na(prov))
    "(none)" else format_provenance(prov, "compact")

  math <- if (view %in% c("math", "both")) {
    evidence_math(gene_symbol, reference = reference)
  } else NULL

  if (isTRUE(verbose)) {
    .ex_header(gene_symbol, audit_class)
    if (view %in% c("evidence", "both")) {
      cat("\nEvidence (plain-English provenance)\n")
      cat("  ", prov_compact, "\n", sep = "")
      cat(prov_text, "\n", sep = "")
    }
    if (view == "both") {
      cat("\n", strrep("-", 60), "\n", sep = "")
    }
    if (view %in% c("math", "both")) {
      .ex_print_math(math)
    }
  }

  invisible(list(gene = gene_symbol,
                 view = view,
                 audit_class = audit_class,
                 provenance = prov_compact,
                 math = math))
}

# ---- formatting helpers (internal) -------------------------------------

.ex_header <- function(gene, audit_class) {
  cat(gene, " — Evidence Math Summary",
      if (!is.na(audit_class)) sprintf("  [%s]", audit_class) else "",
      "\n", sep = "")
}

.ex_print_math <- function(m) {
  .ex_section_trajectory(m$trajectory_fit)
  .ex_section_effect(m$effect_magnitude)
  .ex_section_cohort(m$cohort_consistency)
  .ex_section_coupling(m$rna_protein_coupling)
  .ex_section_serum(m$serum_bridge)
  .ex_section_cell(m$cell_specificity)
  .ex_section_filter(m$filter_survival)
  .ex_section_clinical(m$clinical_role)
}

.fmt_n <- function(x, digits = 3) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) return("NA")
  if (is.integer(x)) return(as.character(x))
  if (is.numeric(x)) sprintf(paste0("%.", digits, "f"), x)
  else as.character(x)
}

.fmt_pct <- function(x) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) return("NA")
  sprintf("%.0f%%", as.numeric(x) * 100)
}

.fmt_signed <- function(x, digits = 2) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) return("NA")
  sprintf(paste0("%+.", digits, "f"), x)
}

.ex_section_trajectory <- function(t) {
  cat("\nTrajectory fit\n")
  cat("  pattern        : ", .fmt_n(t$rna_pattern), "\n", sep = "")
  cat("  rho_best       : ", .fmt_n(t$rho_best, 3), "\n", sep = "")
  cat("  rho_runner_up  : ", .fmt_n(t$rho_runner_up, 3), "\n", sep = "")
  cat("  delta_rho      : ", .fmt_n(t$delta_rho, 3),
      "    (specificity margin; >0.10 = clean)\n", sep = "")
  if (!is.na(t$note))
    cat("  note           : ", t$note, "\n", sep = "")
}

.ex_section_effect <- function(e) {
  cat("\nEffect magnitude\n")
  cat("  ||beta_RNA||_2 : ", .fmt_n(e$rna_beta_norm),
      "    (target stage: ", .fmt_n(e$rna_target_stage), ")\n", sep = "")
  cat("  max|beta_RNA|  : ", .fmt_n(e$rna_beta_max_abs),
      "  at stage ", .fmt_n(e$rna_max_at_stage), "\n", sep = "")
  if (is.na(e$prot_beta_norm)) {
    cat("  protein side   : not in pdactrace_protein_betas\n")
  } else {
    cat("  ||beta_prot||_2: ", .fmt_n(e$prot_beta_norm), "\n", sep = "")
    cat("  max|beta_prot| : ", .fmt_n(e$prot_beta_max_abs),
        "  at stage ", .fmt_n(e$prot_max_at_stage), "\n", sep = "")
  }
}

.ex_section_cohort <- function(c) {
  cat("\nCohort consistency\n")
  cat("  Stouffer Z     : ", .fmt_n(c$stouffer_z),
      "    (padj = ", .fmt_n(c$stouffer_padj, 4), ")\n", sep = "")
  cat("  agreement      : ", .fmt_n(c$cohort_agreement, 2),
      "  fraction of cohorts agreeing\n", sep = "")
  cat("  max meta I2    : ", .fmt_pct(c$max_meta_I2 / 100), "\n", sep = "")
}

.ex_section_coupling <- function(c) {
  cat("\nRNA-protein coupling\n")
  if (!isTRUE(c$prot_in_atlas)) {
    cat("  protein side   : not in pdactrace_protein_betas\n")
    cat("  prot_tier      : ", .fmt_n(c$prot_tier), "\n", sep = "")
    return(invisible(NULL))
  }
  cat("  cosine(beta_RNA, beta_prot) : ", .fmt_n(c$cosine), "\n", sep = "")
  cat("  prot_pattern   : ", .fmt_n(c$prot_pattern), "\n", sep = "")
  cat("  prot_tier      : ", .fmt_n(c$prot_tier), "\n", sep = "")
  cat("  concordant     : ", .fmt_n(c$rnaprot_concordant), "\n", sep = "")
  if (!is.na(c$note))
    cat("  note           : ", c$note, "\n", sep = "")
}

.ex_section_serum <- function(s) {
  cat("\nSerum bridge\n")
  cat("  translation_class      : ", .fmt_n(s$translation_class),
      "\n", sep = "")
  cat("  serum log2FC PDAC vs HC: ", .fmt_signed(s$serum_log2fc_PDAC_vs_HC),
      "\n", sep = "")
  cat("  serum log2FC Pan  vs HC: ", .fmt_signed(s$serum_log2fc_Pan_vs_HC),
      "\n", sep = "")
  cat("  detected in cohorts    : ", .fmt_n(s$serum_n_cohorts_detected),
      "\n", sep = "")
  cat("  phase77 strict         : ", .fmt_n(s$phase77_strict),
      "\n", sep = "")
}

.ex_section_cell <- function(c) {
  cat("\nCell specificity\n")
  cat("  cell_origin_top: ", .fmt_n(c$cell_origin_top), "\n", sep = "")
  cat("  tau index      : ", .fmt_n(c$tau),
      "    (Yanai et al.; ~1 = concentrated in 1 cell type)\n", sep = "")
  cat("  origin padj    : ", .fmt_n(c$cell_origin_padj, 4), "\n", sep = "")
}

.ex_section_filter <- function(f) {
  cat("\nFilter survival\n")
  cat("  passed ", f$passed, " / ", f$total, " steps\n", sep = "")
  cat("  per_step       : ", paste(names(f$per_step),
                                     ifelse(is.na(f$per_step), "NA",
                                            ifelse(f$per_step, "T", "F")),
                                     sep = "=", collapse = ", "),
      "\n", sep = "")
}

.ex_section_clinical <- function(c) {
  cat("\nClinical role\n")
  cat("  resectable_marker : ", .fmt_n(c$resectable_marker), "\n", sep = "")
  cat("  panel_member      : ", .fmt_n(c$panel_member), "\n", sep = "")
}
