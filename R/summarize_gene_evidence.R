#' Human-readable text summary of gene evidence
#'
#' Returns a multi-line, glue-formatted text summary of all evidence
#' available for `gene_symbol` in the pdactrace reference atlas.
#'
#' @param gene_symbol HGNC standard symbol.
#' @param detail Logical. If `TRUE`, append per-stage, per-cohort,
#'   pancreatitis raw-mean, and filter-metric detail blocks.
#' @return A character scalar (single string with `\\n` separators).
#'   Returns `NULL` invisibly with a message if gene is not in atlas.
#' @examples
#'   cat(summarize_gene_evidence("LTBP1"))
#' @export
summarize_gene_evidence <- function(gene_symbol, detail = FALSE) {
  q <- query_gene(gene_symbol)
  if (is.null(q)) return(invisible(NULL))

  # -- RNA layer ----------------------------------------------
  if (is.na(q$rna$pattern)) {
    rna_line <- if (isTRUE(q$rna$excluded_mid_pattern))
      glue::glue(
        "RNA: stage-progressive but excluded by design ",
        "(Mid-onset pattern; LRT padj = {fmt_p(q$rna$lrt_padj)}).")
    else
      glue::glue(
        "RNA: stage-progressive but unclassified at rho >= 0.85 ",
        "(LRT padj = {fmt_p(q$rna$lrt_padj)}; ",
        "best-template rho = {round(q$rna$pattern_rho, 2)}).")
  } else {
    rna_line <- glue::glue(
      "RNA: {q$rna$pattern} (rho = {round(q$rna$pattern_rho, 2)}, ",
      "LRT padj = {fmt_p(q$rna$lrt_padj)}, ",
      "Stouffer Z = {round(q$rna$stouffer_z, 2)}, ",
      "{round(q$rna$cohort_agreement * 100)}% cohort agreement).")
  }

  # -- Protein layer ------------------------------------------
  prot_line <- if (is.na(q$protein$pattern))
    "Tissue protein: no concordant Early-onset pattern."
  else
    glue::glue(
      "Tissue protein: {q$protein$pattern} ",
      "[{q$protein$tier %||% 'untiered'}; ",
      "concordant = {q$protein$concordant}].")

  # -- scRNA layer --------------------------------------------
  scrna_line <- if (is.na(q$scrna$top_celltype)) {
    "scRNA: not in cell-origin reference."
  } else {
    distrib <- q$scrna$distribution[[1]]
    top3 <- if (length(distrib) >= 3) {
      paste(sprintf("%s (%.0f%%)",
                     names(distrib)[1:3],
                     100 * distrib[1:3] / sum(distrib)),
            collapse = ", ")
    } else {
      q$scrna$top_celltype
    }
    glue::glue("Cell origin: {top3}.")
  }

  # -- Serum layer --------------------------------------------
  serum_line <- if (!isTRUE(q$serum$detected)) {
    "Serum: not detected in current MS cohorts."
  } else {
    cls <- q$serum$translation_class
    cls_text <- switch(as.character(cls),
                        "A" = "Class A (concordant)",
                        "B" = "Class B (inverse stromal, manually curated)",
                        "C" = "Class C (decoupled)",
                        "unclassified")
    lfc_pdac <- q$serum$log2fc_PDAC_vs_HC
    lfc_str <- if (is.na(lfc_pdac)) "log2FC NA" else
      sprintf("log2FC PDAC vs HC = %.2f", lfc_pdac)
    p77_str <- if (isTRUE(q$serum$phase77_strict))
      ", phase77 strict candidate" else ""
    plural <- if (q$serum$n_cohorts_detected > 1) "s" else ""
    glue::glue(
      "Serum: detected (n={q$serum$n_cohorts_detected} cohort{plural}), ",
      "{cls_text}, {lfc_str}{p77_str}.")
  }

  # -- Clinical layer -----------------------------------------
  clin_parts <- character()
  if (isTRUE(q$clinical$resectable_marker)) {
    p29 <- q$clinical$resectable_pattern_phase29
    clin_parts <- c(clin_parts,
      glue::glue("phase29 resectable marker (s_pat = {p29 %||% 'NA'})"))
  }
  if (isTRUE(q$clinical$panel_member)) {
    clin_parts <- c(clin_parts, "panel member (LTBP1+SERPINA1 / hybrid)")
  }
  clin_line <- if (length(clin_parts) == 0L)
    "Clinical: no resectable / panel flag."
  else
    glue::glue("Clinical: ", paste(clin_parts, collapse = "; "), ".")

  # -- Filter trace (phase60 7-step) --------------------------
  filt <- q$filter_status
  filt_line <- if (is.null(filt) || all(is.na(unlist(filt)))) {
    "Filter trace: gene not evaluated by phase60 7-step pipeline."
  } else {
    flt_pass <- c(
      if (isTRUE(filt$signal_peptide))    "SignalP" else NULL,
      if (isTRUE(filt$serum_measurable))  "serum-measurable" else NULL,
      if (isTRUE(filt$serum_significant)) "serum-significant" else NULL,
      if (isTRUE(filt$pancreatitis_pdac)) "pan-vs-PDAC" else NULL,
      if (isTRUE(filt$pancreatitis_hc))   "HC-in-middle" else NULL,
      if (isTRUE(filt$direction_match))   "direction-match" else NULL,
      if (isTRUE(filt$final_pass))        "FINAL" else NULL)
    flt_n_pass <- length(flt_pass)
    pass_str <- if (flt_n_pass == 0L) "none" else
      paste(flt_pass, collapse = ", ")
    glue::glue(
      "Filter trace (phase60 7-step): {flt_n_pass}/7 passed ",
      "[{pass_str}].")
  }

  prov_line <- paste0(
    "Evidence:  ", format_provenance(q$provenance, "compact"), "\n",
    "Technical: ", format_provenance(q$provenance, "raw"))

  # -- Detail block (optional) --------------------------------
  detail_lines <- character()
  if (isTRUE(detail)) {
    d <- query_gene_detailed(gene_symbol)
    if (!is.null(d)) {
      # Per-stage one-liner
      ps <- d$per_stage
      ps_summary <- paste(sprintf("%s=%+.2f%s",
                                       ps$stage, ps$log2FC,
                                       ifelse(ps$significant, "*", "")),
                            collapse = ", ")
      detail_lines <- c(detail_lines,
        glue::glue("Per-stage:    {ps_summary} (* Wald padj<0.05)"))

      # Per-cohort one-liner
      pc <- d$per_cohort
      if (nrow(pc) > 0L) {
        cv <- paste(sprintf("%s %s%s", pc$cohort,
                              substr(pc$trend, 1, 3),
                              ifelse(pc$monotonic, "up", ".")),
                      collapse = " | ")
        st_p <- attr(pc, "stouffer_p")
        agr  <- attr(pc, "agreement_pct")
        detail_lines <- c(detail_lines, glue::glue(
          "Per-cohort:   {cv} (Stouffer p={fmt_p(st_p)}, agree {round(100*agr)}%)"))
      }

      # Pancreatitis raw means + t-test
      if (!is.na(d$serum_per_cohort$mean[1])) {
        sp <- d$serum_per_cohort
        pp <- attr(sp, "pdac_vs_pan_pval")
        detail_lines <- c(detail_lines, glue::glue(
          "Pancreatitis: PDAC={round(sp$mean[1],2)}, ",
          "Pan={round(sp$mean[2],2)}, HC={round(sp$mean[3],2)} ",
          "(PDAC-vs-Pan t-pval={fmt_p(pp)})"))
      }

      # Filter underlying numbers (only show fails for cleanness)
      fd <- d$filter_diag
      fail_steps <- fd[!is.na(pass) & pass == FALSE]
      if (nrow(fail_steps) > 0L) {
        detail_lines <- c(detail_lines, glue::glue(
          "Filter fails: {paste(sprintf('%s (%s)', fail_steps$step, fail_steps$underlying_metric), collapse=' | ')}"))
      }
    }
  }

  paste(c(
    glue::glue("=== {gene_symbol} ==="),
    rna_line, prot_line, scrna_line, serum_line, clin_line, filt_line,
    detail_lines,
    prov_line),
    collapse = "\n")
}

# Tiny helper to format p-values
fmt_p <- function(p) {
  if (is.null(p) || is.na(p)) return("NA")
  if (p < 1e-300) return("< 1e-300")
  if (p < 1e-3)   return(formatC(p, format = "e", digits = 1))
  formatC(p, format = "g", digits = 3)
}

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a
