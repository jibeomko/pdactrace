#' Assemble user multi-omics inputs into a per-gene evidence table
#'
#' Combines optional per-layer user inputs into a single per-gene
#' evidence data.table whose column schema is the subset of the
#' bundled `pdactrace_reference` atlas that
#' `compute_audit_score(evidence = ...)` consumes. Layers the user
#' does not supply are filled with `NA`; the downstream feature
#' constructor in `.audit_compute_features()` is NA-tolerant and
#' degrades gracefully (a layer absent → that layer contributes
#' zero to `evidence_strength`; cross-layer concordance falls to
#' the available pair; etc.).
#'
#' @param rna_fit Optional output of [classify_trajectory()] (i.e., a
#'   data.table with `gene_symbol`, `rna_pattern`,
#'   `rna_pattern_rho`, `lrt_padj`).
#' @param prot_fit Optional output of [classify_protein_trajectory()]
#'   (a data.table with `gene_symbol`, `prot_pattern`,
#'   `prot_pattern_rho`).
#' @param scrna_summary Optional data.table with `gene_symbol` and
#'   `cell_origin_top` columns. Cell-of-origin label per gene.
#' @param serum_summary Optional data.table with `gene_symbol`,
#'   `serum_detected` (logical), and `serum_log2fc_PDAC_vs_HC`
#'   (numeric).
#' @param signal_peptide Optional character vector of gene symbols
#'   that carry a signal peptide (UniProt SP-positive). Genes outside
#'   this set are flagged `flt_signal_peptide = FALSE`.
#' @param cross_cohort_agreement Optional named numeric vector,
#'   names = gene symbols, values in `[0, 1]` representing fraction
#'   of cohorts agreeing on direction. When `NULL`, defaults to
#'   `1.0` for genes with `rna_pattern` set and `NA` otherwise.
#' @param max_I2 Optional named numeric vector, names = gene
#'   symbols, values are max meta-analysis I². When `NULL`, defaults
#'   to `0` (no heterogeneity penalty).
#' @return A `data.table` with the column schema consumed by
#'   `compute_audit_score(evidence = ...)`.
#' @examples
#' rna_fit <- data.table::data.table(
#'   gene_symbol = c("FAKE_UP", "FAKE_DOWN"),
#'   rna_pattern = c("Early_Burst_Up", "Early_Loss_Down"),
#'   rna_pattern_rho = c(0.95, 0.92),
#'   lrt_padj = c(1e-4, 1e-3))
#' ev <- assemble_user_evidence(rna_fit = rna_fit)
#' compute_audit_score(evidence = ev)
#' @export
assemble_user_evidence <- function(rna_fit = NULL,
                                    prot_fit = NULL,
                                    scrna_summary = NULL,
                                    serum_summary = NULL,
                                    signal_peptide = NULL,
                                    cross_cohort_agreement = NULL,
                                    max_I2 = NULL) {
  parts <- list(rna_fit, prot_fit, scrna_summary, serum_summary)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0L) {
    stop("At least one input layer must be supplied.")
  }
  all_genes <- unique(unlist(lapply(parts, function(p) p$gene_symbol)))
  all_genes <- all_genes[!is.na(all_genes) & nzchar(all_genes)]
  if (length(all_genes) == 0L) {
    stop("No valid gene_symbol values found across input layers.")
  }

  out <- data.table::data.table(gene_symbol = all_genes)

  # RNA layer ------------------------------------------------------
  if (!is.null(rna_fit)) {
    rdt <- data.table::as.data.table(rna_fit)
    out <- merge(out,
                 rdt[, .(gene_symbol,
                         rna_pattern,
                         rna_pattern_rho,
                         rna_lrt_padj = lrt_padj)],
                 by = "gene_symbol", all.x = TRUE)
  } else {
    out[, rna_pattern := NA_character_]
    out[, rna_pattern_rho := NA_real_]
    out[, rna_lrt_padj := NA_real_]
  }

  # Protein layer --------------------------------------------------
  if (!is.null(prot_fit)) {
    pdt <- data.table::as.data.table(prot_fit)
    out <- merge(out,
                 pdt[, .(gene_symbol, prot_pattern)],
                 by = "gene_symbol", all.x = TRUE)
  } else {
    out[, prot_pattern := NA_character_]
  }

  # scRNA layer ----------------------------------------------------
  if (!is.null(scrna_summary)) {
    sdt <- data.table::as.data.table(scrna_summary)
    if (!"cell_origin_top" %in% names(sdt)) {
      stop("scrna_summary must contain `cell_origin_top` column.")
    }
    out <- merge(out,
                 sdt[, .(gene_symbol, cell_origin_top)],
                 by = "gene_symbol", all.x = TRUE)
  } else {
    out[, cell_origin_top := NA_character_]
  }

  # Serum layer ----------------------------------------------------
  if (!is.null(serum_summary)) {
    udt <- data.table::as.data.table(serum_summary)
    keep <- c("gene_symbol", "serum_detected", "serum_log2fc_PDAC_vs_HC")
    keep <- intersect(keep, names(udt))
    out <- merge(out, udt[, keep, with = FALSE],
                 by = "gene_symbol", all.x = TRUE)
    if (!"serum_detected" %in% names(out)) {
      out[, serum_detected := NA]
    }
    if (!"serum_log2fc_PDAC_vs_HC" %in% names(out)) {
      out[, serum_log2fc_PDAC_vs_HC := NA_real_]
    }
  } else {
    out[, serum_detected := NA]
    out[, serum_log2fc_PDAC_vs_HC := NA_real_]
  }

  # Signal peptide flag -------------------------------------------
  if (!is.null(signal_peptide)) {
    out[, flt_signal_peptide := gene_symbol %in% signal_peptide]
  } else {
    out[, flt_signal_peptide := NA]
  }

  # Tissue-to-serum direction match --------------------------------
  out[, flt_direction_match := .audit_user_dir_match(
    rna_pattern, serum_log2fc_PDAC_vs_HC)]

  # Cross-cohort agreement ----------------------------------------
  if (!is.null(cross_cohort_agreement)) {
    out[, rna_cohort_agreement := cross_cohort_agreement[gene_symbol]]
  } else {
    out[, rna_cohort_agreement := ifelse(!is.na(rna_pattern), 1.0,
                                           NA_real_)]
  }

  # Heterogeneity (I²) --------------------------------------------
  if (!is.null(max_I2)) {
    out[, max_I2_meta := max_I2[gene_symbol]]
  } else {
    out[, max_I2_meta := 0]
  }

  # placeholder columns that .audit_compute_features may reference
  # but that user evidence doesn't compute meta_* columns for. The
  # rna_weak rescue logic in v0.4.0 reads rna_pattern / rna_pattern_rho
  # / max_abs_beta_meta directly (no longer keyed off confidence_tier),
  # so user evidence with NA in these columns will trigger rescue
  # eligibility consistently.
  if (!"max_abs_beta_meta" %in% names(out)) {
    out[, max_abs_beta_meta := NA_real_]
  }
  if (!"rna_pattern_rho" %in% names(out)) {
    out[, rna_pattern_rho := NA_real_]
  }

  attr(out, "source") <- "user-supplied evidence (assemble_user_evidence)"
  out
}

# Internal helper — direction match between RNA pattern and serum log2FC
.audit_user_dir_match <- function(rna_pattern, serum_lfc) {
  out <- rep(NA, length(rna_pattern))
  up_pat <- c("Early_Burst_Up", "Early_Peak")
  down_pat <- c("Early_Loss_Down", "Early_Trough")
  ix <- !is.na(rna_pattern) & !is.na(serum_lfc)
  out[ix & rna_pattern %in% up_pat] <-
    serum_lfc[ix & rna_pattern %in% up_pat] > 0
  out[ix & rna_pattern %in% down_pat] <-
    serum_lfc[ix & rna_pattern %in% down_pat] < 0
  out
}
