#' Algorithm 1 — TRACE-D: Tissue-to-serum Directional Evidence Concordance
#'
#' `compute_trace_d()` deterministically types each gene by the direction
#' concordance between tissue (RNA + Protein) and serum evidence, with
#' an explicit pancreatitis-specificity annotation that separates
#' PDAC-specific signal from shared-inflammation signal. The four
#' classes ("A" direction-preserved, "B" direction-inverted, "C"
#' decoupled, `NA` insufficient tissue evidence) coexist with the
#' legacy `translation_class` column rather than replacing it
#' (back-compatible).
#'
#' This is one of two domain-specific audit algorithms introduced in
#' pdactrace v0.99.19; see [compute_pareto_layers()] for Algorithm 2.
#' Neither algorithm introduces new statistical machinery — the
#' contribution is the domain-specific combination of direction
#' consensus, magnitude thresholds, and pancreatitis-specificity
#' annotation packaged as a deterministic, reproducible audit layer.
#'
#' Pancreatitis-specificity is **annotation only** by design: a
#' Class A candidate flagged `shared_inflammation` is not demoted —
#' downstream consumers apply stricter filters via
#' [list_candidates()] if required.
#'
#' @param atlas Optional atlas (defaults to the bundled
#'   `pdactrace_reference`). Must contain `gene_symbol`, `rna_pattern`,
#'   `prot_pattern`, `max_abs_beta_meta`, `serum_log2fc_PDAC_vs_HC`,
#'   `serum_log2fc_Pan_vs_HC`, `serum_detected`, `flt_signal_peptide`.
#' @param tau_tissue Numeric; minimum `|max_abs_beta_meta|` for a
#'   tissue signal to count. Default `0.5`.
#' @param tau_serum Numeric; minimum `|serum_log2fc_PDAC_vs_HC|` for a
#'   serum signal to count. Default `0.1`.
#' @param weights Named numeric vector of length 3 summing to 1
#'   (tissue, serum, specificity components of `tracd_confidence`).
#'   Defaults to `c(tissue = 0.50, serum = 0.30, specificity = 0.20)`.
#' @return A data.table with one row per gene in `atlas` and the
#'   columns `tracd_tissue_dir`, `tracd_serum_dir`, `tracd_class`,
#'   `tracd_confidence`, `tracd_pancreatitis_overlap_score`,
#'   `tracd_pancreatitis_specificity`, `tracd_tissue_weight`,
#'   `tracd_decision_path`.
#' @examples
#' head(compute_trace_d())
#' compute_trace_d()[gene_symbol %in% c("LGALS3BP", "LTBP1", "GAPDH")]
#' @export
compute_trace_d <- function(atlas = NULL,
                            tau_tissue = 0.5,
                            tau_serum = 0.1,
                            weights = c(tissue = 0.50, serum = 0.30,
                                         specificity = 0.20)) {
  if (length(tau_tissue) != 1L || !is.finite(tau_tissue) ||
        tau_tissue < 0) {
    stop("`tau_tissue` must be a single non-negative numeric.")
  }
  if (length(tau_serum) != 1L || !is.finite(tau_serum) ||
        tau_serum < 0) {
    stop("`tau_serum` must be a single non-negative numeric.")
  }
  needed_weights <- c("tissue", "serum", "specificity")
  if (!all(needed_weights %in% names(weights))) {
    stop("`weights` must be a named numeric with elements: ",
         paste(needed_weights, collapse = ", "))
  }
  if (abs(sum(weights[needed_weights]) - 1) > 1e-6) {
    stop("`weights` must sum to 1 across (tissue, serum, specificity).")
  }

  ref <- .get_reference(atlas)
  required <- c("gene_symbol", "rna_pattern", "prot_pattern",
                "max_abs_beta_meta", "serum_log2fc_PDAC_vs_HC",
                "serum_log2fc_Pan_vs_HC", "serum_detected")
  missing_cols <- setdiff(required, names(ref))
  if (length(missing_cols) > 0L) {
    stop("atlas is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  dt <- data.table::as.data.table(ref)
  n <- nrow(dt)

  # ─── Direction consensus (tissue side) ──────────────────────
  rna_dir <- .audit_direction_of(dt$rna_pattern)
  prot_dir <- .audit_direction_of(dt$prot_pattern)
  rna_present <- !is.na(rna_dir)
  prot_present <- !is.na(prot_dir)
  both_agree <- rna_present & prot_present & (rna_dir == prot_dir)
  tissue_weight <- as.integer(rna_present) + as.integer(prot_present)
  tissue_weight[rna_present & prot_present & !both_agree] <- 0L
  tissue_dir <- rep(NA_character_, n)
  tissue_dir[both_agree] <- rna_dir[both_agree]
  only_rna <- rna_present & !prot_present
  tissue_dir[only_rna] <- rna_dir[only_rna]
  only_prot <- prot_present & !rna_present
  tissue_dir[only_prot] <- prot_dir[only_prot]

  # ─── Tissue / serum signal presence ─────────────────────────
  tissue_mag <- ifelse(is.na(dt$max_abs_beta_meta), 0,
                       dt$max_abs_beta_meta)
  tissue_signal <- !is.na(tissue_dir) & tissue_mag >= tau_tissue
  serum_detected <- !is.na(dt$serum_detected) & dt$serum_detected
  serum_pdac <- ifelse(is.na(dt$serum_log2fc_PDAC_vs_HC), 0,
                       dt$serum_log2fc_PDAC_vs_HC)
  serum_signal <- serum_detected & abs(serum_pdac) >= tau_serum

  # ─── Class assignment ────────────────────────────────────────
  serum_dir <- rep(NA_character_, n)
  serum_dir[serum_signal & serum_pdac > 0] <- "UP"
  serum_dir[serum_signal & serum_pdac < 0] <- "DOWN"

  tracd_class <- rep(NA_character_, n)
  cls_C <- tissue_signal & !serum_signal
  cls_A <- tissue_signal & serum_signal & (serum_dir == tissue_dir)
  cls_B <- tissue_signal & serum_signal & (serum_dir != tissue_dir) &
    !is.na(serum_dir)
  tracd_class[cls_C] <- "C"
  tracd_class[cls_A] <- "A"
  tracd_class[cls_B] <- "B"

  # ─── Pancreatitis specificity (annotation only) ─────────────
  serum_pan <- ifelse(is.na(dt$serum_log2fc_Pan_vs_HC), 0,
                      dt$serum_log2fc_Pan_vs_HC)
  eps <- 1e-6
  overlap_raw <- abs(serum_pan) / pmax(abs(serum_pdac), eps)
  sign_match <- (sign(serum_pdac) == sign(serum_pan)) & serum_pdac != 0
  overlap_score <- rep(NA_real_, n)
  overlap_score[serum_signal & sign_match] <-
    pmin(1, overlap_raw[serum_signal & sign_match])
  overlap_score[serum_signal & !sign_match] <- 0

  specificity <- rep(NA_character_, n)
  specificity[!serum_signal] <- "ambiguous"
  in_pdac_specific <- serum_signal &
    (!is.na(overlap_score) & overlap_score < 0.5)
  specificity[in_pdac_specific] <- "pdac_specific"
  in_shared <- serum_signal &
    (!is.na(overlap_score) & overlap_score >= 0.5) & sign_match
  specificity[in_shared] <- "shared_inflammation"
  specificity[is.na(specificity)] <- "ambiguous"

  # ─── Translation confidence ─────────────────────────────────
  tissue_score <- pmin(1, tissue_mag / 1.0)
  serum_score <- pmin(1, abs(serum_pdac) / 1.0)
  specificity_score <- ifelse(is.na(overlap_score), 0.5,
                              1 - overlap_score)
  confidence <- rep(NA_real_, n)
  confidence[cls_A | cls_B] <-
    weights[["tissue"]]      * tissue_score[cls_A | cls_B] +
    weights[["serum"]]       * serum_score[cls_A | cls_B] +
    weights[["specificity"]] * specificity_score[cls_A | cls_B]
  confidence[cls_C] <- weights[["tissue"]] * tissue_score[cls_C]
  confidence[!tissue_signal] <- NA_real_

  # ─── Audit-trail decision_path string ────────────────────────
  path_tissue <- ifelse(is.na(tissue_dir), "tissue_NA",
                        paste0("tissue_", tissue_dir))
  path_serum <- ifelse(serum_signal,
                       paste0("serum_", serum_dir),
                       ifelse(serum_detected, "serum_subthr",
                              "serum_absent"))
  path_class <- ifelse(is.na(tracd_class), "->NA",
                       paste0("->", tracd_class))
  path_pan <- paste0(",pan_", specificity)
  decision_path <- paste0(path_tissue, "+", path_serum, path_class,
                          path_pan)

  data.table::data.table(
    gene_symbol = dt$gene_symbol,
    tracd_tissue_dir = tissue_dir,
    tracd_serum_dir = serum_dir,
    tracd_class = tracd_class,
    tracd_confidence = confidence,
    tracd_pancreatitis_overlap_score = overlap_score,
    tracd_pancreatitis_specificity = specificity,
    tracd_tissue_weight = tissue_weight,
    tracd_decision_path = decision_path
  )
}
