#' List candidate genes by biology-level filter parameters
#'
#' Multi-parameter filter on the pdactrace reference atlas. Pass `NULL`
#' (default) for any parameter to leave it unconstrained. Returns a
#' ranked `data.table` with key evidence columns.
#'
#' @param onset Character. `"Early"` to restrict to the 4 Early-onset
#'   patterns (the surfaced atlas universe). `NULL` (default) =
#'   unconstrained, which equals onset = "Early" because the atlas
#'   surface only exposes Early × 4 calls (under the v0.4.0 12-template
#'   classifier, Mid / Late / Monotonic best-matches are flagged via
#'   `excluded_*_pattern` columns and not returned by this function).
#' @param tissue_direction Character. `"Up"` keeps `Early_Burst_Up` +
#'   `Early_Peak`; `"Down"` keeps `Early_Loss_Down` + `Early_Trough`.
#'   `NULL` = both.
#' @param pattern Character vector. Specific 4-Early pattern names
#'   (e.g., `c("Early_Burst_Up", "Early_Peak")`). Overrides
#'   `tissue_direction` if both are set.
#' @param translation_class Character. One of `"concordant"` (Class A),
#'   `"inverse"` (Class B), `"decoupled"` (Class C). Aliases `"A"`,
#'   `"B"`, `"C"` are also accepted. `NULL` = any.
#' @param serum_detected Logical. Restrict to genes detected (or not)
#'   in any serum cohort. `NULL` = any.
#' @param signal_peptide Logical. Filter on `flt_signal_peptide`
#'   (phase60 SignalP pass). `NULL` = any.
#' @param panel_member Logical. Member of LTBP1+SERPINA1 / hybrid
#'   panels. `NULL` = any.
#' @param resectable_marker Logical. Member of phase29 25 markers (19
#'   in atlas). `NULL` = any.
#' @param phase77_strict Logical. Member of 22 strict tissue+serum
#'   convergent candidates. `NULL` = any.
#' @param min_audit_class Minimum `audit_class` to retain. One of
#'   `"high_confidence"`, `"supported_uncertain"`, `"penalized"`,
#'   `"low"`, `"excluded"`, or `"ALL"` (no filter). Defaults to
#'   `"ALL"`. Hierarchy (high to low):
#'   high_confidence, supported_uncertain, penalized, low, excluded.
#' @param top_n Optional integer. Return top N rows ranked by
#'   `audit_score` descending.
#' @return `data.table` with `gene_symbol`, key evidence columns, and
#'   `provenance`.
#' @examples
#'   list_candidates(onset = "Early", tissue_direction = "Up")
#'   list_candidates(translation_class = "inverse")
#'   list_candidates(translation_class = "concordant",
#'                   serum_detected = TRUE,
#'                   signal_peptide = TRUE,
#'                   top_n = 20)
#' @export
list_candidates <- function(onset = NULL,
                              tissue_direction = NULL,
                              pattern = NULL,
                              translation_class = NULL,
                              serum_detected = NULL,
                              signal_peptide = NULL,
                              panel_member = NULL,
                              resectable_marker = NULL,
                              phase77_strict = NULL,
                              min_audit_class = c("ALL",
                                                    "high_confidence",
                                                    "supported_uncertain",
                                                    "penalized",
                                                    "low",
                                                    "excluded"),
                              top_n = NULL) {
  min_audit_class <- match.arg(min_audit_class)
  ref <- .get_reference()
  keep <- rep(TRUE, nrow(ref))

  # -- Onset (atlas surface = Early × 4; spec-parity placeholder) --
  if (!is.null(onset)) {
    onset <- match.arg(onset, c("Early"))
    # Restriction implicit because rna_pattern surfaces Early × 4 only.
  }

  # -- Pattern + tissue_direction -----------------------------
  if (!is.null(pattern)) {
    pat_ok <- match.arg(pattern, early_pattern_names(),
                          several.ok = TRUE)
    keep <- keep & ref$rna_pattern %in% pat_ok &
             !is.na(ref$rna_pattern)
  } else if (!is.null(tissue_direction)) {
    tissue_direction <- match.arg(tissue_direction, c("Up", "Down"))
    pats <- if (tissue_direction == "Up")
      c("Early_Burst_Up", "Early_Peak") else
      c("Early_Loss_Down", "Early_Trough")
    keep <- keep & ref$rna_pattern %in% pats &
             !is.na(ref$rna_pattern)
  } else {
    # Default to atlas universe: any 4-Early
    keep <- keep & !is.na(ref$rna_pattern)
  }

  # -- Translation class (with aliases) -----------------------
  if (!is.null(translation_class)) {
    cls_map <- c(concordant = "A", inverse = "B", decoupled = "C",
                 A = "A", B = "B", C = "C")
    cls <- cls_map[translation_class]
    if (is.na(cls)) stop(sprintf(
      "Unknown translation_class: %s. Use 'concordant'/'inverse'/'decoupled'.",
      translation_class))
    keep <- keep & !is.na(ref$translation_class) &
             ref$translation_class == cls
  }

  # -- Logical filters ----------------------------------------
  if (!is.null(serum_detected)) keep <- keep &
    isTRUE_vec_match(ref$serum_detected, serum_detected)
  if (!is.null(signal_peptide))   keep <- keep &
    isTRUE_vec_match(ref$flt_signal_peptide, signal_peptide)
  if (!is.null(panel_member))     keep <- keep &
    isTRUE_vec_match(ref$panel_member, panel_member)
  if (!is.null(resectable_marker)) keep <- keep &
    isTRUE_vec_match(ref$resectable_marker, resectable_marker)
  if (!is.null(phase77_strict))   keep <- keep &
    isTRUE_vec_match(ref$phase77_strict, phase77_strict)

  # -- audit_class filter (v0.3.0) -----------------------------
  class_order <- c("ALL"                 = 0L,
                    "excluded"            = 1L,
                    "low"                 = 2L,
                    "penalized"           = 3L,
                    "supported_uncertain" = 4L,
                    "high_confidence"     = 5L)
  if (min_audit_class != "ALL") {
    min_rank <- class_order[min_audit_class]
    gene_rank <- class_order[ref$audit_class]
    keep <- keep & !is.na(gene_rank) & gene_rank >= min_rank
  }

  hit <- ref[which(keep)]

  cols <- c("gene_symbol", "audit_class", "audit_score",
            "rna_pattern", "rna_pattern_rho",
            "max_abs_beta_meta", "max_I2_meta",
            "rna_lrt_padj", "prot_pattern", "prot_tier",
            "translation_class", "serum_detected",
            "flt_signal_peptide", "flt_final",
            "resectable_marker", "panel_member",
            "evidence_scope", "provenance")
  cols <- intersect(cols, names(hit))
  out <- hit[, cols, with = FALSE]
  # Order by audit_score (descending) — best candidates first
  if ("audit_score" %in% names(out)) {
    data.table::setorder(out, -audit_score, na.last = TRUE)
  } else {
    data.table::setorder(out, -rna_pattern_rho, na.last = TRUE)
  }
  if (!is.null(top_n)) out <- utils::head(out, top_n)
  out
}

#' @keywords internal
isTRUE_vec_match <- function(col, target) {
  if (is.null(col)) return(rep(FALSE, length(col)))
  if (is.null(target)) return(rep(TRUE, length(col)))
  if (isTRUE(target))  return(!is.na(col) & col == TRUE)
  if (isFALSE(target)) return(!is.na(col) & col == FALSE)
  stop("Logical filter must be TRUE / FALSE / NULL.")
}
