#' Per-axis mathematical evidence for one gene
#'
#' Returns the underlying numerical evidence values that fed the
#' frozen v0.3.0 audit decisions, organised by axis. This is the
#' canonical machine-readable form sitting between the plain-English
#' phase-tag layer ([format_provenance()]) and the audit-score
#' decomposition ([explain_score()]). The point is to expose
#' *what the math actually says* -- not just labels, not just a
#' composite score.
#'
#' Every value reported here is either read directly from
#' [pdactrace_reference] (or [pdactrace_protein_betas] for the
#' protein side) or derived by a single named arithmetic operation
#' on those columns. No new fitting, no recomputation. The function
#' name `evidence_math` reflects that it is the mathematical view
#' of the same evidence captured by the verbose / provenance text
#' layers.
#'
#' @section Axes:
#'
#' \describe{
#'   \item{`trajectory_fit`}{Best 12-template Pearson rho, runner-up
#'     rho, and the specificity margin `delta_rho = rho_best -
#'     rho_runner_up`. A delta_rho above ~0.10 indicates clean
#'     template specificity (the gene's profile is closer to its
#'     matched template than to any other).}
#'   \item{`effect_magnitude`}{Euclidean norm and max-absolute of the
#'     per-stage beta vector, on both RNA and protein sides. The
#'     norm uses the 3-vector `(beta_E, beta_M, beta_L)` because
#'     `beta_N = 0` is a reference-level constant by DESeq2 / limma
#'     convention, not a measured effect.}
#'   \item{`cohort_consistency`}{Stouffer Z, BH-adjusted p, fraction
#'     of contributing cohorts agreeing with the meta sign, and the
#'     maximum heterogeneity I-squared across the three contrasts
#'     N-vs-E, M-vs-E, L-vs-E (which feeds the heterogeneity gate).}
#'   \item{`rna_protein_coupling`}{Cosine similarity between the
#'     RNA and tissue-protein beta vectors over `(E, M, L)`.
#'     `+1` means same template + same direction; `-1` means
#'     opposite. Reported alongside the categorical `prot_tier`
#'     and the boolean `rnaprot_concordant` flag.}
#'   \item{`serum_bridge`}{Translation class A / B / C
#'     (A = same-direction RNA<->serum, B = opposite, C =
#'     decoupled), serum log2 fold-changes against healthy and
#'     pancreatitis cohorts, number of detected cohorts, and the
#'     `phase77_strict` membership flag.}
#'   \item{`cell_specificity`}{Top scRNA cell of origin and the
#'     Yanai et al. tau index. tau approaches 1 when expression
#'     is concentrated in one cell type.}
#'   \item{`filter_survival`}{The 7-step phase60 audit as a count
#'     of passed steps plus a per-step boolean vector
#'     (signal_peptide, serum_measurable, serum_significant,
#'     pancreatitis_pdac, pancreatitis_hc, direction_match, final).}
#'   \item{`clinical_role`}{Whether the gene was screened in the
#'     phase29 resectable-marker shortlist and whether it is a
#'     member of the predeclared serum panel.}
#' }
#'
#' @param gene_symbol HGNC gene symbol (length-1 character).
#' @param reference Optional `data.table` to inject in place of
#'   the bundled atlas (used by tests and downstream pipelines).
#' @return Invisibly, a list with one element per axis (named-list
#'   of named scalars). Suitable for `str()` / programmatic use.
#'   The full 12-rho vector for non-Early best-match genes can be
#'   recovered with [score_trajectory()] using the bundled fit.
#' @examples
#' m <- evidence_math("LGALS3BP")
#' m$trajectory_fit
#' m$rna_protein_coupling$cosine
#' m$filter_survival$passed
#' @seealso [explain_gene()] for the formatted text view, and
#'   [compare_genes()] for a multi-gene tidy table.
#' @export
evidence_math <- function(gene_symbol, reference = NULL) {
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

  prot <- .em_protein_betas_row(target)

  out <- list(
    gene                 = gene_symbol,
    trajectory_fit       = .em_trajectory_fit(row),
    effect_magnitude     = .em_effect_magnitude(row, prot),
    cohort_consistency   = .em_cohort_consistency(row),
    rna_protein_coupling = .em_rna_protein_coupling(row, prot),
    serum_bridge         = .em_serum_bridge(row),
    cell_specificity     = .em_cell_specificity(row),
    filter_survival      = .em_filter_survival(row),
    clinical_role        = .em_clinical_role(row))

  attr(out, "reference_version") <- as.character(row$last_updated)
  out
}

# ---- per-axis builders --------------------------------------------------

.em_trajectory_fit <- function(row) {
  pat <- as.character(row$rna_pattern)
  rho <- as.numeric(row$rna_pattern_rho)
  ru  <- as.numeric(row$rna_pattern_rho_runner_up)
  drho <- if (is.na(rho) || is.na(ru)) NA_real_ else rho - ru
  note <- if (is.na(pat))
    "non-Early best match -- use score_trajectory() for the full 12-rho vector"
    else NA_character_
  list(rna_pattern   = pat,
       rho_best      = rho,
       rho_runner_up = ru,
       delta_rho     = drho,
       note          = note)
}

.em_effect_magnitude <- function(row, prot) {
  rna_b <- c(as.numeric(row$rna_beta_E),
             as.numeric(row$rna_beta_M),
             as.numeric(row$rna_beta_L))
  rna_norm <- if (any(is.na(rna_b))) NA_real_ else sqrt(sum(rna_b^2))
  rna_max  <- if (any(is.na(rna_b))) NA_real_ else max(abs(rna_b))
  rna_at   <- if (any(is.na(rna_b))) NA_character_
              else c("E", "M", "L")[which.max(abs(rna_b))]
  rna_target <- as.character(row$rna_pattern)
  rna_target_stage <- if (is.na(rna_target)) NA_character_
                       else if (grepl("^Early", rna_target)) "E"
                       else if (grepl("^Mid",   rna_target)) "M"
                       else if (grepl("^Late",  rna_target)) "L"
                       else "all"

  if (is.null(prot)) {
    prot_norm <- NA_real_; prot_max <- NA_real_; prot_at <- NA_character_
  } else {
    prot_b <- c(prot$prot_beta_E, prot$prot_beta_M, prot$prot_beta_L)
    prot_norm <- if (any(is.na(prot_b))) NA_real_ else sqrt(sum(prot_b^2))
    prot_max  <- if (any(is.na(prot_b))) NA_real_ else max(abs(prot_b))
    prot_at   <- if (any(is.na(prot_b))) NA_character_
                 else c("E", "M", "L")[which.max(abs(prot_b))]
  }

  list(rna_beta_norm     = rna_norm,
       rna_beta_max_abs  = rna_max,
       rna_max_at_stage  = rna_at,
       rna_target_stage  = rna_target_stage,
       prot_beta_norm    = prot_norm,
       prot_beta_max_abs = prot_max,
       prot_max_at_stage = prot_at)
}

.em_cohort_consistency <- function(row) {
  i2_vec <- c(as.numeric(row$meta_NvE_I2),
              as.numeric(row$meta_MvE_I2),
              as.numeric(row$meta_LvE_I2))
  max_i2 <- if (all(is.na(i2_vec))) NA_real_ else max(i2_vec, na.rm = TRUE)
  list(stouffer_z       = as.numeric(row$rna_stouffer_z),
       stouffer_p       = as.numeric(row$rna_stouffer_p),
       stouffer_padj    = as.numeric(row$rna_stouffer_padj),
       cohort_agreement = as.numeric(row$rna_cohort_agreement),
       max_meta_I2      = max_i2)
}

.em_rna_protein_coupling <- function(row, prot) {
  rna_b <- c(as.numeric(row$rna_beta_E),
             as.numeric(row$rna_beta_M),
             as.numeric(row$rna_beta_L))
  if (is.null(prot)) {
    return(list(cosine            = NA_real_,
                prot_pattern      = NA_character_,
                prot_tier         = as.character(row$prot_tier),
                rnaprot_concordant = as.logical(row$rnaprot_concordant),
                prot_in_atlas     = FALSE,
                note              = "gene absent from pdactrace_protein_betas"))
  }
  prot_b <- c(prot$prot_beta_E, prot$prot_beta_M, prot$prot_beta_L)
  cos_val <- NA_real_
  cos_note <- NA_character_
  if (!any(is.na(rna_b)) && !any(is.na(prot_b))) {
    nr <- sqrt(sum(rna_b^2)); np <- sqrt(sum(prot_b^2))
    if (nr == 0 || np == 0) {
      cos_note <- "zero beta vector"
    } else {
      cos_val <- as.numeric(sum(rna_b * prot_b) / (nr * np))
    }
  } else {
    cos_note <- "missing beta entries"
  }
  list(cosine             = cos_val,
       prot_pattern       = as.character(prot$prot_pattern_12),
       prot_pattern_rho   = as.numeric(prot$prot_pattern_rho),
       prot_tier          = as.character(row$prot_tier),
       rnaprot_concordant = as.logical(row$rnaprot_concordant),
       prot_in_atlas      = TRUE,
       note               = cos_note)
}

.em_serum_bridge <- function(row) {
  list(translation_class        = as.character(row$translation_class),
       serum_log2fc_PDAC_vs_HC  = as.numeric(row$serum_log2fc_PDAC_vs_HC),
       serum_log2fc_Pan_vs_HC   = as.numeric(row$serum_log2fc_Pan_vs_HC),
       serum_n_cohorts_detected = as.integer(row$serum_n_cohorts_detected),
       phase77_strict           = as.logical(row$phase77_strict))
}

.em_cell_specificity <- function(row) {
  list(cell_origin_top  = as.character(row$cell_origin_top),
       tau              = as.numeric(row$cell_specificity_tau),
       cell_origin_padj = as.numeric(row$cell_origin_padj))
}

.em_filter_survival <- function(row) {
  steps <- c(
    signal_peptide      = as.logical(row$flt_signal_peptide),
    serum_measurable    = as.logical(row$flt_serum_measurable),
    serum_significant   = as.logical(row$flt_serum_significant),
    pancreatitis_pdac   = as.logical(row$flt_pancreatitis_pdac),
    pancreatitis_hc     = as.logical(row$flt_pancreatitis_hc),
    direction_match     = as.logical(row$flt_direction_match),
    final               = as.logical(row$flt_final))
  passed <- sum(steps, na.rm = TRUE)
  list(passed   = as.integer(passed),
       total    = length(steps),
       per_step = steps)
}

.em_clinical_role <- function(row) {
  list(resectable_marker = isTRUE(as.logical(row$resectable_marker)),
       panel_member      = isTRUE(as.logical(row$panel_member)))
}

# ---- protein-betas accessor (cached, keyed by gene_symbol) -------------

.em_protein_betas_row <- function(target) {
  e <- new.env()
  utils::data("pdactrace_protein_betas", package = "pdactrace", envir = e)
  pb <- e$pdactrace_protein_betas
  if (!data.table::is.data.table(pb)) {
    pb <- data.table::as.data.table(pb)
  }
  hit <- pb[gene_symbol == target]
  if (nrow(hit) == 0L) return(NULL)
  as.list(hit[1L])
}
