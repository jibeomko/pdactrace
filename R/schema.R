#' pdactrace reference DB schema (v0.4.0, 12-template catalog with Early × 4 surface)
#'
#' Single source-of-truth for the column structure of
#' `pdactrace_reference`. v0.4.0 widened the trajectory classifier to
#' a 12-template competitive catalog (Early × 4 + Mid × 4 + Late × 2 +
#' Monotonic × 2). The atlas surface still restricts visible
#' `rna_pattern` / `prot_pattern` calls to **Early × 4 only** — Mid /
#' Late / Monotonic best-matches are flagged via `excluded_mid_pattern`,
#' `excluded_late_pattern`, `excluded_monotonic_pattern` for transparent
#' provenance.
#'
#' Design scope: pdactrace is a *resectable-stage marker discovery*
#' atlas. Including all 12 templates in the matching step is a
#' deliberate design choice — a gene must out-compete 11 alternatives
#' (Mid + Late + Monotonic best-matches) before it is surfaced as Early.
#' This widens the negative-evidence pool (e.g., GAPDH and CDH13 are
#' correctly demoted to Monotonic_Up under v0.4.0) and strengthens the
#' resectable-stage relevance of surfaced Early calls.
#'
#' Layers:
#' * **identifiers** - gene_symbol, ensembl_id, entrez_id
#' * **rna**         - DESeq2 LRT padj + per-stage Wald padj (v0.4.0) +
#'                     4 stage coefficients + 12-template best-match
#'                     (surfaced only if Early × 4) + Pearson rho +
#'                     Stouffer + 3 `excluded_*_pattern` flags
#' * **protein**     - tissue protein 12-template (surfaced if Early × 4) +
#'                     Tier1/Tier2 + RNA-protein concordance flag
#' * **scrna**       - dominant cell origin + full distribution +
#'                     hypergeom enrichment padj
#' * **serum**       - detected flag + log2FC PDAC/Pan/HC +
#'                     translation_class A/B/C
#' * **clinical**    - resectable_marker / panel_member flags
#' * **provenance**  - source phase script names + last_updated +
#'                     evidence_scope
#' * **t3_ready**    - detectability / direction labels + features +
#'                     model_trainable + model_card_ref (NA in T2.5)
#'
#' @format A `data.table` with one row per gene and the following
#'   columns (T2.5 = populated / T3 = NA placeholder):
#' \describe{
#'   \item{gene_symbol}{HGNC standard symbol (chr) (T2.5)}
#'   \item{ensembl_id}{ENSG identifier (chr) (T2.5)}
#'   \item{entrez_id}{NCBI gene id (int) (T2.5)}
#'   \item{rna_lrt_padj}{DESeq2 LRT BH-adjusted p-value (num) (T2.5)}
#'   \item{rna_beta_N, rna_beta_E, rna_beta_M, rna_beta_L}{
#'     log2FC stage coefficients vs Normal (num) (T2.5)}
#'   \item{rna_pattern}{Surfaced atlas call: Early × 4 best-match label
#'     (chr) or NA if the 12-template best-match was Mid / Late /
#'     Monotonic / Unclassified (T2.5)}
#'   \item{rna_pattern_rho}{Pearson rho with best 12-template match
#'     (num) (T2.5)}
#'   \item{rna_pattern_rho_runner_up}{rho with second-best template
#'     (num) (T2.5)}
#'   \item{rna_stouffer_z}{multi-cohort Stouffer Z (num) (T2.5)}
#'   \item{rna_cohort_agreement}{fraction of cohorts agreeing (num) (T2.5)}
#'   \item{excluded_mid_pattern}{TRUE if 12-template best-match was Mid_*
#'     (excluded from surface; recorded for transparency) (lgl) (T2.5)}
#'   \item{excluded_late_pattern}{TRUE if 12-template best-match was Late_*
#'     (v0.4.0) (lgl) (T2.5)}
#'   \item{excluded_monotonic_pattern}{TRUE if 12-template best-match was
#'     Monotonic_* (v0.4.0) (lgl) (T2.5)}
#'   \item{prot_pattern}{Surfaced atlas call: protein Early × 4 label or NA
#'     (chr) (T2.5)}
#'   \item{prot_tier}{Tier1 gold / Tier2 silver / NA (chr) (T2.5)}
#'   \item{rnaprot_concordant}{RNA pattern == Protein pattern (lgl) (T2.5)}
#'   \item{cell_origin_top}{dominant scRNA cell type (chr) (T2.5)}
#'   \item{cell_origin_distrib}{list-column of full % distribution
#'     (list) (T2.5)}
#'   \item{cell_origin_padj}{hypergeometric enrichment padj (num) (T2.5)}
#'   \item{serum_detected}{detected in >=1 of 3 serum cohorts (lgl) (T2.5)}
#'   \item{serum_n_cohorts_detected}{0-3 (int) (T2.5)}
#'   \item{serum_log2fc_PDAC_vs_HC}{log2FC PDAC vs HC (num) (T2.5)}
#'   \item{serum_log2fc_Pan_vs_HC}{log2FC Pancreatitis vs HC (num) (T2.5)}
#'   \item{translation_class}{"A" / "B" / "C" / NA (chr) (T2.5)}
#'   \item{phase77_strict}{member of 22 strict candidates (lgl) (T2.5)}
#'   \item{resectable_marker}{member of phase29 25 markers (lgl) (T2.5)}
#'   \item{panel_member}{member of LTBP1+SERPINA1 / CA19-9 hybrid (lgl) (T2.5)}
#'   \item{provenance}{comma-separated source phase scripts (chr) (T2.5)}
#'   \item{last_updated}{atlas snapshot date (Date) (T2.5)}
#'   \item{evidence_scope}{"tissue_only" / "tissue_serum" /
#'     "panel_validated" (chr) (T2.5)}
#'   \item{detectability_label}{T3 supervised label (lgl) (NA in T2.5)}
#'   \item{detectability_features}{T3 feature list-column (list) (NA in T2.5)}
#'   \item{serum_direction_label}{T3 direction label (Class A/C only)
#'     (chr) (NA in T2.5)}
#'   \item{model_trainable}{T3 training-set inclusion flag (lgl)
#'     (NA in T2.5)}
#'   \item{audit_*}{v0.3.0 3-axis + 2-gate audit score, five-label
#'     `audit_class` (`high_confidence`, `supported_uncertain`,
#'     `penalized`, `excluded`, `low`), Monte Carlo uncertainty, and
#'     MC `confidence_class` (mixed) (populated)}
#'   \item{model_card_ref}{T3 model-card identifier (chr) (NA in T2.5)}
#' }
#' @keywords internal
"pdactrace_reference"

#' Schema column specification
#'
#' Returns the canonical column list with type and T2.5/T3 status.
#' Used by `data-raw/build_reference.R` to allocate columns and by
#' `tests/testthat/test-schema.R` to verify completeness.
#'
#' @return data.table with columns: name, type, layer, t25_status.
#' @examples
#' head(schema_spec())
#' @export
schema_spec <- function() {
  data.table::data.table(
    name = c(
      "gene_symbol", "ensembl_id", "entrez_id",
      # rna (26): per-stage beta + per-stage LRT-padj/Wald-padj/lfcSE +
      # pattern + Stouffer + per-cohort + 3 exclusion flags (Mid/Late/Mono).
      # NOTE: rna_padj_E/M/L are DESeq2 LRT padj copies (LRT-mode results()
      # returns the same omnibus LRT padj for every per-contrast call). Use
      # rna_lrt_padj as the significance gate. v0.4.0 added per-contrast Wald
      # padj via nbinomWaldTest() refit (rna_wald_padj_E/M/L) and the 12-
      # template catalog (E×4 + M×4 + L×2 + Mono×2); rna_pattern surfaces
      # only Early × 4, with Mid/Late/Monotonic flagged via excluded_*_pattern.
      "rna_lrt_padj", "rna_beta_N", "rna_beta_E", "rna_beta_M", "rna_beta_L",
      "rna_padj_E", "rna_padj_M", "rna_padj_L",
      # v0.4.0: genuine per-contrast Wald padj from nbinomWaldTest refit.
      "rna_wald_padj_E", "rna_wald_padj_M", "rna_wald_padj_L",
      "rna_lfcSE_E", "rna_lfcSE_M", "rna_lfcSE_L",
      "rna_pattern", "rna_pattern_rho", "rna_pattern_rho_runner_up",
      "rna_stouffer_z", "rna_stouffer_p", "rna_stouffer_padj",
      "rna_cohort_agreement", "rna_per_cohort_trend",
      "rna_per_cohort_monotonic", "excluded_mid_pattern",
      "excluded_late_pattern", "excluded_monotonic_pattern",
      # protein (3)
      "prot_pattern", "prot_tier", "rnaprot_concordant",
      # scrna (4)
      "cell_origin_top", "cell_origin_distrib", "cell_origin_padj",
      "cell_specificity_tau",
      # serum (6)
      "serum_detected", "serum_n_cohorts_detected",
      "serum_log2fc_PDAC_vs_HC", "serum_log2fc_Pan_vs_HC",
      "translation_class", "phase77_strict",
      # clinical (3)
      "resectable_marker", "resectable_pattern_phase29", "panel_member",
      # filter_status (7)
      "flt_signal_peptide", "flt_serum_measurable", "flt_serum_significant",
      "flt_pancreatitis_pdac", "flt_pancreatitis_hc", "flt_direction_match",
      "flt_final",
      # annotation (8): pool stats + pancreatitis stats + serum raw means + t-test
      "ann_pool_logfc", "ann_pool_padj",
      "ann_pan_vs_hc_logfc", "ann_pan_vs_hc_pval", "ann_pan_excluded_phase60",
      "ann_pdac_mean", "ann_pan_mean", "ann_hc_mean",
      "ann_pdac_vs_pan_pval",
      # provenance (3)
      "provenance", "last_updated", "evidence_scope",
      # meta_analysis (17 - v0.2.0 random-effects meta + composite tier)
      "meta_NvE_beta", "meta_NvE_se", "meta_NvE_pval", "meta_NvE_padj",
      "meta_NvE_I2", "meta_NvE_k",
      "meta_MvE_beta", "meta_MvE_pval", "meta_MvE_padj",
      "meta_MvE_I2", "meta_MvE_k",
      "meta_LvE_beta", "meta_LvE_pval", "meta_LvE_padj",
      "meta_LvE_I2", "meta_LvE_k",
      "meta_cohort_divergent",
      # tier (2 - audit-framework dependencies; the v0.2.0
      # confidence_tier / early_onset_score / heterogeneity_factor
      # columns were removed in v0.4.0 because v0.3.0 audit_class
      # supersedes them and dual classification systems caused
      # reviewer confusion (see NEWS.md).
      "max_abs_beta_meta", "max_I2_meta",
      # audit (27 - v0.3.0 3+2 framework: 3 axes + 2 gates + 4 classes,
      #              7-feature internals retained for supplement S1)
      "audit_score_layer", "audit_score_direction", "audit_score_early",
      "audit_score_serum", "audit_score_rescue",
      "audit_evidence_strength", "audit_biological_coherence",
      "audit_translational_relevance",
      "audit_leakage_gate", "audit_heterogeneity_gate",
      "audit_positive_score", "audit_leakage_mult",
      "audit_heterogeneity_mult", "audit_score_raw", "audit_score",
      "audit_class",
      "audit_is_housekeeping", "audit_is_plasma_high_abundance",
      "audit_rescue_eligible",
      "audit_score_median", "audit_score_lo95", "audit_score_hi95",
      "audit_uncertainty_width",
      "audit_rank_median", "audit_rank_lo95", "audit_rank_hi95",
      "audit_confidence_class",
      # t3_ready (3)
      "serum_direction_label", "direction_model_trainable",
      "direction_model_card_ref"),
    type = c(
      "character", "character", "integer",
      # rna (26)
      "numeric", "numeric", "numeric", "numeric", "numeric",
      "numeric", "numeric", "numeric",
      "numeric", "numeric", "numeric",  # rna_wald_padj_E/M/L (v0.4.0)
      "numeric", "numeric", "numeric",
      "character", "numeric", "numeric",
      "numeric", "numeric", "numeric",
      "numeric", "list", "list", "logical",
      "logical", "logical",  # excluded_late_pattern, excluded_monotonic_pattern (v0.4.0)
      # protein (3)
      "character", "character", "logical",
      # scrna (4)
      "character", "list", "numeric", "numeric",
      # serum (6)
      "logical", "integer", "numeric", "numeric",
      "character", "logical",
      # clinical (3)
      "logical", "character", "logical",
      # filter_status (7)
      "logical", "logical", "logical", "logical", "logical", "logical",
      "logical",
      # annotation (9)
      "numeric", "numeric",
      "numeric", "numeric", "logical",
      "numeric", "numeric", "numeric",
      "numeric",
      # provenance (3)
      "character", "Date", "character",
      # meta_analysis (17)
      "numeric", "numeric", "numeric", "numeric", "numeric", "integer",
      "numeric", "numeric", "numeric", "numeric", "integer",
      "numeric", "numeric", "numeric", "numeric", "integer",
      "logical",
      # tier (2)
      "numeric", "numeric",
      # audit (27)
      "numeric", "numeric", "numeric", "numeric", "numeric",
      "numeric", "numeric", "numeric",
      "numeric", "numeric",
      "numeric", "numeric", "numeric", "numeric", "numeric",
      "character",
      "logical", "logical", "logical",
      "numeric", "numeric", "numeric", "numeric",
      "numeric", "numeric", "numeric", "character",
      # t3_ready (3)
      "character", "logical", "character"),
    layer = c(
      rep("identifier", 3),
      rep("rna", 26),
      rep("protein", 3),
      rep("scrna", 4),
      rep("serum", 6),
      rep("clinical", 3),
      rep("filter_status", 7),
      rep("annotation", 9),
      rep("provenance", 3),
      rep("meta_analysis", 17),
      rep("tier", 2),
      rep("audit", 27),
      rep("t3_ready", 3)),
    t25_status = c(
      rep("populated", 3),
      rep("populated", 26),
      rep("populated", 3),
      rep("populated", 4),
      rep("populated", 6),
      rep("populated", 3),
      rep("populated", 7),
      rep("populated", 9),
      rep("populated", 3),
      rep("populated", 17),
      rep("populated", 2),
      rep("populated", 27),
      rep("na_placeholder", 3))
  )
}

#' Early-onset pattern names (the 4 atlas-surfaced patterns)
#'
#' The 4 patterns surfaced by `pdactrace_reference$rna_pattern`. Under
#' v0.4.0 the trajectory classifier matches against 12 templates total
#' (Early × 4 + Mid × 4 + Late × 2 + Monotonic × 2) but only Early-onset
#' best-matches are surfaced — non-Early best-matches are flagged via
#' `excluded_mid_pattern`, `excluded_late_pattern`,
#' `excluded_monotonic_pattern`.
#'
#' @return character(4)
#' @examples
#' early_pattern_names()
#' early_patterns()  # short alias
#' @export
early_pattern_names <- function() {
  c("Early_Burst_Up", "Early_Loss_Down", "Early_Peak", "Early_Trough")
}

#' @rdname early_pattern_names
#' @export
early_patterns <- function() early_pattern_names()

#' Mid-onset pattern names (excluded from atlas surface)
#'
#' Returned for transparency. These patterns are *not* surfaced in
#' `pdactrace_reference$rna_pattern`; genes whose canonical phase33
#' 4-Early-pattern call falls in this set will have `rna_pattern = NA`
#' and `excluded_mid_pattern = TRUE`.
#'
#' @return character(4)
#' @examples
#' mid_pattern_names_excluded()
#' mid_patterns()  # short alias
#' @export
mid_pattern_names_excluded <- function() {
  c("Mid_Peak", "Mid_Trough", "Mid_Plateau_Up", "Mid_Plateau_Down")
}

#' @rdname mid_pattern_names_excluded
#' @export
mid_patterns <- function() mid_pattern_names_excluded()
