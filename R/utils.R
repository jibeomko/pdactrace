#' Internal: vectorised isTRUE that handles NA gracefully
#'
#' Returns a logical vector of the same length as `x`. Values that are
#' `NA` are returned as `FALSE`. Use this when filtering data.table
#' rows by a logical column that may contain `NA`.
#'
#' @param x A logical (or coercible) vector.
#' @return Logical vector.
#' @keywords internal
isTRUE_vec <- function(x) {
  if (is.null(x)) return(logical(0))
  !is.na(x) & x == TRUE
}

utils::globalVariables(c(
  ":=", ".", ".SD",
  # align_patient_profile() column names
  "vote_share", "weighted_dist", "cor_to_stage_axis",
  "cor_pval", "cor_lo95", "cor_hi95", "n_genes_used",
  "audit_confidence_class", "audit_heterogeneity_mult",
  "audit_is_housekeeping", "audit_is_plasma_high_abundance",
  "audit_leakage_mult", "audit_positive_score", "audit_rank_hi95",
  "audit_rank_lo95", "audit_rank_median", "audit_rescue_eligible",
  "audit_score", "audit_score_direction", "audit_score_early",
  "audit_score_hi95", "audit_score_layer", "audit_score_lo95",
  "audit_score_median", "audit_score_raw", "audit_score_rescue",
  "audit_score_serum", "audit_uncertainty_width", "cell_origin_top",
  "celltype", "ci_hi", "ci_lo", "class_label", "class_route",
  "cohort", "confidence_class",
  "cross_cohort_concord", "cross_layer_concord", "dir_prot",
  "dir_rna", "evidence_tier", "flag",
  "flt_direction_match", "flt_final", "flt_pancreatitis_hc",
  "flt_pancreatitis_pdac", "flt_serum_measurable",
  "flt_serum_significant", "flt_signal_peptide", "gene",
  "gene_symbol", "group", "het_mult", "include_primary_eval",
  "include_secondary_eval", "is_early", "is_highlight", "is_hk",
  "is_meta", "is_plasma_hi", "label", "layer", "layer_count",
  "layer_prot", "layer_rna", "layer_scrna", "layer_serum",
  "leakage_mult", "log2FC", "lrt_padj", "lrt_sig_factor",
  "lrt_significant", "max_I2_meta", "mean_expression", "monotonic",
  "n_phase60_pass", "n_ref", "n_target", "other_layer_count",
  "panel_member", "pass", "pct_of_total", "pdactrace_call",
  "phase77_strict", "positive_score", "prot_pattern", "rank_hi95",
  "rank_lo95", "rank_median", "rescue_eligible",
  "resectable_marker", "rna_beta_E", "rna_beta_L", "rna_beta_M",
  "rna_cohort_agreement", "rna_lrt_padj", "rna_pattern",
  "rna_pattern_rho", "rna_pattern_rho_runner_up", "rna_weak",
  "score", "score_95ci", "score_direction", "score_early",
  "score_layer", "score_rescue", "score_serum", "se",
  "serum_detected", "serum_log2fc_PDAC_vs_HC", "sig_label",
  "stage", "status", "step", "tissue_dir", "tissue_effect",
  "tissue_signed", "translation_class", "trend", "trend_factor",
  "uncertainty_width", "underlying_metric", "v0.2_tier", "value",
  "variable", "x_pos"
))
