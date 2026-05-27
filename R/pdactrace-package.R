#' pdactrace . PDAC-TRACE: queryable stage-aware tissue-to-serum biomarker atlas
#'
#' **Central message**: *PDAC tissue biomarker is not always a serum-up
#' biomarker* - tissue-derived candidates can preserve, invert, or
#' decouple when projected into serum.
#'
#' pdactrace integrates public PDAC multi-omics data - bulk RNA-seq,
#' tissue proteomics, scRNA cell origin, serum proteomics, and
#' pancreatitis context - under a **12-template competitive trajectory
#' framework** (Normal -> Early -> Mid -> Late) where genes compete
#' against 11 alternative shapes (Early × 4 + Mid × 4 + Late × 2 +
#' Monotonic × 2) and only Early-onset best-matches are surfaced as
#' resectable-stage candidates. The current release adds a **3-axis
#' + 2-gate audit framework** (v0.3.0) and a **user-data extension
#' API + per-stage Wald padj** (v0.4.0).
#'
#' @section 12-template framework with Early × 4 surface:
#' Each gene's z-scored 4-point N/E/M/L trajectory is matched against
#' 12 pre-declared templates by Pearson-rho argmax. A template call is
#' retained only when the best-match rho passes `rho_cutoff`
#' (`0.85` by default in [classify_trajectory()]). Only the 4
#' Early-onset best-matches are surfaced in
#' `pdactrace_reference$rna_pattern`:
#' * `Early_Burst_Up` - Normal low -> Early up -> sustained
#' * `Early_Loss_Down` - Normal high -> Early down -> sustained
#' * `Early_Peak` - Early peak -> Mid/Late decline
#' * `Early_Trough` - Early trough -> Mid/Late recovery
#'
#' Non-Early best-matches (Mid × 4 / Late × 2 / Monotonic × 2) are
#' **excluded from the surface** but flagged via
#' `excluded_mid_pattern`, `excluded_late_pattern`,
#' `excluded_monotonic_pattern` for transparent provenance. Including
#' the 8 non-Early templates in the matching step is the v0.4.0
#' design discipline: a gene must beat 11 alternatives before being
#' surfaced — this widens the negative-evidence pool and correctly
#' demotes leaky candidates (e.g., GAPDH, CDH13 → Monotonic_Up).
#'
#' @section Core lookup API:
#' * [query_gene()]               - single-gene full evidence dump
#' * [query_panel()]              - multi-gene join
#' * [list_candidates()]          - criterion-based candidate filter
#' * [summarize_gene_evidence()]  - human-readable text summary
#' * [plot_gene_evidence()]       - multi-panel composite figure
#'
#' @section Trajectory scoring (user RNAseq input):
#' * [fit_stage_de()]         - DESeq2 LRT wrapper for stage-aware DE
#' * [classify_trajectory()]  - best-match against the 12-template catalog
#' * [score_trajectory()]     - gene-level rho vector vs all 12 templates
#' * [align_patient_profile()] - sample-level alignment of a single patient's
#'   log2(tumor/normal) profile against the atlas's stage-trajectory axes
#'   (rho-style readout, not a stage prediction; v0.4.1)
#'
#' @section Filter audit + visualization:
#' * [trace_filters()]            - phase60 7-step filter audit + class route
#' * [plot_filter_trace()]        - pass/fail step bar across genes
#' * [plot_panel_heatmap()]       - gene x evidence comparison heatmap
#' * [plot_candidate_landscape()] - tissue x serum scatter, Class A/B colored
#'
#' @section Class B inverse policy:
#' `translation_class == "B"` is treated as a **rare, manually
#' curated failure mode** of tissue-to-serum translation, **not** a
#' predicted class. Even in v0.2.0+, Class B will be flagged via
#' nearest-reference evidence rather than predicted.
#'
#' @section Citation:
#' Ko J. *pdactrace: a queryable stage-aware PDAC tissue-to-serum
#' biomarker reference atlas.* (2026)
#'
#' @docType package
#' @name pdactrace-package
#' @aliases pdactrace
#' @importFrom utils data head
#' @importFrom methods setGeneric setMethod is
#' @importFrom stats complete.cases
"_PACKAGE"

# Required by data.table when := and .SD are used inside package code.
# See vignette('datatable-importing') in data.table.
.datatable.aware <- TRUE

# Declare NSE variables used inside data.table expressions to satisfy
# R CMD check's static analysis.
utils::globalVariables(c(
  "audit_class", "audit_score_raw_v3", "audit_score_v3",
  "biological_coherence", "evidence_strength",
  "heterogeneity_gate", "leakage_gate",
  "positive_score_v3", "translational_relevance",
  "max_abs_beta_meta",
  # plot_gene_hexagon polygon-ring + spoke helpers
  "ring", "x", "y", "x_end", "y_end",
  "label_x", "label_y", "color",
  # compare_candidates / explain_score NSE
  "rna_pattern", "cell_origin_top", "audit_score",
  "redundancy_with",
  # data.table column-list NSE
  "..beta_cols", "..keep_cols", "..cols", "..feat_cols",
  "..flt_cols", ".N",
  # plot_template_atlas / plot_gene_template NSE + ggplot aes columns
  "stage", "z", "template", "n_cohort",
  "z_mean", "z_sd", "z_lo", "z_hi",
  # .template_aggregate per-stage z-score NSE
  "z_N", "z_E", "z_M", "z_L", "template_argmax",
  # .ft_serum_strip ggplot aes NSE
  "contrast", "class_lbl", "log2fc",
  # evidence_math / compare_genes NSE
  "gene_symbol", "axis", "metric", "value", "key",
  # evidence_features / anchor_similarity / evidence_model NSE
  "feature", "coef_value", "anchor_similarity", "evidence_tier",
  "include_primary_eval", "include_secondary_eval", "gene",
  "prot_beta_E", "prot_beta_M", "prot_beta_L", "prot_pattern_rho",
  "predicted_prob",
  # plot_filter_trace serum strip helper
  "serum_log2fc_Pan_vs_HC"
))
