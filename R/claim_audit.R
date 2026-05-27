#' Classify tissue-to-blood biomarker claim tiers
#'
#' `classify_claim_tier()` separates evidence strength from claim
#' strength. A gene can have strong stage-aware tissue evidence while
#' still receiving a conservative blood-biomarker claim such as
#' `translation_unknown` or `confounded`.
#'
#' @param atlas Optional atlas-like data table. Defaults to the bundled
#'   `pdactrace_reference`.
#' @param genes Optional character vector of genes to return. Ranking and
#'   thresholds are computed on the supplied atlas, then filtered.
#' @param tau_tissue Minimum absolute tissue effect used by
#'   [compute_trace_d()].
#' @param tau_serum Minimum absolute serum log2 fold-change used by
#'   [compute_trace_d()] and the serum-signal flag. Claim-tier
#'   classification calls TRACE-D in `legacy_translation = "fallback"`
#'   mode so the bundled atlas' historical `translation_class` remains
#'   visible when serum log2FC is unavailable.
#' @return A data.table with one row per gene and claim audit columns.
#' @examples
#' classify_claim_tier(genes = c("LGALS3BP", "LTBP1", "ALB", "GAPDH"))
#' @export
classify_claim_tier <- function(atlas = NULL,
                                genes = NULL,
                                tau_tissue = 0.5,
                                tau_serum = 0.1) {
  ref <- data.table::as.data.table(.get_reference(atlas))
  required <- c("gene_symbol", "rna_pattern", "prot_pattern",
                "max_abs_beta_meta", "serum_detected",
                "serum_n_cohorts_detected", "serum_log2fc_PDAC_vs_HC",
                "serum_log2fc_Pan_vs_HC", "translation_class",
                "flt_signal_peptide",
                "audit_is_housekeeping",
                "audit_is_plasma_high_abundance")
  missing_cols <- setdiff(required, names(ref))
  if (length(missing_cols) > 0L) {
    stop("atlas is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  trace <- compute_trace_d(ref, tau_tissue = tau_tissue,
                           tau_serum = tau_serum,
                           legacy_translation = "fallback")
  score <- compute_audit_score(NULL, evidence = ref)
  score <- score[, .(gene_symbol,
                     audit_score_det = audit_score,
                     audit_class_det = audit_class,
                     leakage_gate,
                     heterogeneity_gate)]

  dt <- merge(ref, trace, by = "gene_symbol", all.x = TRUE)
  dt <- merge(dt, score, by = "gene_symbol", all.x = TRUE)

  tissue_mag <- ifelse(is.na(dt$max_abs_beta_meta), 0,
                       dt$max_abs_beta_meta)
  rna_supported <- !is.na(dt$rna_pattern)
  prot_supported <- !is.na(dt$prot_pattern)
  tissue_supported <- rna_supported | prot_supported |
    tissue_mag >= tau_tissue

  serum_observed <- (!is.na(dt$serum_detected) & dt$serum_detected) |
    (!is.na(dt$serum_n_cohorts_detected) &
       dt$serum_n_cohorts_detected > 0L) |
    (!is.na(dt$translation_class) & nzchar(dt$translation_class))
  serum_fc <- ifelse(is.na(dt$serum_log2fc_PDAC_vs_HC), 0,
                     dt$serum_log2fc_PDAC_vs_HC)
  serum_signal <- serum_observed &
    (abs(serum_fc) >= tau_serum |
       (!is.na(dt$translation_class) & nzchar(dt$translation_class)) |
       (!is.na(dt$tracd_class) & nzchar(dt$tracd_class)))
  serum_concordant <- dt$tracd_class == "A" | dt$translation_class == "A"
  serum_concordant[is.na(serum_concordant)] <- FALSE

  exportable_plausible <- (!is.na(dt$flt_signal_peptide) &
                             dt$flt_signal_peptide) |
    serum_observed
  exportable_plausible[is.na(exportable_plausible)] <- FALSE

  hard_excluded <- (!is.na(dt$leakage_gate) & dt$leakage_gate == 0) |
    (!is.na(dt$audit_class_det) & dt$audit_class_det == "excluded") |
    (!is.na(dt$audit_is_housekeeping) & dt$audit_is_housekeeping)
  hard_excluded[is.na(hard_excluded)] <- FALSE

  plasma_risk <- (!is.na(dt$leakage_gate) & dt$leakage_gate < 1) |
    (!is.na(dt$audit_class_det) & dt$audit_class_det == "penalized") |
    (!is.na(dt$audit_is_plasma_high_abundance) &
       dt$audit_is_plasma_high_abundance)
  plasma_risk[is.na(plasma_risk)] <- FALSE

  shared_inflammation <- dt$tracd_pancreatitis_specificity ==
    "shared_inflammation"
  shared_inflammation[is.na(shared_inflammation)] <- FALSE
  confounded <- !hard_excluded & (plasma_risk | shared_inflammation)

  confounder_risk <- data.table::fcase(
    hard_excluded, "hard_exclusion",
    plasma_risk & shared_inflammation, "plasma_and_inflammation",
    plasma_risk, "plasma_high_abundance",
    shared_inflammation, "shared_inflammation",
    default = "none")

  translation_status <- data.table::fcase(
    dt$tracd_class == "A" | dt$translation_class == "A",
    "direction_preserved",
    dt$tracd_class == "B" | dt$translation_class == "B",
    "direction_inverted",
    dt$tracd_class == "C" | dt$translation_class == "C",
    "decoupled_or_unobserved",
    serum_signal, "serum_observed_no_tissue_direction",
    tissue_supported, "translation_unknown",
    default = "insufficient_tissue_evidence")

  claim_tier <- data.table::fcase(
    hard_excluded, "excluded",
    confounded, "confounded",
    serum_concordant, "serum_concordant",
    serum_signal, "serum_observed",
    tissue_supported & exportable_plausible, "exportable_plausible",
    tissue_supported, "translation_unknown",
    default = "insufficient_tissue_evidence")

  claim_strength <- data.table::fcase(
    claim_tier == "serum_concordant", 6L,
    claim_tier == "serum_observed", 5L,
    claim_tier == "exportable_plausible", 4L,
    claim_tier == "translation_unknown", 3L,
    claim_tier == "confounded", 2L,
    claim_tier == "excluded", 1L,
    default = 0L)

  claim_reason <- paste0(
    "tissue=", ifelse(tissue_supported, "supported", "insufficient"),
    ";export=", ifelse(exportable_plausible, "plausible", "unknown"),
    ";serum=", ifelse(serum_signal, "signal",
                      ifelse(serum_observed, "observed_subthreshold",
                             "not_observed")),
    ";translation=", translation_status,
    ";risk=", confounder_risk)

  out <- data.table::data.table(
    gene_symbol = dt$gene_symbol,
    claim_tier = claim_tier,
    claim_strength = claim_strength,
    tissue_supported = tissue_supported,
    exportable_plausible = exportable_plausible,
    serum_observed = serum_observed,
    serum_signal = serum_signal,
    serum_concordant = serum_concordant,
    translation_status = translation_status,
    confounder_risk = confounder_risk,
    tracd_class = dt$tracd_class,
    tracd_confidence = dt$tracd_confidence,
    audit_class = dt$audit_class_det,
    audit_score = dt$audit_score_det,
    claim_reason = claim_reason)

  if (!is.null(genes)) {
    keep_genes <- .audit_genes(genes, reference = ref)
    out <- out[gene_symbol %in% keep_genes]
    out <- out[match(keep_genes, out$gene_symbol)]
  }
  out
}

#' Summarize the tissue-to-blood translation gap
#'
#' @inheritParams classify_claim_tier
#' @return A data.table with claim-tier counts and percentages.
#' @examples
#' summarize_translation_gap()
#' @export
summarize_translation_gap <- function(atlas = NULL,
                                      tau_tissue = 0.5,
                                      tau_serum = 0.1) {
  tiers <- classify_claim_tier(atlas = atlas, tau_tissue = tau_tissue,
                               tau_serum = tau_serum)
  total <- nrow(tiers)
  tiers[, .(n = .N,
            pct = .N / total,
            tissue_supported = sum(tissue_supported, na.rm = TRUE),
            serum_observed = sum(serum_observed, na.rm = TRUE),
            serum_concordant = sum(serum_concordant, na.rm = TRUE)),
        by = claim_tier][order(-n)]
}

#' Evaluate ranking robustness across plausible audit weights
#'
#' Samples weight vectors over the three audit axes and reports how
#' often each gene remains inside top-N sets. This is a sensitivity
#' analysis, not a new trained model.
#'
#' @param atlas Optional atlas-like data table.
#' @param genes Optional genes to return after atlas-wide ranking.
#' @param n_draws Number of sampled weight vectors.
#' @param top_n Integer vector of top-N cutoffs.
#' @param weight_floor Minimum weight assigned to each axis. The
#'   residual mass is sampled uniformly on the simplex.
#' @param seed Random seed.
#' @return A data.table with rank intervals and top-N stability columns.
#' @examples
#' \donttest{
#' run_weight_robustness(genes = c("LGALS3BP", "LTBP1"), n_draws = 50)
#' }
#' @export
run_weight_robustness <- function(atlas = NULL,
                                  genes = NULL,
                                  n_draws = 1000L,
                                  top_n = c(50L, 100L, 500L),
                                  weight_floor = 0.05,
                                  seed = 20260527L) {
  n_draws <- as.integer(n_draws)
  if (length(n_draws) != 1L || is.na(n_draws) || n_draws < 10L) {
    stop("`n_draws` must be a single integer >= 10.")
  }
  top_n <- sort(unique(as.integer(top_n)))
  if (length(top_n) == 0L || any(is.na(top_n)) || any(top_n < 1L)) {
    stop("`top_n` must contain positive integers.")
  }
  if (length(weight_floor) != 1L || !is.finite(weight_floor) ||
      weight_floor < 0 || weight_floor >= 1 / 3) {
    stop("`weight_floor` must be in [0, 1/3).")
  }

  ref <- data.table::as.data.table(.get_reference(atlas))
  score <- compute_audit_score(NULL, evidence = ref)
  X <- as.matrix(score[, .(evidence_strength,
                           biological_coherence,
                           translational_relevance)])
  gate <- score$leakage_gate * score$heterogeneity_gate
  n <- nrow(score)
  top_n <- pmin(top_n, n)

  withr::local_seed(seed)
  rank_mat <- matrix(NA_integer_, nrow = n, ncol = n_draws)
  top_counts <- matrix(0L, nrow = n, ncol = length(top_n))
  score_mean <- numeric(n)
  weights_seen <- matrix(NA_real_, nrow = n_draws, ncol = 3L)
  colnames(weights_seen) <- c("evidence_strength",
                              "biological_coherence",
                              "translational_relevance")
  residual <- 1 - 3 * weight_floor

  for (i in seq_len(n_draws)) {
    w_raw <- stats::rexp(3L)
    w <- weight_floor + residual * w_raw / sum(w_raw)
    weights_seen[i, ] <- w
    raw <- as.numeric(X %*% w) * gate
    mx <- max(raw, na.rm = TRUE)
    draw_score <- if (is.finite(mx) && mx > 0) raw / mx else rep(0, n)
    score_mean <- score_mean + draw_score / n_draws
    ord <- order(-draw_score, na.last = TRUE)
    rank <- integer(n)
    rank[ord] <- seq_along(ord)
    rank_mat[, i] <- rank
    for (j in seq_along(top_n)) {
      top_idx <- ord[seq_len(top_n[j])]
      top_counts[top_idx, j] <- top_counts[top_idx, j] + 1L
    }
  }

  out <- data.table::data.table(
    gene_symbol = score$gene_symbol,
    weight_score_mean = score_mean,
    weight_rank_median = apply(rank_mat, 1L, stats::median, na.rm = TRUE),
    weight_rank_lo95 = apply(rank_mat, 1L, stats::quantile,
                             probs = 0.025, na.rm = TRUE),
    weight_rank_hi95 = apply(rank_mat, 1L, stats::quantile,
                             probs = 0.975, na.rm = TRUE))
  for (j in seq_along(top_n)) {
    nm <- paste0("weight_top", top_n[j], "_stability")
    out[, (nm) := top_counts[, j] / n_draws]
  }
  data.table::setorder(out, weight_rank_median, -weight_score_mean)
  attr(out, "weights") <- weights_seen
  attr(out, "weight_floor") <- weight_floor

  if (!is.null(genes)) {
    keep_genes <- .audit_genes(genes, reference = ref)
    out <- out[gene_symbol %in% keep_genes]
    out <- out[match(keep_genes, out$gene_symbol)]
  }
  out
}

#' Classify Pareto evidence support
#'
#' Convenience wrapper around [compute_pareto_layers()] that converts
#' numeric layers into manuscript-facing evidence classes.
#'
#' @param atlas Optional atlas-like data table.
#' @param genes Optional genes to return.
#' @param top_n Candidate pool size passed to [compute_pareto_layers()].
#' @return A data.table with Pareto layer, rank, and class.
#' @examples
#' compute_pareto_class(genes = c("LGALS3BP", "LTBP1"), top_n = 200)
#' @export
compute_pareto_class <- function(atlas = NULL,
                                 genes = NULL,
                                 top_n = 2000L) {
  ref <- data.table::as.data.table(.get_reference(atlas))
  layers <- compute_pareto_layers(atlas = ref, top_n = top_n)
  score_vals <- ref$audit_score
  pool_idx <- order(-score_vals, na.last = NA)
  pool_idx <- pool_idx[seq_len(min(as.integer(top_n), length(pool_idx)))]
  in_pool <- ref$gene_symbol %in% ref$gene_symbol[pool_idx]
  out <- merge(data.table::data.table(gene_symbol = ref$gene_symbol,
                                      in_pareto_pool = in_pool),
               layers, by = "gene_symbol", all.x = TRUE)
  gate_excluded <- out$in_pareto_pool & out$pareto_excluded_by_gate
  gate_excluded[is.na(gate_excluded)] <- FALSE
  out[, pareto_class := data.table::fcase(
    !in_pareto_pool, "outside_top_n",
    gate_excluded, "gate_excluded",
    pareto_layer == 1L, "pareto_front",
    pareto_layer <= 3L, "pareto_supported",
    !is.na(pareto_layer), "pareto_lower",
    default = "outside_top_n")]
  if (!is.null(genes)) {
    keep_genes <- .audit_genes(genes, reference = ref)
    out <- out[gene_symbol %in% keep_genes]
    out <- out[match(keep_genes, out$gene_symbol)]
  }
  out[]
}
