#' Build a transparent evidence graph for one gene
#'
#' `build_evidence_graph()` exposes the audit-score evidence as a small
#' node/edge graph. It is intentionally descriptive: pdactrace v0.3.0 does
#' not train a graph neural network, and the graph is used to make the
#' hand-engineered evidence aggregation auditable.
#'
#' @param gene_symbol HGNC gene symbol.
#' @return A list with `nodes`, `edges`, and `score` data.tables, or `NULL`
#'   invisibly if the gene is not in the atlas.
#' @examples
#' g <- build_evidence_graph("LTBP1")
#' head(g$nodes)
#' g$score
#' @export
build_evidence_graph <- function(gene_symbol) {
  gene <- .audit_one_gene(gene_symbol)
  ref <- .get_reference()
  idx <- which(ref[["gene_symbol"]] == gene)
  row <- ref[idx]
  if (nrow(row) == 0L) {
    message(sprintf("No evidence for '%s' in pdactrace atlas.", gene))
    return(invisible(NULL))
  }

  feat <- extract_graph_features(gene)
  unc <- tryCatch(propagate_uncertainty(gene),
                  error = function(e) data.table::data.table())

  nodes <- data.table::data.table(
    node_id = c("gene", "layer_presence", "direction_agreement",
                "early_pattern", "serum_bridge", "rescue_signal",
                "leakage_gate", "heterogeneity_gate", "audit_score"),
    layer = c("gene", "aggregate", "aggregate", "rna", "serum",
              "cross_layer", "artifact", "cohort", "readout"),
    node_type = c("target", rep("positive_feature", 5),
                  "multiplier", "multiplier", "readout"),
    value_numeric = c(
      NA_real_,
      feat$score_layer,
      feat$score_direction,
      feat$score_early,
      feat$score_serum,
      feat$score_rescue,
      feat$leakage_mult,
      feat$het_mult,
      feat$audit_score
    ),
    weight = c(NA_real_, 0.20, 0.20, 0.20, 0.10, 0.10,
               NA_real_, NA_real_, NA_real_),
    value_label = c(
      gene,
      sprintf("%d/4 layers", feat$layer_count),
      sprintf("RNA/protein concordant: %s",
              ifelse(isTRUE(feat$cross_layer_concord), "yes", "no")),
      as.character(row$rna_pattern),
      ifelse(isTRUE(row$serum_detected), "serum detected", "not detected"),
      ifelse(isTRUE(feat$rescue_eligible), "eligible", "not eligible"),
      .audit_leakage_label(feat),
      .audit_heterogeneity_label(row$max_I2_meta, feat$het_mult),
      sprintf("audit_score %.3f", feat$audit_score)
    ),
    available = c(TRUE, TRUE, TRUE, !is.na(row$rna_pattern),
                  !is.na(row$serum_detected), TRUE, TRUE, TRUE, TRUE)
  )

  edges <- data.table::data.table(
    from = c("layer_presence", "direction_agreement", "early_pattern",
             "serum_bridge", "rescue_signal", "leakage_gate",
             "heterogeneity_gate"),
    to = rep("audit_score", 7),
    edge_type = c(rep("positive_support", 5), "artifact_gate",
                  "heterogeneity_gate"),
    multiplier = c(rep(NA_real_, 5), feat$leakage_mult, feat$het_mult)
  )

  score <- compute_audit_score(gene)
  if (nrow(unc) > 0L) {
    score <- merge(score, unc, by = "gene_symbol", all.x = TRUE)
  }

  out <- list(gene_symbol = gene, nodes = nodes, edges = edges,
              score = score)
  class(out) <- c("pdactrace_evidence_graph", "list")
  out
}

#' Extract graph-aware audit features
#'
#' Returns the frozen v0.3.0 audit feature components used by
#' `compute_audit_score()`. These are deterministic summaries of the
#' evidence graph, not learned embeddings.
#'
#' @param gene_symbol Character vector of HGNC gene symbols.
#' @return A data.table with one row per matched gene.
#' @examples
#' extract_graph_features(c("LGALS3BP", "LTBP1", "GAPDH"))
#' @export
extract_graph_features <- function(gene_symbol) {
  genes <- .audit_genes(gene_symbol)
  ref <- .get_reference()
  all_feat <- .audit_feature_table(ref)
  out <- all_feat[all_feat$gene_symbol %in% genes]
  out[match(genes[genes %in% out$gene_symbol], out$gene_symbol)]
}

#' Compute the frozen pdactrace audit score (3-axis + 2-gate)
#'
#' Aggregates the seven implementation features into three transparent
#' evidence axes (`evidence_strength`, `biological_coherence`,
#' `translational_relevance`) and applies two pre-specified reliability
#' gates (`leakage_gate`, `heterogeneity_gate`) before computing
#' `audit_score = positive_score * leakage_gate * heterogeneity_gate`.
#' Each gene receives one of four `audit_class` labels:
#' `high_confidence`, `supported_uncertain`, `penalized`, or `excluded`.
#'
#' This is not a supervised biomarker predictor; the formula is
#' pre-specified and external anchors are used only post-freeze for
#' enrichment evaluation. The seven implementation features are kept as
#' supplementary detail; see `extract_graph_features()`.
#'
#' @param gene_symbol Optional character vector of HGNC gene symbols.
#'   When `NULL` (default) returns the full atlas.
#' @param evidence Optional user-supplied evidence table (e.g.,
#'   produced by [assemble_user_evidence()]) to score in place of the
#'   bundled atlas. When `NULL` (default), the frozen v0.3.0 reference
#'   atlas is used; the framework discipline (frozen weights, gates)
#'   applies identically to either input.
#' @return A data.table with `evidence_strength`,
#'   `biological_coherence`, `translational_relevance`,
#'   `leakage_gate`, `heterogeneity_gate`, `positive_score`,
#'   `audit_score`, and `audit_class`.
#' @examples
#' compute_audit_score(c("LGALS3BP", "LTBP1", "ALB", "GAPDH"))
#' @export
compute_audit_score <- function(gene_symbol = NULL, evidence = NULL) {
  ref <- .get_reference(evidence)
  feat <- .audit_feature_table(ref)

  # 3-axis aggregates from internal 7-feature components
  feat[, evidence_strength := pmin(1, score_layer + 0.5 * score_rescue)]
  feat[, biological_coherence := (score_direction + score_early) / 2]
  feat[, translational_relevance := score_serum]
  feat[, leakage_gate := leakage_mult]
  feat[, heterogeneity_gate := het_mult]

  # 3+2 frozen formula (locked weights)
  feat[, positive_score_v3 := 0.40 * evidence_strength +
                                0.35 * biological_coherence +
                                0.25 * translational_relevance]
  feat[, audit_score_raw_v3 := positive_score_v3 *
                                 leakage_gate * heterogeneity_gate]
  mx <- max(feat$audit_score_raw_v3, na.rm = TRUE)
  if (is.finite(mx) && mx > 0) {
    feat[, audit_score_v3 := audit_score_raw_v3 / mx]
  } else {
    feat[, audit_score_v3 := 0]
  }

  # 4-class assignment (locked rules; gates take precedence over score)
  feat[, audit_class := data.table::fcase(
    leakage_gate == 0, "excluded",
    leakage_gate < 1, "penalized",
    heterogeneity_gate < 1, "supported_uncertain",
    audit_score_v3 >= 0.5, "high_confidence",
    audit_score_v3 >= 0.3, "supported_uncertain",
    default = "low")]

  out <- feat[, .(gene_symbol,
                  evidence_strength,
                  biological_coherence,
                  translational_relevance,
                  leakage_gate,
                  heterogeneity_gate,
                  positive_score = positive_score_v3,
                  audit_score = audit_score_v3,
                  audit_class)]
  if (is.null(gene_symbol)) return(out)
  genes <- .audit_genes(gene_symbol, reference = evidence)
  res <- out[gene_symbol %in% genes]
  res[match(genes, res$gene_symbol)]
}

#' Propagate audit-score uncertainty
#'
#' By default this returns the stored Monte Carlo uncertainty summary from
#' the atlas. If `n_mc` is supplied, the uncertainty is recomputed over the
#' full atlas using the frozen scoring rule and perturbing evidence values
#' only.
#'
#' @param gene_symbol Character vector of HGNC gene symbols.
#' @param n_mc Optional integer number of Monte Carlo iterations.
#' @param seed Random seed used when `n_mc` is supplied.
#' @param evidence Optional user-supplied evidence table (see
#'   [assemble_user_evidence()]). When `NULL` (default), the bundled
#'   reference atlas is used. When supplied without `n_mc`, the
#'   stored MC summary is unavailable and `n_mc` must be set.
#' @return A data.table with score CI, rank CI, uncertainty width, and class.
#' @examples
#' propagate_uncertainty(c("LTBP1", "THBS2", "GAPDH"))
#' @export
propagate_uncertainty <- function(gene_symbol, n_mc = NULL, seed = 42,
                                   evidence = NULL) {
  genes <- .audit_genes(gene_symbol, reference = evidence)
  ref <- .get_reference(evidence)

  if (is.null(n_mc)) {
    stored <- .audit_uncertainty_table(ref)
    if (is.null(stored)) {
      stop("Stored audit uncertainty columns are unavailable. ",
           "Call propagate_uncertainty(..., n_mc = 500) to recompute.")
    }
    out <- stored[stored$gene_symbol %in% genes]
    return(out[match(genes[genes %in% out$gene_symbol], out$gene_symbol)])
  }

  n_mc <- as.integer(n_mc)
  if (length(n_mc) != 1L || is.na(n_mc) || n_mc < 10L) {
    stop("n_mc must be a single integer >= 10.")
  }
  mc <- .audit_mc_table(ref, n_mc = n_mc, seed = seed)
  out <- mc[mc$gene_symbol %in% genes]
  out[match(genes[genes %in% out$gene_symbol], out$gene_symbol)]
}

#' Evaluate external-anchor enrichment in top-ranked genes
#'
#' External anchors are used only for evaluation. They are not used by the
#' frozen v0.3.0 scoring rule.
#'
#' @param top_n Integer vector of top-N cutoffs.
#' @param tier Anchor tier: `primary`, `secondary`, `exploratory`, or `all`.
#' @param score_col Ranking column. Defaults to deterministic `audit_score`.
#' @param evidence Optional user-supplied evidence table (see
#'   [assemble_user_evidence()]). When `NULL` (default), the bundled
#'   reference atlas is used.
#' @param anchors Optional user-supplied anchor table with columns
#'   `gene`, `include_primary_eval`, and `include_secondary_eval`.
#'   When `NULL` (default), the bundled `pdactrace_external_anchors`
#'   table is used.
#' @param ... Forwarded to `evaluate_anchor_enrichment()` when invoked
#'   through the `anchor_enrichment()` short alias.
#' @return A data.table with hypergeometric enrichment results.
#' @examples
#' evaluate_anchor_enrichment(top_n = 100, tier = "secondary")
#' @export
evaluate_anchor_enrichment <- function(top_n = c(50, 100, 200, 500, 1000),
                                       tier = c("primary", "secondary",
                                                "exploratory", "all"),
                                       score_col = c("audit_score",
                                                     "audit_score_median",
                                                     "pareto_rank",
                                                     "pareto_stability_top1",
                                                     "tracd_confidence"),
                                       evidence = NULL,
                                       anchors = NULL) {
  tier <- match.arg(tier, several.ok = TRUE)
  score_col <- match.arg(score_col)
  if ("all" %in% tier) tier <- c("primary", "secondary", "exploratory")

  ref <- .get_reference(evidence)
  unc <- .audit_uncertainty_table(ref)
  score_v3 <- compute_audit_score(NULL, evidence = evidence)
  if (score_col == "audit_score_median" && !is.null(unc)) {
    score_dt <- merge(score_v3[, .(gene_symbol)],
                       unc[, .(gene_symbol, audit_score_median)],
                       by = "gene_symbol", all.x = TRUE)
  } else if (score_col == "pareto_rank") {
    if (!"pareto_rank" %in% names(ref)) {
      stop("Atlas does not carry pareto_rank. ",
           "Rebuild via data-raw/attach_pareto_columns.R or supply ",
           "an atlas containing pareto_rank.")
    }
    score_dt <- merge(score_v3[, .(gene_symbol)],
                       data.table::as.data.table(ref)[
                         , .(gene_symbol, pareto_rank)],
                       by = "gene_symbol", all.x = TRUE)
  } else if (score_col == "pareto_stability_top1") {
    if (!"pareto_stability_top1" %in% names(ref)) {
      stop("Atlas does not carry pareto_stability_top1. ",
           "Rebuild via data-raw/attach_pareto_columns.R.")
    }
    score_dt <- merge(score_v3[, .(gene_symbol)],
                       data.table::as.data.table(ref)[
                         , .(gene_symbol, pareto_stability_top1)],
                       by = "gene_symbol", all.x = TRUE)
  } else if (score_col == "tracd_confidence") {
    if (!"tracd_confidence" %in% names(ref)) {
      stop("Atlas does not carry tracd_confidence. ",
           "Rebuild via data-raw/attach_trace_d_columns.R or supply ",
           "an atlas containing tracd_confidence.")
    }
    score_dt <- merge(score_v3[, .(gene_symbol)],
                       data.table::as.data.table(ref)[
                         , .(gene_symbol, tracd_confidence)],
                       by = "gene_symbol", all.x = TRUE)
  } else {
    score_dt <- score_v3[, .(gene_symbol, audit_score)]
  }
  data.table::setnames(score_dt, names(score_dt)[2L], "score")
  if (score_col == "pareto_rank") {
    # Lower pareto_rank is better; NA (excluded) ranks last.
    data.table::setorder(score_dt, score, na.last = TRUE)
  } else {
    data.table::setorder(score_dt, -score, na.last = TRUE)
  }

  anchor_dt <- if (is.null(anchors)) {
    .audit_anchor_table()
  } else {
    if (!data.table::is.data.table(anchors)) {
      anchors <- data.table::as.data.table(anchors)
    }
    required_cols <- c("gene", "include_primary_eval",
                       "include_secondary_eval")
    if (!all(required_cols %in% names(anchors))) {
      stop("`anchors` must contain columns: ",
           paste(required_cols, collapse = ", "))
    }
    anchors
  }
  tier_sets <- list(
    primary = anchor_dt[include_primary_eval == TRUE, gene],
    secondary = anchor_dt[include_secondary_eval == TRUE, gene],
    exploratory = anchor_dt$gene
  )
  tier_labels <- c(primary = "PRIMARY (T1_validated, direct)",
                   secondary = "SECONDARY (T1+T2)",
                   exploratory = "EXPLORATORY (all)")

  n_atlas <- nrow(score_dt)
  top_n <- sort(unique(as.integer(top_n)))
  if (any(is.na(top_n)) || any(top_n < 1L)) {
    stop("top_n must contain positive integers.")
  }
  top_n <- pmin(top_n, n_atlas)

  res <- list()
  k <- 0L
  for (tr in tier) {
    anchors_declared <- unique(tier_sets[[tr]])
    anchors_in <- intersect(anchors_declared, score_dt$gene_symbol)
    n_in <- length(anchors_in)
    for (n in top_n) {
      top_set <- score_dt$gene_symbol[seq_len(n)]
      hits <- sum(top_set %in% anchors_in)
      expect <- n * n_in / n_atlas
      fold <- hits / max(expect, 1e-12)
      pval <- stats::phyper(hits - 1, n_in, n_atlas - n_in, n,
                            lower.tail = FALSE)
      k <- k + 1L
      res[[k]] <- data.table::data.table(
        tier = tier_labels[[tr]],
        tier_key = tr,
        top_n = n,
        n_anchor_declared = length(anchors_declared),
        n_anchor_in_atlas = n_in,
        hits = hits,
        expect = expect,
        fold = fold,
        pval = pval
      )
    }
  }
  data.table::rbindlist(res)
}

#' @rdname evaluate_anchor_enrichment
#' @export
anchor_enrichment <- function(...) evaluate_anchor_enrichment(...)

#' Summarize audit-score case-study genes
#'
#' Combines v0.2 tier, v0.3 score, uncertainty, anchor provenance, and
#' artifact/rescue calls for manuscript case-study panels.
#'
#' @param genes Character vector of HGNC gene symbols.
#' @return A data.table with one row per matched gene.
#' @examples
#' case_study(c("THBS2", "LTBP1", "ALB", "GAPDH"))
#' @export
case_study <- function(genes) {
  genes <- .audit_genes(genes)
  ref <- .get_reference()
  feat <- extract_graph_features(genes)
  score <- compute_audit_score(genes)
  unc <- tryCatch(propagate_uncertainty(genes),
                  error = function(e) data.table::data.table(gene_symbol = genes))
  anchors <- .audit_anchor_table()
  anchor_tier <- anchors[, .(anchor_tier = paste(unique(evidence_tier),
                                                 collapse = ";")),
                         by = gene]

  rows <- ref[ref$gene_symbol %in% genes]
  rows <- rows[match(genes[genes %in% rows$gene_symbol], rows$gene_symbol)]
  out <- data.table::data.table(
    gene_symbol = rows$gene_symbol,
    rna_pattern = rows$rna_pattern,
    prot_pattern = rows$prot_pattern,
    scrna_origin = rows$cell_origin_top,
    serum_support = rows$translation_class,
    max_I2_meta = rows$max_I2_meta
  )
  out <- merge(out, score, by = "gene_symbol", all.x = TRUE)
  out <- merge(out,
               feat[, .(gene_symbol, is_hk, is_plasma_hi, rescue_eligible)],
               by = "gene_symbol", all.x = TRUE)
  out <- merge(out, unc, by = "gene_symbol", all.x = TRUE)
  out <- merge(out, anchor_tier, by.x = "gene_symbol", by.y = "gene",
               all.x = TRUE)
  out[is.na(anchor_tier), anchor_tier := "none"]
  out[, score_95ci := sprintf("[%.2f, %.2f]",
                              audit_score_lo95, audit_score_hi95)]
  out[, pdactrace_call := .audit_call(anchor_tier, is_hk, is_plasma_hi,
                                      rescue_eligible)]
  out[order(match(gene_symbol, genes))]
}

#' @export
print.pdactrace_evidence_graph <- function(x, ...) {
  cat(sprintf("%s evidence graph: %d nodes, %d edges\n",
              x$gene_symbol, nrow(x$nodes), nrow(x$edges)))
  if ("audit_score" %in% names(x$score)) {
    cat(sprintf("audit_score: %.3f\n", x$score$audit_score[1]))
  }
  invisible(x)
}

# Internal helpers ---------------------------------------------------------

.audit_one_gene <- function(gene_symbol) {
  if (!is.character(gene_symbol) || length(gene_symbol) != 1L ||
      is.na(gene_symbol) || !nzchar(gene_symbol)) {
    stop("gene_symbol must be a single non-empty character string.")
  }
  gene_symbol
}

.audit_genes <- function(genes, reference = NULL) {
  if (!is.character(genes) || length(genes) == 0L) {
    stop("genes must be a non-empty character vector.")
  }
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  if (length(genes) == 0L) stop("No valid gene symbols supplied.")
  ref <- .get_reference(reference)
  missing <- setdiff(genes, ref$gene_symbol)
  if (length(missing) > 0L) {
    message(sprintf("Skipped %d gene(s) not in reference: %s",
                    length(missing), paste(missing, collapse = ", ")))
  }
  intersect(genes, ref$gene_symbol)
}

.audit_feature_table <- function(ref) {
  prefixed <- c("audit_score_layer", "audit_score_direction",
                "audit_score_early", "audit_score_serum",
                "audit_score_rescue", "audit_positive_score",
                "audit_leakage_mult", "audit_heterogeneity_mult",
                "audit_score_raw", "audit_score",
                "audit_is_housekeeping",
                "audit_is_plasma_high_abundance",
                "audit_rescue_eligible")
  if (all(prefixed %in% names(ref))) {
    out <- data.table::as.data.table(ref)[, .(
      gene_symbol,
      score_layer = audit_score_layer,
      score_direction = audit_score_direction,
      score_early = audit_score_early,
      score_serum = audit_score_serum,
      score_rescue = audit_score_rescue,
      positive_score = audit_positive_score,
      leakage_mult = audit_leakage_mult,
      het_mult = audit_heterogeneity_mult,
      audit_score_raw,
      audit_score,
      is_hk = audit_is_housekeeping,
      is_plasma_hi = audit_is_plasma_high_abundance,
      rescue_eligible = audit_rescue_eligible
    )]
  } else {
    out <- .audit_compute_features(ref)
  }

  ref_dt <- data.table::as.data.table(ref)
  out[, layer_count :=
        as.integer(!is.na(ref_dt$rna_pattern)) +
        as.integer(!is.na(ref_dt$prot_pattern)) +
        as.integer(!is.na(ref_dt$cell_origin_top)) +
        as.integer(!is.na(ref_dt$serum_detected) & ref_dt$serum_detected)]
  out[, cross_layer_concord := .audit_direction_of(ref_dt$rna_pattern) ==
        .audit_direction_of(ref_dt$prot_pattern) &
        !is.na(.audit_direction_of(ref_dt$rna_pattern)) &
        !is.na(.audit_direction_of(ref_dt$prot_pattern))]
  out
}

.audit_compute_features <- function(ref) {
  dt <- data.table::copy(data.table::as.data.table(ref))
  early <- early_pattern_names()
  dt[, layer_rna := !is.na(rna_pattern)]
  dt[, layer_prot := !is.na(prot_pattern)]
  dt[, layer_scrna := !is.na(cell_origin_top)]
  dt[, layer_serum := !is.na(serum_detected) & serum_detected == TRUE]
  dt[, score_layer := (as.integer(layer_rna) + as.integer(layer_prot) +
                         as.integer(layer_scrna) +
                         as.integer(layer_serum)) / 4]
  dt[, dir_rna := .audit_direction_of(rna_pattern)]
  dt[, dir_prot := .audit_direction_of(prot_pattern)]
  dt[, cross_layer_concord := !is.na(dir_rna) & !is.na(dir_prot) &
       dir_rna == dir_prot]
  dt[, cross_cohort_concord := pmax(0, pmin(1, .audit_na0(rna_cohort_agreement)))]
  dt[, score_direction := 0.5 * cross_cohort_concord +
       0.5 * as.integer(cross_layer_concord)]
  dt[, is_early := !is.na(rna_pattern) & rna_pattern %in% early]
  dt[, lrt_sig_factor := pmin(1, -log10(pmax(rna_lrt_padj, 1e-10,
                                             na.rm = TRUE)) / 4)]
  dt[, score_early := ifelse(is_early & !is.na(rna_pattern_rho),
                             rna_pattern_rho * .audit_na0(lrt_sig_factor),
                             0)]
  dt[, score_early := pmax(0, pmin(1, score_early))]
  dt[, score_serum := 0.4 * as.integer(layer_serum) +
       0.3 * as.integer(!is.na(flt_signal_peptide) &
                          flt_signal_peptide == TRUE) +
       0.3 * as.integer(!is.na(flt_direction_match) &
                          flt_direction_match == TRUE)]
  # rna_weak: gene whose RNA evidence alone is insufficient. Triggers
  # rescue eligibility when other layers (Prot + scRNA + serum) compensate.
  # v0.4.0: replaced legacy confidence_tier check with direct RNA-evidence
  # criteria (low effect size OR weak template fit OR non-Early surface).
  dt[, rna_weak := is.na(rna_pattern) |
       (!is.na(max_abs_beta_meta) & max_abs_beta_meta < 0.585) |
       (!is.na(rna_pattern_rho) & rna_pattern_rho < 0.85)]
  dt[, other_layer_count := as.integer(layer_prot) +
       as.integer(layer_scrna) + as.integer(layer_serum)]
  dt[, rescue_eligible := rna_weak & other_layer_count >= 2L &
       cross_layer_concord]
  dt[, score_rescue := as.integer(rescue_eligible)]
  dt[, positive_score := 0.20 * score_layer +
       0.20 * score_direction +
       0.20 * score_early +
       0.10 * score_serum +
       0.10 * score_rescue]

  dt[, is_hk := .audit_existing_flag(dt, "audit_is_housekeeping")]
  dt[, is_plasma_hi := .audit_existing_flag(dt,
                                            "audit_is_plasma_high_abundance")]
  if ("audit_leakage_mult" %in% names(dt)) {
    dt[, leakage_mult := audit_leakage_mult]
  } else {
    dt[, leakage_mult := data.table::fcase(is_hk, 0.00,
                                           is_plasma_hi, 0.50,
                                           default = 1.00)]
  }
  dt[, het_mult := data.table::fcase(
    is.na(max_I2_meta), 1.00,
    max_I2_meta < 50, 1.00,
    max_I2_meta < 70, 1.00,
    max_I2_meta < 90, 0.70,
    default = 0.30
  )]
  dt[, audit_score_raw := positive_score * leakage_mult * het_mult]
  mx <- max(dt$audit_score_raw, na.rm = TRUE)
  dt[, audit_score := audit_score_raw / mx]
  dt[, .(gene_symbol, score_layer, score_direction, score_early,
         score_serum, score_rescue, positive_score, leakage_mult, het_mult,
         audit_score_raw, audit_score, is_hk, is_plasma_hi, rescue_eligible)]
}

.audit_uncertainty_table <- function(ref) {
  cols <- c("audit_score_median", "audit_score_lo95", "audit_score_hi95",
            "audit_uncertainty_width", "audit_rank_median",
            "audit_rank_lo95", "audit_rank_hi95", "audit_confidence_class")
  if (!all(cols %in% names(ref))) return(NULL)
  data.table::as.data.table(ref)[, .(
    gene_symbol,
    audit_score_median,
    audit_score_lo95,
    audit_score_hi95,
    uncertainty_width = audit_uncertainty_width,
    rank_median = audit_rank_median,
    rank_lo95 = audit_rank_lo95,
    rank_hi95 = audit_rank_hi95,
    confidence_class = audit_confidence_class
  )]
}

.audit_mc_table <- function(ref, n_mc, seed) {
  # Scope the seed to this helper only (BiocCheck-compliant; withr
  # restores the global RNG state when the function exits).
  withr::local_seed(seed)
  dt <- data.table::copy(data.table::as.data.table(ref))
  feat <- .audit_feature_table(dt)
  n_genes <- nrow(dt)
  score_mc <- matrix(NA_real_, nrow = n_genes, ncol = n_mc)
  rank_mc <- matrix(NA_real_, nrow = n_genes, ncol = n_mc)

  rho_b <- dt$rna_pattern_rho
  agree_b <- dt$rna_cohort_agreement
  i2_b <- dt$max_I2_meta
  lrt_b <- pmin(1, -log10(pmax(dt$rna_lrt_padj, 1e-10,
                               na.rm = TRUE)) / 4)
  is_early_b <- !is.na(dt$rna_pattern) & dt$rna_pattern %in% early_pattern_names()
  cross_b <- as.integer(feat$cross_layer_concord)
  sl_b <- feat$score_layer
  ss_b <- feat$score_serum
  sr_b <- feat$score_rescue
  leak_b <- feat$leakage_mult

  rho_se <- 0.05
  i2_se <- 8
  for (iter in seq_len(n_mc)) {
    rho_p <- ifelse(is.na(rho_b), NA_real_,
                    tanh(atanh(pmin(0.999, pmax(-0.999, rho_b))) +
                           stats::rnorm(n_genes, 0, rho_se)))
    rho_p <- pmax(-1, pmin(1, rho_p))
    k <- pmax(0, pmin(4, round(.audit_na0(agree_b) * 4)))
    agree_p <- stats::rbeta(n_genes, k + 1, 4 - k + 1)
    i2_p <- ifelse(is.na(i2_b), NA_real_,
                   pmax(0, pmin(100, i2_b +
                                  stats::rnorm(n_genes, 0, i2_se))))
    score_early <- ifelse(is_early_b & !is.na(rho_p),
                          pmax(0, pmin(1, rho_p * .audit_na0(lrt_b))), 0)
    score_direction <- 0.5 * agree_p + 0.5 * cross_b
    pos <- 0.20 * sl_b + 0.20 * score_direction +
      0.20 * score_early + 0.10 * ss_b + 0.10 * sr_b
    het <- data.table::fcase(is.na(i2_p), 1.00,
                             i2_p < 50, 1.00,
                             i2_p < 70, 1.00,
                             i2_p < 90, 0.70,
                             default = 0.30)
    raw <- pos * leak_b * het
    score_mc[, iter] <- raw / max(raw, na.rm = TRUE)
    rank_mc[, iter] <- data.table::frank(-score_mc[, iter],
                                         na.last = "keep")
  }

  out <- data.table::data.table(gene_symbol = dt$gene_symbol)
  out[, audit_score_median := apply(score_mc, 1, stats::median,
                                    na.rm = TRUE)]
  out[, audit_score_lo95 := apply(score_mc, 1, stats::quantile,
                                  probs = 0.025, na.rm = TRUE)]
  out[, audit_score_hi95 := apply(score_mc, 1, stats::quantile,
                                  probs = 0.975, na.rm = TRUE)]
  out[, uncertainty_width := audit_score_hi95 - audit_score_lo95]
  out[, rank_median := apply(rank_mc, 1, stats::median, na.rm = TRUE)]
  out[, rank_lo95 := apply(rank_mc, 1, stats::quantile, probs = 0.025,
                           na.rm = TRUE)]
  out[, rank_hi95 := apply(rank_mc, 1, stats::quantile, probs = 0.975,
                           na.rm = TRUE)]
  out[, confidence_class := data.table::fcase(
    feat$is_hk | (audit_score_hi95 == 0 & audit_score_lo95 == 0),
    "excluded",
    audit_score_lo95 >= 0.5, "stable_high",
    audit_score_hi95 >= 0.5 & audit_score_lo95 < 0.5, "high_uncertain",
    audit_score_lo95 >= 0.3, "medium",
    default = "low"
  )]
  out
}

.audit_anchor_table <- function() {
  e <- new.env()
  ok <- tryCatch({
    data("pdactrace_external_anchors", package = "pdactrace", envir = e)
    TRUE
  }, error = function(e) FALSE)
  if (ok && exists("pdactrace_external_anchors", envir = e)) {
    return(data.table::as.data.table(e$pdactrace_external_anchors))
  }
  csv <- file.path("data-raw", "external_positives_anchors_v2.csv")
  if (file.exists(csv)) return(data.table::fread(csv))
  stop("External anchor table is unavailable.")
}

.audit_direction_of <- function(pattern) {
  out <- rep(NA_character_, length(pattern))
  out[pattern %in% c("Early_Burst_Up", "Early_Peak")] <- "UP"
  out[pattern %in% c("Early_Loss_Down", "Early_Trough")] <- "DOWN"
  out
}

.audit_na0 <- function(x) {
  x[is.na(x)] <- 0
  x
}

.audit_existing_flag <- function(dt, col) {
  if (col %in% names(dt)) return(!is.na(dt[[col]]) & dt[[col]] == TRUE)
  rep(FALSE, nrow(dt))
}

.audit_call <- function(anchor_tier, is_hk, is_plasma_hi, rescue_eligible) {
  out <- rep("novel_nomination", length(anchor_tier))
  out[anchor_tier != "none"] <- "confirmed_anchor"
  out[isTRUE_vec(is_hk) | isTRUE_vec(is_plasma_hi)] <- "rejected_artifact"
  rescued <- isTRUE_vec(rescue_eligible) &
    !(isTRUE_vec(is_hk) | isTRUE_vec(is_plasma_hi)) &
    anchor_tier == "none"
  out[rescued] <- "rescued_candidate"
  out
}

.audit_leakage_label <- function(feat) {
  if (isTRUE(feat$is_hk)) return("housekeeping hard gate")
  if (isTRUE(feat$is_plasma_hi)) return("plasma high-abundance half gate")
  "no leakage flag"
}

.audit_heterogeneity_label <- function(i2, mult) {
  if (is.na(i2)) return("I2 unavailable")
  sprintf("max I2 %.1f, multiplier %.2f", i2, mult)
}
