#' Compute weight-free Pareto layers on the three evidence axes
#'
#' `compute_pareto_layers()` implements layered non-dominated sorting on
#' `audit_evidence_strength`, `audit_biological_coherence`, and
#' `audit_translational_relevance`. The two reliability gates
#' (`audit_leakage_gate`, `audit_heterogeneity_gate`) are applied as a
#' hard pre-filter rather than as weights, mirroring the
#' "axes graded, gates reliability filter" structure of the frozen 3+2
#' audit framework. Ranking is layer-first, with crowding distance
#' (NSGA-II convention) as a within-layer tiebreaker and the frozen
#' deterministic `audit_score` as a final deterministic tiebreaker.
#'
#' This function is complementary to [compute_audit_score()] and does
#' not modify the frozen weighted formula. Use it when reviewers ask
#' whether the audit ranking is driven by the scalar weights
#' 0.40 / 0.35 / 0.25 — Pareto layers expose a weight-free alternative.
#'
#' @param atlas Optional user-supplied atlas (defaults to the bundled
#'   `pdactrace_reference`). When non-NULL it must contain the columns
#'   referenced by `axes`, `gate_cols`, and `score_col`.
#' @param axes Character(3) of the axis column names. Defaults to the
#'   three frozen axes.
#' @param gate_cols Character(2) of the gate column names. Defaults to
#'   `c("audit_leakage_gate", "audit_heterogeneity_gate")`.
#' @param score_col Column used both to restrict the candidate pool to
#'   the top `top_n` rows and as a deterministic tiebreaker after
#'   crowding distance. Defaults to `"audit_score"`.
#' @param gate_filter Logical; when `TRUE` (default), genes with any
#'   gate value `< 1` are excluded from the Pareto candidate pool and
#'   recorded with `pareto_excluded_by_gate = TRUE`.
#' @param top_n Integer; restrict the Pareto candidate pool to the top
#'   `top_n` rows by `score_col`. Defaults to `2000L` which comfortably
#'   covers the BiB headline top-100 anchor enrichment.
#' @return A data.table with one row per gene in `atlas` and the
#'   columns `pareto_layer` (integer; `NA` for excluded genes),
#'   `crowding_distance` (numeric; `Inf` for boundary points,
#'   `NA` for excluded), `pareto_rank` (integer 1..N within the
#'   candidate pool, `NA` for excluded), and `pareto_excluded_by_gate`
#'   (logical).
#' @examples
#' head(compute_pareto_layers(top_n = 200L))
#' @export
compute_pareto_layers <- function(atlas = NULL,
                                  axes = c("audit_evidence_strength",
                                            "audit_biological_coherence",
                                            "audit_translational_relevance"),
                                  gate_cols = c("audit_leakage_gate",
                                                 "audit_heterogeneity_gate"),
                                  score_col = "audit_score",
                                  gate_filter = TRUE,
                                  top_n = 2000L) {
  if (length(axes) != 3L || !is.character(axes)) {
    stop("`axes` must be a character vector of length 3.")
  }
  if (length(gate_cols) != 2L || !is.character(gate_cols)) {
    stop("`gate_cols` must be a character vector of length 2.")
  }
  top_n <- as.integer(top_n)
  if (length(top_n) != 1L || is.na(top_n) || top_n < 1L) {
    stop("`top_n` must be a single positive integer.")
  }

  ref <- .get_reference(atlas)
  required <- c("gene_symbol", axes, gate_cols, score_col)
  missing_cols <- setdiff(required, names(ref))
  if (length(missing_cols) > 0L) {
    stop("atlas is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  dt <- data.table::as.data.table(ref)
  out <- data.table::data.table(gene_symbol = dt$gene_symbol,
                                pareto_layer = NA_integer_,
                                crowding_distance = NA_real_,
                                pareto_rank = NA_integer_,
                                pareto_excluded_by_gate = TRUE)

  score_vals <- dt[[score_col]]
  pool_idx <- order(-score_vals, na.last = NA)
  pool_idx <- pool_idx[seq_len(min(top_n, length(pool_idx)))]
  if (length(pool_idx) == 0L) return(out)

  pool_dt <- dt[pool_idx]
  row_id <- pool_idx
  if (gate_filter) {
    gate_ok <- (pool_dt[[gate_cols[1L]]] >= 1) &
      (pool_dt[[gate_cols[2L]]] >= 1)
    gate_ok[is.na(gate_ok)] <- FALSE
    pool_dt <- pool_dt[gate_ok]
    row_id <- row_id[gate_ok]
  }
  if (nrow(pool_dt) == 0L) return(out)

  X <- as.matrix(pool_dt[, axes, with = FALSE])
  layers <- .fast_nondominated_sort_3d(X)
  crowd <- .crowding_distance_3d(X, layers)
  tiebreak <- pool_dt[[score_col]]

  # Order: layer ascending, crowding descending, audit_score descending.
  ord <- order(layers, -crowd, -tiebreak)
  ranks <- integer(length(ord))
  ranks[ord] <- seq_along(ord)

  data.table::set(out, i = row_id, j = "pareto_layer", value = layers)
  data.table::set(out, i = row_id, j = "crowding_distance", value = crowd)
  data.table::set(out, i = row_id, j = "pareto_rank", value = ranks)
  data.table::set(out, i = row_id, j = "pareto_excluded_by_gate",
                  value = FALSE)
  out
}

#' Evaluate Monte-Carlo Pareto stability under perturbed evidence
#'
#' `evaluate_pareto_stability()` recomputes the Pareto layer assignment
#' across `n_draws` Monte Carlo perturbations of the underlying
#' evidence (pattern-correlation rho, cohort agreement, and
#' between-cohort I^2), using the same perturbation kernel as
#' [propagate_uncertainty()]. For each gene it returns the proportion
#' of draws in which the gene appears in Pareto layer 1 of the
#' gate-passing candidate pool. This quantifies axis-agnostic ranking
#' robustness — a complement to the score-space CI returned by
#' [propagate_uncertainty()].
#'
#' This function is the MC analogue of [compute_pareto_layers()]; it
#' does not modify the frozen audit score in any draw.
#'
#' @param atlas Optional atlas (see [compute_pareto_layers()]).
#' @param n_draws Integer number of Monte Carlo draws. Defaults to
#'   `1000L`, matching the bundled MC summary.
#' @param axes Character(3) axis columns; see [compute_pareto_layers()].
#' @param gate_cols Character(2) gate columns.
#' @param score_col Column used to restrict the candidate pool to the
#'   top `top_n` rows by deterministic ranking, per draw.
#' @param top_n Integer pool size; see [compute_pareto_layers()].
#' @param seed Random seed.
#' @return A data.table with one row per gene in `atlas`, with columns
#'   `pareto_layer_median`, `pareto_stability_top1`,
#'   `pareto_layer_lo95`, `pareto_layer_hi95`,
#'   `pareto_top10_pct_stability`. Genes outside the gate-passing pool
#'   in every draw receive `NA` median/CI but a stability of `0`.
#' @examples
#' \donttest{
#'   evaluate_pareto_stability(n_draws = 50L)
#' }
#' @export
evaluate_pareto_stability <- function(atlas = NULL,
                                       n_draws = 1000L,
                                       axes = c("audit_evidence_strength",
                                                 "audit_biological_coherence",
                                                 "audit_translational_relevance"),
                                       gate_cols = c("audit_leakage_gate",
                                                      "audit_heterogeneity_gate"),
                                       score_col = "audit_score",
                                       top_n = 2000L,
                                       seed = 20250506L) {
  n_draws <- as.integer(n_draws)
  if (length(n_draws) != 1L || is.na(n_draws) || n_draws < 10L) {
    stop("`n_draws` must be a single integer >= 10.")
  }
  ref <- .get_reference(atlas)
  required <- c("gene_symbol", axes, gate_cols, score_col,
                # MC inputs (locked perturbation kernel from .audit_mc_table)
                "rna_pattern", "rna_pattern_rho", "rna_cohort_agreement",
                "rna_lrt_padj", "max_I2_meta")
  missing_cols <- setdiff(required, names(ref))
  if (length(missing_cols) > 0L) {
    stop("atlas is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  dt <- data.table::as.data.table(ref)
  n_genes <- nrow(dt)

  # Pre-compute the static feature components reused across draws.
  # These are deterministic per gene under the locked perturbation kernel.
  feat <- .audit_feature_table(dt)
  early_set <- early_pattern_names()
  rho_b <- dt$rna_pattern_rho
  agree_b <- dt$rna_cohort_agreement
  i2_b <- dt$max_I2_meta
  lrt_b <- pmin(1, -log10(pmax(dt$rna_lrt_padj, 1e-10, na.rm = TRUE)) / 4)
  is_early_b <- !is.na(dt$rna_pattern) & dt$rna_pattern %in% early_set
  cross_b <- as.integer(feat$cross_layer_concord)
  score_layer_b <- feat$score_layer
  score_serum_b <- feat$score_serum
  score_rescue_b <- feat$score_rescue
  leak_b <- feat$leakage_mult

  # Static axes that do NOT change with perturbation. Locked formula
  # mirrors compute_audit_score() lines 141 and 143.
  evidence_strength_static <- pmin(1, score_layer_b + 0.5 * score_rescue_b)
  translational_relevance_static <- score_serum_b

  withr::local_seed(seed)

  in_layer1 <- matrix(0L, nrow = n_genes, ncol = n_draws)
  layer_mat <- matrix(NA_integer_, nrow = n_genes, ncol = n_draws)
  pool_size <- min(as.integer(top_n), n_genes)
  rho_se <- 0.05
  i2_se <- 8

  for (iter in seq_len(n_draws)) {
    rho_p <- ifelse(is.na(rho_b), NA_real_,
                    tanh(atanh(pmin(0.999, pmax(-0.999, rho_b))) +
                           stats::rnorm(n_genes, 0, rho_se)))
    rho_p <- pmax(-1, pmin(1, rho_p))
    k <- pmax(0, pmin(4, round(.audit_na0(agree_b) * 4)))
    agree_p <- stats::rbeta(n_genes, k + 1, 4 - k + 1)
    i2_p <- ifelse(is.na(i2_b), NA_real_,
                   pmax(0, pmin(100, i2_b +
                                  stats::rnorm(n_genes, 0, i2_se))))
    score_early_p <- ifelse(is_early_b & !is.na(rho_p),
                            pmax(0, pmin(1, rho_p * .audit_na0(lrt_b))), 0)
    score_direction_p <- 0.5 * agree_p + 0.5 * cross_b
    biological_coherence_p <- (score_direction_p + score_early_p) / 2
    het_p <- data.table::fcase(is.na(i2_p), 1.00,
                                i2_p < 50, 1.00,
                                i2_p < 70, 1.00,
                                i2_p < 90, 0.70,
                                default = 0.30)

    # Pool selection per draw uses the perturbed score for consistency
    # with the deterministic top-N gate restriction.
    pos_p <- 0.40 * evidence_strength_static +
      0.35 * biological_coherence_p +
      0.25 * translational_relevance_static
    raw_p <- pos_p * leak_b * het_p
    mx_p <- max(raw_p, na.rm = TRUE)
    score_p <- if (is.finite(mx_p) && mx_p > 0) raw_p / mx_p else
      rep(0, n_genes)

    pool_idx <- order(-score_p, na.last = NA)
    pool_idx <- pool_idx[seq_len(min(pool_size, length(pool_idx)))]
    if (length(pool_idx) == 0L) next
    gate_ok <- (leak_b[pool_idx] >= 1) & (het_p[pool_idx] >= 1)
    gate_ok[is.na(gate_ok)] <- FALSE
    keep_idx <- pool_idx[gate_ok]
    if (length(keep_idx) == 0L) next

    X <- cbind(
      evidence_strength_static[keep_idx],
      biological_coherence_p[keep_idx],
      translational_relevance_static[keep_idx]
    )
    layers <- .fast_nondominated_sort_3d(X)
    layer_mat[keep_idx, iter] <- layers
    in_layer1[keep_idx, iter] <- as.integer(layers == 1L)
  }

  top10_pct_cut <- max(1L, as.integer(ceiling(pool_size * 0.1)))
  # "Top 10% by Pareto rank within draw" = the lowest layers whose
  # cumulative size covers the top10 cut. Rank ties go to "min" so
  # genes tied at the same Pareto layer share a rank.
  rank_mat <- apply(layer_mat, 2L, function(col) {
    keep <- !is.na(col)
    r <- rep(NA_integer_, length(col))
    if (any(keep)) {
      r[keep] <- data.table::frank(col[keep], ties.method = "min")
    }
    r
  })
  if (!is.matrix(rank_mat)) {
    rank_mat <- matrix(rank_mat, nrow = n_genes, ncol = n_draws)
  }
  top10_mat <- !is.na(rank_mat) & rank_mat <= top10_pct_cut

  layer_median <- apply(layer_mat, 1L, function(x) {
    if (all(is.na(x))) NA_integer_ else
      as.integer(stats::median(x, na.rm = TRUE))
  })
  layer_lo95 <- apply(layer_mat, 1L, function(x) {
    if (all(is.na(x))) NA_integer_ else
      as.integer(stats::quantile(x, 0.025, na.rm = TRUE))
  })
  layer_hi95 <- apply(layer_mat, 1L, function(x) {
    if (all(is.na(x))) NA_integer_ else
      as.integer(stats::quantile(x, 0.975, na.rm = TRUE))
  })
  stability_top1 <- rowSums(in_layer1) / n_draws
  stability_top10 <- rowSums(top10_mat) / n_draws

  data.table::data.table(
    gene_symbol = dt$gene_symbol,
    pareto_layer_median = layer_median,
    pareto_stability_top1 = stability_top1,
    pareto_layer_lo95 = layer_lo95,
    pareto_layer_hi95 = layer_hi95,
    pareto_top10_pct_stability = stability_top10
  )
}

# -- Internal helpers -------------------------------------------------

#' Fast non-dominated sort on 3-axis points
#'
#' Deb 2000 fast NDS, O(M N^2). Returns the integer layer assignment
#' for each row of `X` (a numeric matrix with three columns).
#'
#' Higher values are better in all three columns. Ties are allowed —
#' a gene that ties another on every axis does not dominate it.
#'
#' @param X Numeric matrix with three columns and N rows.
#' @return Integer vector of length N (layer indices, 1 = best).
#' @keywords internal
.fast_nondominated_sort_3d <- function(X) {
  if (!is.matrix(X) || ncol(X) != 3L) {
    stop(".fast_nondominated_sort_3d expects a 3-column matrix.")
  }
  N <- nrow(X)
  if (N == 0L) return(integer(0))
  if (N == 1L) return(1L)

  np <- integer(N)
  S <- vector("list", N)
  x1 <- X[, 1L]; x2 <- X[, 2L]; x3 <- X[, 3L]
  for (i in seq_len(N)) {
    geq <- (x1[i] >= x1) & (x2[i] >= x2) & (x3[i] >= x3)
    gt  <- (x1[i] >  x1) | (x2[i] >  x2) | (x3[i] >  x3)
    leq <- (x1[i] <= x1) & (x2[i] <= x2) & (x3[i] <= x3)
    lt  <- (x1[i] <  x1) | (x2[i] <  x2) | (x3[i] <  x3)
    geq[i] <- FALSE
    gt[i] <- FALSE
    leq[i] <- FALSE
    lt[i] <- FALSE
    S[[i]] <- which(geq & gt)   # i dominates these
    np[i] <- sum(leq & lt)       # this many dominate i
  }

  layers <- integer(N)
  current <- which(np == 0L)
  L <- 1L
  while (length(current) > 0L) {
    layers[current] <- L
    next_front <- integer(0)
    for (i in current) {
      js <- S[[i]]
      if (length(js) > 0L) {
        np[js] <- np[js] - 1L
        promoted <- js[np[js] == 0L]
        if (length(promoted) > 0L) next_front <- c(next_front, promoted)
      }
    }
    current <- unique(next_front)
    L <- L + 1L
  }
  layers
}

#' Crowding distance for 3-axis non-dominated layers
#'
#' Implements the NSGA-II crowding distance per layer: the sum of
#' axis-normalized neighbour gaps. Boundary points receive `Inf`.
#'
#' @param X Numeric matrix with three columns and N rows.
#' @param layers Integer vector of layer indices (output of
#'   [.fast_nondominated_sort_3d()]).
#' @return Numeric vector of length N.
#' @keywords internal
.crowding_distance_3d <- function(X, layers) {
  N <- nrow(X)
  if (N == 0L) return(numeric(0))
  d <- numeric(N)
  for (L in unique(layers)) {
    idx <- which(layers == L)
    if (length(idx) <= 2L) {
      d[idx] <- Inf
      next
    }
    for (axis in seq_len(3L)) {
      vals <- X[idx, axis]
      ord <- order(vals)
      sorted <- vals[ord]
      span <- sorted[length(sorted)] - sorted[1L]
      d[idx[ord[1L]]] <- Inf
      d[idx[ord[length(ord)]]] <- Inf
      if (span > 0 && length(idx) > 2L) {
        for (i in seq.int(2L, length(idx) - 1L)) {
          d[idx[ord[i]]] <- d[idx[ord[i]]] +
            (sorted[i + 1L] - sorted[i - 1L]) / span
        }
      }
    }
  }
  d
}
