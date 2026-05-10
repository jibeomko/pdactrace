#' Descriptive feature-space similarity to bundled external anchors
#'
#' Computes how close each gene's evidence-feature vector is to the
#' centroid of a curated positive-anchor reference set (cosine or
#' inverse-Euclidean similarity in the z-scored feature space). This
#' is the **descriptive** layer that gives a continuous,
#' ML-flavoured prioritisation signal **without training a
#' supervised classifier** -- the bundled anchor labels are used
#' only to define the centroid, never as training labels for a
#' shipped predictor.
#'
#' This complements [evaluate_anchor_enrichment()]
#' (post-hoc top-N hypergeometric evaluation) and the deterministic
#' `audit_score`. The three are independent: `audit_score` ranks
#' by frozen rule, `anchor_similarity` ranks by feature-space
#' proximity to known biomarkers, and `evaluate_anchor_enrichment()`
#' is the evaluation harness that consumes either ranking via its
#' `score_col` argument.
#'
#' @section Method:
#' \enumerate{
#'   \item Build the per-gene feature matrix via
#'     [make_evidence_features()] with `scale = "z"` (mean 0, sd 1
#'     per feature over complete cases).
#'   \item Resolve the anchor set via the same internal helper
#'     [evaluate_anchor_enrichment()] uses, filtered by `tier`.
#'   \item Anchor centroid `c_k = mean over selected anchors of
#'     feature_k` (NAs skipped per feature).
#'   \item For each atlas gene g, compute `sim(g) = <z_g, c> /
#'     (||z_g|| * ||c||)` (cosine, default) or
#'     `sim(g) = 1 / (1 + sqrt(sum_k (z_g_k - c_k)^2))` (Euclidean).
#'   \item Genes with too few non-`NA` features (`< 3`) get
#'     `anchor_similarity = NA`.
#' }
#'
#' @param genes Optional character vector of gene symbols. `NULL`
#'   (default) scores every gene in the reference.
#' @param anchors Optional `data.table` with columns `gene` and
#'   `evidence_tier` (matching the bundled
#'   [pdactrace_external_anchors] schema). `NULL` uses the bundled
#'   anchor table.
#' @param method One of `"cosine"` (default) or `"euclidean"`.
#' @param tier Anchor tier filter: `"primary"` (T1_validated; the
#'   smallest, highest-evidence set; default), `"secondary"`
#'   (T1 + T2_literature_db), or `"all"` (also includes the single
#'   exploratory anchor). Matches the partition used by
#'   [evaluate_anchor_enrichment()].
#' @param reference Optional `data.table` to inject in place of the
#'   bundled atlas (for tests).
#' @return A `data.table` with one row per requested gene and
#'   columns `gene_symbol`, `anchor_similarity`, `anchor_n`
#'   (number of anchor rows the centroid was built from), and
#'   `anchor_tier`. Rows are sorted by `anchor_similarity`
#'   descending with `NA` last.
#' @examples
#' sim <- score_anchor_similarity(tier = "primary")
#' head(sim, 10)
#' @seealso [make_evidence_features()],
#'   [evaluate_anchor_enrichment()], [model_card()].
#' @export
score_anchor_similarity <- function(genes  = NULL,
                                     anchors = NULL,
                                     method = c("cosine", "euclidean"),
                                     tier   = c("primary",
                                                 "secondary",
                                                 "all"),
                                     reference = NULL) {
  method <- match.arg(method)
  tier   <- match.arg(tier)

  feats <- make_evidence_features(genes = genes,
                                   reference = reference,
                                   scale = "z")
  feat_cols <- setdiff(names(feats), "gene_symbol")
  X <- as.matrix(feats[, ..feat_cols])

  anc_dt <- .as_anchor_table(anchors, tier = tier)
  anc_genes <- anc_dt$gene
  anc_in <- anc_genes[anc_genes %in% feats$gene_symbol]
  if (length(anc_in) == 0L) {
    stop("No anchor genes are present in the requested feature ",
         "frame; cannot compute centroid.", call. = FALSE)
  }

  Xa <- X[match(anc_in, feats$gene_symbol), , drop = FALSE]
  centroid <- colMeans(Xa, na.rm = TRUE)
  if (any(!is.finite(centroid))) {
    centroid[!is.finite(centroid)] <- NA_real_
  }

  sim <- vapply(seq_len(nrow(X)), function(i) {
    .anchor_similarity_row(X[i, ], centroid, method = method)
  }, numeric(1L))

  out <- data.table::data.table(
    gene_symbol       = feats$gene_symbol,
    anchor_similarity = sim,
    anchor_n          = length(anc_in),
    anchor_tier       = tier)
  data.table::setorder(out, -anchor_similarity, na.last = TRUE)
  data.table::setattr(out, "method", method)
  out[]
}

# ---- internal -----------------------------------------------------------

.anchor_similarity_row <- function(z, c, method) {
  ok <- is.finite(z) & is.finite(c)
  if (sum(ok) < 3L) return(NA_real_)
  zz <- z[ok]; cc <- c[ok]
  if (method == "cosine") {
    nz <- sqrt(sum(zz^2)); nc <- sqrt(sum(cc^2))
    if (nz == 0 || nc == 0) return(NA_real_)
    return(sum(zz * cc) / (nz * nc))
  }
  # euclidean -> bounded similarity in (0, 1]
  d <- sqrt(sum((zz - cc)^2))
  1 / (1 + d)
}

# Same anchor-table accessor that audit_score.R uses, with the same
# tier semantics. We re-implement the tier filter here so this file
# can stand alone if R/audit_score.R is loaded later.
.as_anchor_table <- function(anchors = NULL,
                              tier = c("primary", "secondary", "all")) {
  tier <- match.arg(tier)
  dt <- if (!is.null(anchors)) {
    if (!data.table::is.data.table(anchors))
      data.table::as.data.table(anchors) else anchors
  } else {
    e <- new.env()
    ok <- tryCatch({
      utils::data("pdactrace_external_anchors", package = "pdactrace",
                   envir = e)
      TRUE
    }, error = function(err) FALSE)
    if (!ok) stop("`pdactrace_external_anchors` is unavailable.",
                  call. = FALSE)
    data.table::as.data.table(e$pdactrace_external_anchors)
  }
  if (tier == "primary") {
    flag <- as.logical(dt$include_primary_eval)
    dt[!is.na(flag) & flag]
  } else if (tier == "secondary") {
    flag <- as.logical(dt$include_secondary_eval)
    dt[!is.na(flag) & flag]
  } else {
    dt
  }
}
