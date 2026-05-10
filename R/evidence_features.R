#' Build a flat numeric feature matrix for downstream modelling
#'
#' Returns a wide numeric `data.table` (genes by features) suitable
#' for off-the-shelf supervised or descriptive modelling. Computes
#' the same arithmetic that [evidence_math()] exposes per-axis, but
#' as a single tabular shape so callers can hand it directly to
#' `glmnet::cv.glmnet()`, `ranger::ranger()`, scikit-learn, etc.
#'
#' Feature columns are deterministic functions of the bundled
#' [pdactrace_reference] and [pdactrace_protein_betas] objects; no
#' new fitting is performed. The full column list is grouped below
#' by evidence axis.
#'
#' @section Feature columns:
#'
#' \describe{
#'   \item{Trajectory fit}{`trajectory_rho_best`,
#'     `trajectory_delta_rho` (= `rho_best - rho_runner_up`).}
#'   \item{Effect magnitude}{`rna_beta_norm`
#'     (= sqrt(sum(rna_beta_E,M,L^2))), `rna_beta_max_abs`,
#'     `prot_beta_norm`, `prot_beta_max_abs`.}
#'   \item{Cohort robustness}{`cohort_agreement`,
#'     `stouffer_abs_z`, `stouffer_neglog10_fdr`,
#'     `max_meta_I2`.}
#'   \item{RNA-protein convergence}{`rna_protein_cosine`
#'     (cosine of (E,M,L) beta vectors),
#'     `rna_protein_concordant` (0/1),
#'     `prot_pattern_rho`.}
#'   \item{Serum bridge}{`serum_detected` (0/1),
#'     `serum_abs_log2fc`, `serum_n_cohorts_detected`,
#'     `phase77_strict` (0/1).}
#'   \item{Specificity}{`cell_specificity_tau`,
#'     `filter_pass_fraction` (passed / 7).}
#'   \item{Leakage flags}{`leakage_housekeeping` (0/1),
#'     `leakage_plasma_high_abundance` (0/1).}
#' }
#'
#' Plus one identifier column (`gene_symbol`). Logical inputs are
#' coerced to integer 0/1 so the resulting frame is uniformly
#' numeric except for the key.
#'
#' @section Scaling:
#' `scale = "z"` standardizes each feature column to mean 0, sd 1
#' over the **complete-case** subset (NAs are preserved as NA and
#' do not contribute to the moment estimates). `scale = "robust"`
#' uses median / MAD instead. `scale = "none"` (default) returns
#' raw values.
#'
#' @section NA policy:
#' Genes missing from [pdactrace_protein_betas] get `NA` for the
#' protein-side columns (`prot_beta_*`, `rna_protein_cosine`,
#' `prot_pattern_rho`). Genes with `NA` `rna_pattern` (no Early
#' best-match on the atlas surface) get `NA` for
#' `trajectory_rho_best` and `trajectory_delta_rho`. Set
#' `drop_na_rows = TRUE` to keep only complete-case rows; default
#' `FALSE` preserves the full atlas row count.
#'
#' @param genes Optional character vector of gene symbols. `NULL`
#'   (default) returns features for every gene in the reference.
#' @param reference Optional `data.table` to inject in place of the
#'   bundled atlas (used by tests and downstream pipelines).
#' @param protein_betas Optional `data.table` to inject in place of
#'   the bundled protein-betas table (used by tests).
#' @param scale One of `"none"` (default), `"z"`, or `"robust"`.
#' @param impute One of `"none"` (default), `"mean"`, or `"zero"`.
#'   Most genes are missing from `pdactrace_protein_betas` and / or
#'   from the bundled serum cohorts; the resulting `NA`s leave most
#'   genes incomplete-case. Use `"mean"` to fill `NA` with the
#'   column mean (over complete cases) -- standard pre-fit hygiene
#'   for `glmnet`. Use `"zero"` to fill `NA` with 0 -- conservative
#'   "no evidence" assumption, useful when the absence of an
#'   evidence layer is itself meaningful.
#' @param drop_na_rows Logical. If `TRUE`, drop rows containing
#'   any `NA` in feature columns. Default `FALSE`. Applied **after**
#'   `impute`, so `impute = "mean" + drop_na_rows = FALSE` is
#'   typically what you want for ML fitting.
#' @return A `data.table` keyed on `gene_symbol` with one row per
#'   gene and one column per feature.
#' @examples
#' feats <- make_evidence_features()
#' dim(feats)
#' head(feats)
#' z <- make_evidence_features(scale = "z", drop_na_rows = TRUE)
#' summary(z[, -1L])
#' @seealso [evidence_math()] for the per-axis nested-list view,
#'   [score_anchor_similarity()] for descriptive feature-space
#'   scoring against bundled anchors, [fit_user_evidence_model()]
#'   for the user-supplied supervised wrapper.
#' @export
make_evidence_features <- function(genes = NULL,
                                    reference = NULL,
                                    protein_betas = NULL,
                                    scale = c("none", "z", "robust"),
                                    impute = c("none", "mean", "zero"),
                                    drop_na_rows = FALSE) {
  scale  <- match.arg(scale)
  impute <- match.arg(impute)
  ref <- .get_reference(reference)
  pb  <- .ef_protein_betas(protein_betas)

  if (!is.null(genes)) {
    g <- genes
    ref <- ref[gene_symbol %in% g]
    pb  <- pb[gene_symbol %in% g]
  }

  # ---- merge protein betas onto the atlas keyed by gene_symbol ----
  pb_slim <- pb[, .(gene_symbol,
                     prot_beta_E_pb = prot_beta_E,
                     prot_beta_M_pb = prot_beta_M,
                     prot_beta_L_pb = prot_beta_L,
                     prot_pattern_rho_pb = prot_pattern_rho)]
  dt <- merge(ref, pb_slim, by = "gene_symbol",
               all.x = TRUE, sort = FALSE)

  # ---- vectorized per-axis features ------------------------------
  rna_E <- as.numeric(dt$rna_beta_E)
  rna_M <- as.numeric(dt$rna_beta_M)
  rna_L <- as.numeric(dt$rna_beta_L)
  rna_norm <- sqrt(rna_E^2 + rna_M^2 + rna_L^2)
  rna_maxabs <- pmax(abs(rna_E), abs(rna_M), abs(rna_L))

  prot_E <- as.numeric(dt$prot_beta_E_pb)
  prot_M <- as.numeric(dt$prot_beta_M_pb)
  prot_L <- as.numeric(dt$prot_beta_L_pb)
  prot_norm <- sqrt(prot_E^2 + prot_M^2 + prot_L^2)
  prot_maxabs <- pmax(abs(prot_E), abs(prot_M), abs(prot_L))

  # cosine of (E,M,L) RNA vs protein
  dot <- rna_E * prot_E + rna_M * prot_M + rna_L * prot_L
  cos_rp <- dot / (rna_norm * prot_norm)
  cos_rp[!is.finite(cos_rp)] <- NA_real_

  # Stouffer transforms
  st_z <- as.numeric(dt$rna_stouffer_z)
  st_padj <- as.numeric(dt$rna_stouffer_padj)
  st_neglog <- -log10(pmax(st_padj, .Machine$double.eps))

  # max meta I2 across the three contrasts
  max_i2 <- pmax(as.numeric(dt$meta_NvE_I2),
                  as.numeric(dt$meta_MvE_I2),
                  as.numeric(dt$meta_LvE_I2),
                  na.rm = TRUE)

  # filter pass fraction (7 boolean steps)
  flt_steps <- c("flt_signal_peptide", "flt_serum_measurable",
                  "flt_serum_significant", "flt_pancreatitis_pdac",
                  "flt_pancreatitis_hc", "flt_direction_match",
                  "flt_final")
  flt_mat <- vapply(flt_steps, function(nm) {
    v <- as.logical(dt[[nm]])
    ifelse(is.na(v), 0L, as.integer(v))
  }, integer(nrow(dt)))
  filter_pass_fraction <- rowSums(flt_mat) / length(flt_steps)

  # serum
  serum_log2fc <- as.numeric(dt$serum_log2fc_PDAC_vs_HC)

  out <- data.table::data.table(
    gene_symbol                  = dt$gene_symbol,
    trajectory_rho_best          = as.numeric(dt$rna_pattern_rho),
    trajectory_delta_rho         = as.numeric(dt$rna_pattern_rho) -
                                    as.numeric(dt$rna_pattern_rho_runner_up),
    rna_beta_norm                = rna_norm,
    rna_beta_max_abs             = rna_maxabs,
    prot_beta_norm               = prot_norm,
    prot_beta_max_abs            = prot_maxabs,
    cohort_agreement             = as.numeric(dt$rna_cohort_agreement),
    stouffer_abs_z               = abs(st_z),
    stouffer_neglog10_fdr        = st_neglog,
    max_meta_I2                  = max_i2,
    rna_protein_cosine           = cos_rp,
    rna_protein_concordant       = .ef_logical01(dt$rnaprot_concordant),
    prot_pattern_rho             = as.numeric(dt$prot_pattern_rho_pb),
    serum_detected               = .ef_logical01(dt$serum_detected),
    serum_abs_log2fc             = abs(serum_log2fc),
    serum_n_cohorts_detected     = as.integer(dt$serum_n_cohorts_detected),
    phase77_strict               = .ef_logical01(dt$phase77_strict),
    cell_specificity_tau         = as.numeric(dt$cell_specificity_tau),
    filter_pass_fraction         = filter_pass_fraction,
    leakage_housekeeping         = .ef_logical01(dt$audit_is_housekeeping),
    leakage_plasma_high_abundance =
      .ef_logical01(dt$audit_is_plasma_high_abundance))

  # Pmax with NAs may emit Inf/-Inf; canonicalize.
  num_cols <- setdiff(names(out), "gene_symbol")
  for (j in num_cols) {
    v <- out[[j]]
    v[!is.finite(v)] <- NA_real_
    data.table::set(out, j = j, value = v)
  }

  if (impute %in% c("mean", "zero")) {
    out <- .ef_impute(out, num_cols, type = impute)
  }
  if (isTRUE(drop_na_rows)) {
    keep <- stats::complete.cases(out[, num_cols, with = FALSE])
    out <- out[keep]
  }

  if (scale == "z") {
    out <- .ef_scale(out, num_cols, type = "z")
  } else if (scale == "robust") {
    out <- .ef_scale(out, num_cols, type = "robust")
  }
  data.table::setattr(out, "feature_set_version", "v1.0")
  data.table::setattr(out, "scale", scale)
  data.table::setattr(out, "impute", impute)
  out[]
}

# ---- internal helpers --------------------------------------------------

.ef_protein_betas <- function(pb = NULL) {
  if (!is.null(pb)) {
    return(if (data.table::is.data.table(pb)) pb
            else data.table::as.data.table(pb))
  }
  e <- new.env()
  utils::data("pdactrace_protein_betas", package = "pdactrace", envir = e)
  data.table::as.data.table(e$pdactrace_protein_betas)
}

.ef_logical01 <- function(x) {
  v <- suppressWarnings(as.logical(x))
  ifelse(is.na(v), NA_integer_, as.integer(v))
}

.ef_impute <- function(dt, cols, type = c("mean", "zero")) {
  type <- match.arg(type)
  for (j in cols) {
    v <- dt[[j]]
    if (!is.numeric(v)) next
    na_idx <- is.na(v)
    if (!any(na_idx)) next
    fill <- if (type == "zero") 0 else {
      cc <- v[!na_idx]
      if (length(cc) == 0L) next
      mean(cc)
    }
    v[na_idx] <- fill
    data.table::set(dt, j = j, value = v)
  }
  dt
}

.ef_scale <- function(dt, cols, type = c("z", "robust")) {
  type <- match.arg(type)
  for (j in cols) {
    v <- dt[[j]]
    if (!is.numeric(v)) next
    cc <- v[stats::complete.cases(v)]
    if (length(cc) < 2L) next
    if (type == "z") {
      m <- mean(cc); s <- stats::sd(cc)
      if (is.finite(s) && s > 0) {
        v <- (v - m) / s
      }
    } else {
      m <- stats::median(cc); s <- stats::mad(cc, constant = 1)
      if (is.finite(s) && s > 0) {
        v <- (v - m) / s
      }
    }
    data.table::set(dt, j = j, value = v)
  }
  dt
}
