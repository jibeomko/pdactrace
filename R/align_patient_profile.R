#' Align a single patient's omics profile against atlas stage trajectories
#'
#' Reports how a single patient's tumor-vs-matched-normal log2 fold-change
#' profile aligns with the frozen pdactrace atlas's stage-trajectory axes
#' (Early / Mid / Late), expressed as a Pearson correlation against each
#' cohort-adjusted stage-effect vector and a complementary audit-weighted
#' vote share. **The result is an alignment readout, not a stage prediction.**
#' No supervised model has been fit — the atlas is frozen and the function
#' is deterministic given (input, atlas).
#'
#' Two complementary signals are reported, never collapsed into one decision:
#' \describe{
#'   \item{Primary: `cor_to_stage_axis`}{Audit-weighted Pearson rho between
#'     the patient's per-gene log2FC vector and each of the three stage
#'     axes `rna_beta_E`, `rna_beta_M`, `rna_beta_L` (with Fisher-z 95% CI).
#'     Invariant to per-patient global shifts (library-size / age / FFPE
#'     drift). Has a closed-form null distribution.}
#'   \item{Secondary: `vote_share`}{Audit-weighted fraction of genes whose
#'     closest atlas β is each stage's β. Less rigorous but interpretable
#'     for non-statistical readers; surfaces disagreement-with-rho as a
#'     diagnostic.}
#' }
#'
#' Misalignment with all three stage axes is a legitimate output and does
#' not imply Normal — it indicates the patient does not resemble any of the
#' atlas's PDAC stage profiles. Limitations: rho is unstable below ~50
#' overlapping genes; the protein layer reports only categorical
#' concordance because the atlas does not carry per-stage numeric protein
#' effect sizes (use `layer_combine = "rna_only"` to skip the protein arm).
#'
#' @param rna_logfc Named numeric vector. Names are HGNC gene symbols;
#'   values are log2(patient_tumor / matched_normal). Missing values are
#'   dropped. Required.
#' @param prot_logfc Optional named numeric vector with the same format
#'   for the protein layer. When supplied, returns a categorical
#'   concordance table in `$prot` (the atlas does not carry per-stage
#'   numeric protein effect sizes, so a numeric alignment is not possible).
#' @param weight_by Per-gene weight in both rho and vote_share. One of
#'   `"audit_score"` (default — uses `audit_score` from the bundled atlas)
#'   or `"uniform"` (equal weights).
#' @param min_audit_score Minimum atlas `audit_score` to include in the
#'   gene dictionary. Default `0.3` matches the
#'   `low` -> `supported_uncertain` boundary in the audit framework.
#' @param top_n_genes Optional integer. Cap genes used at the top N by
#'   `gene_selection`. `NULL` = use all eligible genes.
#' @param gene_selection If `top_n_genes` is set, rank by either
#'   `"highest_audit"` (default) or `"highest_abs_beta"` (largest
#'   max-stage effect size).
#' @param min_genes Warn-only floor. If the number of overlapping genes is
#'   below this (default 50), a warning is issued. The function hard-errors
#'   below 10 (rho is meaningless on so few points).
#' @param layer_combine Either `"separate"` (default — RNA and protein
#'   reported separately, never combined into one β-space metric) or
#'   `"rna_only"` (skip the protein arm even if `prot_logfc` is supplied).
#' @param reference Optional atlas data.table override (advanced; for
#'   testability). Default `NULL` uses the bundled `pdactrace_reference`.
#' @param ... Forwarded to `align_patient_profile()` when invoked through
#'   the `align_patient()` short alias.
#'
#' @return A `list` of class `pdactrace_patient_alignment` with elements
#'   `$rna` (4-row data.table with stage / rho / CI / vote_share / weighted
#'   distance / n_genes_used), `$prot` (data.table or NULL), `$summary`
#'   (one-line character), and `$attrs` (named list of parameters and
#'   provenance metadata).
#'
#' @examples
#'   # Synthetic Late-stage patient profile
#'   set.seed(7)
#'   ref <- pdactrace:::.get_reference()
#'   pool <- ref[!is.na(rna_beta_L) & audit_score >= 0.3][seq_len(300)]
#'   patient <- setNames(
#'     pool$rna_beta_L + rnorm(nrow(pool), sd = 0.3),
#'     pool$gene_symbol)
#'   aln <- align_patient_profile(patient, top_n_genes = 200)
#'   print(aln)
#'   aln$rna
#' @export
align_patient_profile <- function(rna_logfc,
                                    prot_logfc       = NULL,
                                    weight_by        = c("audit_score",
                                                          "uniform"),
                                    min_audit_score  = 0.3,
                                    top_n_genes      = 500L,
                                    gene_selection   = c("highest_audit",
                                                           "highest_abs_beta"),
                                    min_genes        = 50L,
                                    layer_combine    = c("separate",
                                                           "rna_only"),
                                    reference        = NULL) {
  # ── 1. Argument validation ─────────────────────────────────
  weight_by      <- match.arg(weight_by)
  layer_combine  <- match.arg(layer_combine)
  gene_selection <- match.arg(gene_selection)
  if (!is.numeric(rna_logfc) || is.null(names(rna_logfc))) {
    stop("`rna_logfc` must be a named numeric vector ",
         "(names = HGNC gene symbols).")
  }
  rna_logfc <- rna_logfc[!is.na(rna_logfc)]

  # ── 2. Pull atlas (frozen reference) ───────────────────────
  ref <- .get_reference(reference)

  # ── 3. Restrict to eligible atlas genes ────────────────────
  pool <- ref[!is.na(rna_beta_E) & !is.na(rna_beta_M) &
                !is.na(rna_beta_L) & !is.na(audit_score) &
                audit_score >= min_audit_score,
                .(gene_symbol, rna_beta_E, rna_beta_M, rna_beta_L,
                  audit_score)]

  # ── 4. Intersect with patient input; record drop counts ───
  in_atlas  <- intersect(names(rna_logfc), pool$gene_symbol)
  n_dropped <- length(rna_logfc) - length(in_atlas)
  pat       <- rna_logfc[in_atlas]
  pool      <- pool[gene_symbol %in% in_atlas]
  data.table::setkey(pool, gene_symbol)
  pool      <- pool[in_atlas]   # align ordering

  # ── 5. Top-N gene selection ────────────────────────────────
  if (!is.null(top_n_genes) && nrow(pool) > top_n_genes) {
    rank_col <- if (gene_selection == "highest_audit") {
      -pool$audit_score
    } else {
      -pmax(abs(pool$rna_beta_E), abs(pool$rna_beta_M),
              abs(pool$rna_beta_L))
    }
    keep_ix <- order(rank_col, pool$gene_symbol)[seq_len(top_n_genes)]
    pool    <- pool[keep_ix]
    pat     <- pat[pool$gene_symbol]
  }

  # ── 6. Floor checks ────────────────────────────────────────
  n_used <- length(pat)
  if (n_used < 10L) {
    stop("Need >=10 atlas-overlapping genes for meaningful alignment; ",
         "got ", n_used, ".")
  }
  if (n_used < min_genes) {
    warning("Only ", n_used, " genes used (< min_genes=", min_genes,
            "); rho may be unstable.")
  }

  # ── 7. Audit-weighted Pearson rho per stage ────────────────
  w <- if (weight_by == "audit_score") pool$audit_score
       else rep(1, n_used)
  cor_E <- .weighted_cor(pat, pool$rna_beta_E, w)
  cor_M <- .weighted_cor(pat, pool$rna_beta_M, w)
  cor_L <- .weighted_cor(pat, pool$rna_beta_L, w)

  # ── 8. Audit-weighted vote share ───────────────────────────
  d <- cbind(N = abs(pat),  # rna_beta_N is identically 0
             E = abs(pat - pool$rna_beta_E),
             M = abs(pat - pool$rna_beta_M),
             L = abs(pat - pool$rna_beta_L))
  winner <- max.col(-d, ties.method = "first")
  votes  <- vapply(seq_len(4), function(k) sum(w[winner == k]) / sum(w),
                     numeric(1))
  weighted_dist <- vapply(seq_len(4),
                            function(k) sum(w * d[, k]) / sum(w),
                            numeric(1))

  # ── 9. Assemble RNA table ──────────────────────────────────
  rna_dt <- data.table::data.table(
    stage             = c("Normal", "Early", "Mid", "Late"),
    cor_to_stage_axis = c(NA_real_, cor_E$rho, cor_M$rho, cor_L$rho),
    cor_pval          = c(NA_real_, cor_E$p,   cor_M$p,   cor_L$p),
    cor_lo95          = c(NA_real_, cor_E$lo,  cor_M$lo,  cor_L$lo),
    cor_hi95          = c(NA_real_, cor_E$hi,  cor_M$hi,  cor_L$hi),
    vote_share        = votes,
    weighted_dist     = weighted_dist,
    n_genes_used      = n_used)

  # ── 10. Optional protein layer (categorical only) ─────────
  prot_dt <- NULL
  if (!is.null(prot_logfc) && layer_combine != "rna_only") {
    prot_dt <- .align_protein_categorical(prot_logfc, ref,
                                            min_audit_score)
  }

  # ── 11. Summary + structured return ────────────────────────
  summary_str <- .format_alignment_summary(rna_dt)
  out <- list(
    rna     = rna_dt,
    prot    = prot_dt,
    summary = summary_str,
    attrs   = list(
      weight_by                    = weight_by,
      min_audit_score              = min_audit_score,
      top_n_genes                  = top_n_genes,
      gene_selection               = gene_selection,
      layer_combine                = layer_combine,
      n_genes_input                = length(rna_logfc),
      n_genes_dropped_not_in_atlas = n_dropped,
      n_genes_used                 = n_used,
      atlas_version                = tryCatch(
        list_atlas_metadata()$version,
        error = function(e) NA_character_)
    )
  )
  class(out) <- c("pdactrace_patient_alignment", "list")
  out
}

#' @rdname align_patient_profile
#' @export
align_patient <- function(...) align_patient_profile(...)

#' @export
print.pdactrace_patient_alignment <- function(x, ...) {
  cat(x$summary, "\n\n")
  print(x$rna, row.names = FALSE)
  if (!is.null(x$prot)) {
    cat("\nProtein concordance (categorical):\n")
    print(x$prot, row.names = FALSE)
  }
  cat(sprintf("\n[atlas v%s; %d/%d input genes used]\n",
              x$attrs$atlas_version,
              x$attrs$n_genes_used,
              x$attrs$n_genes_input))
  invisible(x)
}

# ── Helpers ──────────────────────────────────────────────────

# Audit-weighted Pearson rho with Fisher-z 95% CI and two-sided p-value.
.weighted_cor <- function(x, y, w) {
  w   <- w / sum(w)
  mx  <- sum(w * x);            my <- sum(w * y)
  vx  <- sum(w * (x - mx) ^ 2); vy <- sum(w * (y - my) ^ 2)
  cov <- sum(w * (x - mx) * (y - my))
  rho <- cov / sqrt(vx * vy)
  rho <- max(-0.999999, min(0.999999, rho))   # clamp for atanh
  n   <- length(x)
  z   <- atanh(rho)
  se  <- 1 / sqrt(n - 3)
  ci_lo <- tanh(z - 1.96 * se)
  ci_hi <- tanh(z + 1.96 * se)
  p   <- 2 * stats::pnorm(-abs(z) * sqrt(n - 3))
  list(rho = rho, p = p, lo = ci_lo, hi = ci_hi)
}

# Categorical protein concordance: count how many of the patient's most-up
# / most-down genes match the atlas's surfaced prot_pattern label per
# Early-onset family. The atlas does not carry per-stage numeric protein
# effect sizes, so a numeric alignment in the same metric space as RNA is
# intentionally not implemented (would require fabricating data).
.align_protein_categorical <- function(prot_logfc, ref, min_audit_score) {
  prot_logfc <- prot_logfc[!is.na(prot_logfc)]
  early_pats <- c("Early_Burst_Up", "Early_Loss_Down",
                   "Early_Peak", "Early_Trough")
  pool <- ref[!is.na(prot_pattern) & prot_pattern %in% early_pats &
                !is.na(audit_score) & audit_score >= min_audit_score,
                .(gene_symbol, prot_pattern, audit_score)]
  in_atlas <- intersect(names(prot_logfc), pool$gene_symbol)
  if (length(in_atlas) == 0L) {
    return(data.table::data.table(
      prot_pattern = early_pats,
      n_genes      = 0L,
      n_concordant = 0L,
      concordance  = NA_real_))
  }
  pat  <- prot_logfc[in_atlas]
  pool <- pool[gene_symbol %in% in_atlas]
  data.table::setkey(pool, gene_symbol)
  pool <- pool[in_atlas]
  # Direction-up == patient logfc > 0; direction-down == < 0
  dir_up   <- pat > 0
  dir_down <- pat < 0
  out <- vapply(early_pats, function(p) {
    is_p <- pool$prot_pattern == p
    n    <- sum(is_p)
    if (n == 0) return(c(n = 0, n_conc = 0, conc = NA_real_))
    expect_up <- p %in% c("Early_Burst_Up", "Early_Peak")
    conc <- if (expect_up) sum(is_p & dir_up) else sum(is_p & dir_down)
    c(n = n, n_conc = conc, conc = conc / n)
  }, numeric(3))
  data.table::data.table(
    prot_pattern = early_pats,
    n_genes      = as.integer(out["n", ]),
    n_concordant = as.integer(out["n_conc", ]),
    concordance  = out["conc", ])
}

.format_alignment_summary <- function(rna_dt) {
  stage_idx <- which.max(rna_dt$cor_to_stage_axis[2:4]) + 1L
  best_stage <- rna_dt$stage[stage_idx]
  best_rho   <- rna_dt$cor_to_stage_axis[stage_idx]
  best_p     <- rna_dt$cor_pval[stage_idx]
  vote_idx   <- which.max(rna_dt$vote_share)
  vote_stage <- rna_dt$stage[vote_idx]
  vote_share <- rna_dt$vote_share[vote_idx]
  agree <- if (best_stage == vote_stage) "agrees" else "disagrees"
  sprintf(
    paste0("Best-match by rho: %s (rho=%.2f, p=%.2g, %d genes); ",
            "voting %s: %s (vote_share=%.2f)."),
    best_stage, best_rho, best_p, rna_dt$n_genes_used[1],
    agree, vote_stage, vote_share)
}
