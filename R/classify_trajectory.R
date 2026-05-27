#' Classify each gene's trajectory against the 12-template catalog
#'
#' Takes the output of [fit_stage_de()] and matches each LRT-significant
#' gene's z-scored 4-point profile (beta_N, beta_E, beta_M, beta_L) against
#' the **canonical 12-template catalog** (Early × 4 + Mid × 4 + Late × 2 +
#' Monotonic × 2; v0.4.0). The best-match template is assigned if
#' rho >= `rho_cutoff`. The atlas surface (`rna_pattern` in
#' `pdactrace_reference`) restricts visible calls to **Early × 4 only** —
#' Mid / Late / Monotonic best-matches are recorded internally but flagged
#' via `excluded_*_pattern` columns rather than surfaced. Including the 8
#' non-Early templates in the *matching* step is a deliberate design choice
#' for honest competition: a gene must beat 11 alternatives before it is
#' surfaced as Early.
#'
#' @param fit A `data.table` returned by [fit_stage_de()].
#' @param rho_cutoff Numeric. Minimum Pearson rho for assignment.
#'   Default 0.85 (matches phase33 canonical).
#' @param sig_only Logical. If `TRUE` (default), restrict to
#'   `lrt_significant == TRUE` rows.
#' @return A `data.table` with the input plus columns:
#'   * `rna_pattern` - character: best 12-template match or NA. Note:
#'     this returns the *raw* best-match label (any of the 12 templates).
#'     The atlas-surfaced version in `pdactrace_reference$rna_pattern`
#'     is filtered to Early × 4 only; non-Early calls are flagged via
#'     `excluded_mid_pattern` / `excluded_late_pattern` /
#'     `excluded_monotonic_pattern`.
#'   * `rna_pattern_rho` - numeric: Pearson rho with best template
#'   * `rna_pattern_rho_runner_up` - second-best rho (gap detection)
#' @examples
#' fit <- data.table::data.table(
#'   gene_symbol = c("G1", "G2"),
#'   beta_N = 0, beta_E = c(2, -2), beta_M = c(0.5, -0.5),
#'   beta_L = 0, lrt_padj = c(0.01, 0.02),
#'   lrt_significant = TRUE)
#' classify_trajectory(fit, sig_only = FALSE)
#'
#' \donttest{
#'   fit <- fit_stage_de(my_counts, my_stage, my_cohort)
#'   pat <- classify_trajectory(fit, rho_cutoff = 0.85)
#'   table(pat$rna_pattern)
#' }
#' @export
classify_trajectory <- function(fit, rho_cutoff = 0.85, sig_only = TRUE) {
  e <- new.env()
  data("default_templates", package = "pdactrace", envir = e)
  templates <- e$default_templates
  tpl_mat <- do.call(rbind, templates)

  rows <- if (sig_only) which(fit$lrt_significant) else seq_len(nrow(fit))
  if (length(rows) == 0L) {
    message("No LRT-significant genes; returning input unchanged.")
    return(fit)
  }

  pats   <- character(length(rows))
  rhos   <- numeric(length(rows))
  rhos2  <- numeric(length(rows))
  for (i in seq_along(rows)) {
    j <- rows[i]
    z <- c(fit$beta_N[j], fit$beta_E[j], fit$beta_M[j], fit$beta_L[j])
    if (any(is.na(z)) || stats::sd(z) == 0) {
      pats[i] <- NA_character_
      rhos[i] <- NA_real_
      rhos2[i] <- NA_real_
      next
    }
    z <- (z - mean(z)) / stats::sd(z)
    rho_v <- apply(tpl_mat, 1, function(t) stats::cor(z, t))
    ord <- order(rho_v, decreasing = TRUE)
    pats[i] <- if (rho_v[ord[1]] >= rho_cutoff)
      names(rho_v)[ord[1]] else NA_character_
    rhos[i]  <- as.numeric(rho_v[ord[1]])
    rhos2[i] <- as.numeric(rho_v[ord[2]])
  }

  out <- data.table::copy(fit)
  out[, rna_pattern := NA_character_]
  out[, rna_pattern_rho := NA_real_]
  out[, rna_pattern_rho_runner_up := NA_real_]
  out[rows, rna_pattern              := pats]
  out[rows, rna_pattern_rho          := rhos]
  out[rows, rna_pattern_rho_runner_up := rhos2]
  attr(out, "rho_cutoff")  <- rho_cutoff
  attr(out, "scope")       <- paste0(
    "12-template competitive matching (Ex4 + Mx4 + Lx2 + Monox2); ",
    "atlas surface restricts visible calls to Early x 4")
  attr(out, "n_templates") <- 12L
  out
}

#' Score one or more genes against all 12 templates
#'
#' Returns full rho vector (12 values) per gene. Useful for inspecting
#' how clearly a gene maps to a single pattern vs ambiguous/borderline,
#' and for diagnosing whether a non-Early best-match excluded the gene
#' from the Early × 4 atlas surface.
#'
#' @param pat A `data.table` from [classify_trajectory()] *or* a fit
#'   from [fit_stage_de()] (auto-classified if no `rna_pattern` col).
#' @param gene Character vector. One or more gene symbols.
#' @return `data.table` with `gene_symbol` and 12 rho columns
#'   (one per template): `rho_Early_Burst_Up`, `rho_Early_Loss_Down`,
#'   `rho_Early_Peak`, `rho_Early_Trough`, `rho_Mid_Plateau_Up`,
#'   `rho_Mid_Plateau_Down`, `rho_Mid_Peak`, `rho_Mid_Trough`,
#'   `rho_Late_Burst_Up`, `rho_Late_Loss_Down`, `rho_Monotonic_Up`,
#'   `rho_Monotonic_Down`.
#' @examples
#' fit <- data.table::data.table(
#'   gene_symbol = "G1", beta_N = 0, beta_E = 2,
#'   beta_M = 0.5, beta_L = 0, lrt_padj = 0.01,
#'   lrt_significant = TRUE)
#' pat <- classify_trajectory(fit, sig_only = FALSE)
#' score_trajectory(pat, "G1")
#' @export
score_trajectory <- function(pat, gene) {
  e <- new.env()
  data("default_templates", package = "pdactrace", envir = e)
  templates <- e$default_templates
  tpl_mat <- do.call(rbind, templates)

  hit <- pat[pat$gene_symbol %in% gene, ]
  if (nrow(hit) == 0L) {
    message("No matching genes in input.")
    return(data.table::data.table())
  }
  out <- data.table::data.table(gene_symbol = hit$gene_symbol)
  for (tpl in rownames(tpl_mat)) {
    out[, (paste0("rho_", tpl)) := NA_real_]
  }
  for (i in seq_len(nrow(hit))) {
    z <- c(hit$beta_N[i], hit$beta_E[i], hit$beta_M[i], hit$beta_L[i])
    if (any(is.na(z)) || stats::sd(z) == 0) next
    z <- (z - mean(z)) / stats::sd(z)
    rho_v <- apply(tpl_mat, 1, function(t) stats::cor(z, t))
    for (tpl in names(rho_v)) {
      data.table::set(out, i, paste0("rho_", tpl),
                        as.numeric(rho_v[tpl]))
    }
  }
  out
}
