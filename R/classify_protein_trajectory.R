#' Classify each protein's trajectory against the 12-template catalog
#'
#' Protein-side wrapper of [classify_trajectory()]: takes the output of
#' [fit_stage_de_protein()], runs the same shared template-matching
#' engine (12 templates: Early x 4 + Mid x 4 + Late x 2 + Monotonic x 2),
#' and renames the output column `rna_pattern` to `prot_pattern` so
#' [assemble_user_evidence()] can consume both layers without colliding
#' on column names. As with the RNA side, the atlas surface restricts
#' visible `prot_pattern` calls to Early x 4; non-Early best-matches are
#' recorded internally but not surfaced.
#'
#' @param fit A `data.table` returned by [fit_stage_de_protein()].
#' @param rho_cutoff Numeric. Minimum Pearson rho for assignment.
#'   Default 0.85 (matches phase34 canonical).
#' @param sig_only Logical. If `TRUE` (default), restrict to
#'   `lrt_significant == TRUE` rows.
#' @param ... Forwarded to `classify_protein_trajectory()` when invoked
#'   through the `classify_prot_trajectory()` short alias.
#' @return A `data.table` with columns `gene_symbol`, `prot_pattern`,
#'   `prot_pattern_rho`, `prot_pattern_rho_runner_up`. The remaining
#'   columns of the input fit are preserved.
#' @examples
#' prot_fit <- data.table::data.table(
#'   gene_symbol = c("P1", "P2"),
#'   beta_N = 0, beta_E = c(2, -2), beta_M = c(0.5, -0.5),
#'   beta_L = 0, lrt_padj = c(0.01, 0.02),
#'   lrt_significant = TRUE)
#' classify_protein_trajectory(prot_fit, sig_only = FALSE)
#'
#' \donttest{
#'   prot_fit <- fit_stage_de_protein(my_intensity, my_stage, my_cohort)
#'   prot_pat <- classify_protein_trajectory(prot_fit)
#' }
#' @export
classify_protein_trajectory <- function(fit, rho_cutoff = 0.85,
                                         sig_only = TRUE) {
  out <- classify_trajectory(fit, rho_cutoff = rho_cutoff,
                              sig_only = sig_only)
  # If classify_trajectory returned the input unchanged (no LRT-sig
  # rows), the rna_pattern columns won't exist; skip the rename in
  # that case.
  data.table::setnames(out,
    old = c("rna_pattern", "rna_pattern_rho",
            "rna_pattern_rho_runner_up"),
    new = c("prot_pattern", "prot_pattern_rho",
            "prot_pattern_rho_runner_up"),
    skip_absent = TRUE)
  attr(out, "scope") <- paste0(
    "Protein 12-template competitive matching ",
    "(Ex4 + Mx4 + Lx2 + Monox2); atlas surface = Early x 4 only")
  attr(out, "n_templates") <- 12L
  out
}

#' @rdname classify_protein_trajectory
#' @export
classify_prot_trajectory <- function(...) classify_protein_trajectory(...)
