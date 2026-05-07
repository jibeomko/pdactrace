#' Decompose a gene's audit score into its three axes and two gates
#'
#' Returns a transparent breakdown of how the frozen v0.3.0 audit rule
#' arrived at one gene's `audit_score` and `audit_class`. The point is
#' to make the score auditable in plain English: *"LTBP1 lands in
#' `supported_uncertain` because its `heterogeneity_gate` = 0.7 —
#' the cohort-level meta-analysis I2 exceeds 70% on at least one
#' contrast."*
#'
#' The audit rule is a frozen weighted sum:
#'
#' ```text
#' positive_score = 0.40 * evidence_strength
#'                + 0.35 * biological_coherence
#'                + 0.25 * translational_relevance
#' audit_score    = positive_score * leakage_gate * heterogeneity_gate
#' ```
#'
#' This function reports each component, the gate state, the audit
#' label that the gate transitions trigger, and a short prose
#' explanation derived from those numbers. It does NOT recompute any
#' new evidence — it reads the bundled atlas and presents what the
#' frozen rule already produced. Verbose dispatch is intended for
#' interactive console use; the underlying `data.table` is also
#' returned for programmatic use.
#'
#' @param gene_symbol HGNC gene symbol (length-1 character).
#' @param verbose Logical. If `TRUE` (default), prints a one-paragraph
#'   plain-language summary to the console. The structured breakdown
#'   is returned invisibly either way.
#' @return Invisibly returns a list with elements:
#'   * `gene` — the resolved HGNC symbol.
#'   * `audit_class`, `audit_score`, `positive_score` — top-line scalars.
#'   * `axes` — `data.table` (3 rows): name, weight, value, contribution.
#'   * `gates` — `data.table` (2 rows): name, value, triggered_by.
#'   * `explanation` — character(1), one-paragraph rationale.
#' @examples
#' explain_score("LGALS3BP")
#' explain_score("LTBP1")
#' explain_score("GAPDH")
#' @export
explain_score <- function(gene_symbol, verbose = TRUE) {
  if (!is.character(gene_symbol) || length(gene_symbol) != 1L) {
    stop("`gene_symbol` must be a length-1 character string.",
         call. = FALSE)
  }
  target <- gene_symbol
  ref <- .get_reference()
  row <- ref[gene_symbol == target]
  if (nrow(row) == 0L) {
    stop(sprintf("Gene '%s' is not in the bundled atlas.",
                 gene_symbol), call. = FALSE)
  }

  # Pull the five components and gates ------------------------------
  ev  <- as.numeric(row$audit_evidence_strength)
  bc  <- as.numeric(row$audit_biological_coherence)
  tr  <- as.numeric(row$audit_translational_relevance)
  lg  <- as.numeric(row$audit_leakage_gate)
  hg  <- as.numeric(row$audit_heterogeneity_gate)
  ps  <- as.numeric(row$audit_positive_score)
  sc  <- as.numeric(row$audit_score)
  cls <- as.character(row$audit_class)

  axes <- data.table::data.table(
    axis        = c("evidence_strength", "biological_coherence",
                    "translational_relevance"),
    weight      = c(0.40, 0.35, 0.25),
    value       = c(ev, bc, tr),
    contribution = c(0.40 * ev, 0.35 * bc, 0.25 * tr))

  # Identify what triggered each gate -------------------------------
  is_hk      <- isTRUE(row$audit_is_housekeeping)
  is_high_pl <- isTRUE(row$audit_is_plasma_high_abundance)
  leakage_trigger <- if (is_hk) "housekeeping flag (gate = 0)"
                     else if (is_high_pl)
                       "plasma-high-abundance flag (gate = 0.5)"
                     else "clean (gate = 1)"

  max_i2 <- as.numeric(row$max_I2_meta)
  het_trigger <- if (is.na(max_i2)) "no meta I2 available (gate = 1)"
                 else if (max_i2 < 70) sprintf(
                   "max meta I2 = %.0f%%, below 70%% (gate = 1)", max_i2)
                 else if (max_i2 < 90) sprintf(
                   "max meta I2 = %.0f%% in [70%%, 90%%) (gate = 0.7)",
                   max_i2)
                 else sprintf(
                   "max meta I2 = %.0f%%, >=90%% (gate = 0.3)", max_i2)

  gates <- data.table::data.table(
    gate       = c("leakage_gate", "heterogeneity_gate"),
    value      = c(lg, hg),
    triggered_by = c(leakage_trigger, het_trigger))

  # Compose plain-English explanation -------------------------------
  axis_text <- sprintf(
    "%s = %.2f * %.2f + %.2f * %.2f + %.2f * %.2f = %.3f",
    "positive_score", 0.40, ev, 0.35, bc, 0.25, tr, ps)
  gate_text <- sprintf(
    "audit_score = %.3f * leakage_gate(%.2f) * heterogeneity_gate(%.2f) = %.3f",
    ps, lg, hg, sc)

  reason <- switch(cls,
    high_confidence = sprintf(
      "Clean on both gates (leakage = %.2f, heterogeneity = %.2f) and a strong %.2f positive_score.",
      lg, hg, ps),
    supported_uncertain = sprintf(
      "Positive score is solid (%.2f) but the heterogeneity_gate = %.2f penalises cohort divergence (%s).",
      ps, hg, het_trigger),
    penalized = sprintf(
      "Positive score = %.2f, but the leakage_gate = %.2f (%s) discounts the score.",
      ps, lg, leakage_trigger),
    excluded = sprintf(
      "leakage_gate = 0 hard-zeros the score regardless of trajectory shape (%s).",
      leakage_trigger),
    low = sprintf(
      "Positive score = %.2f is below the high_confidence boundary; layered evidence is thin.",
      ps),
    sprintf("audit_class = '%s'.", cls))

  explanation <- paste(
    sprintf("%s lands in `%s` with audit_score = %.3f.",
            gene_symbol, cls, sc),
    axis_text, gate_text, reason, sep = "\n  ")

  if (isTRUE(verbose)) {
    cat(explanation, "\n")
  }
  out <- list(
    gene = gene_symbol,
    audit_class = cls,
    audit_score = sc,
    positive_score = ps,
    axes = axes,
    gates = gates,
    explanation = explanation)
  invisible(out)
}
