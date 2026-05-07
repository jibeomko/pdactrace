#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════
# data-raw/build_reference.R
#
# Build pdactrace_reference.rda from 8 source phase CSV outputs.
#
# Base: phase33 LRT-significant genes (n ≈ 10,113).
# Scope: 4-Early-only by design (Mid genes flagged via
#        `excluded_mid_pattern = TRUE`, rna_pattern set to NA).
# T3 columns (direction-related): NA placeholder; populated by
# v0.2.0 future extension. Detectability columns are *also* NA at
# this build step — they are populated later by the Day-8 trained
# detectability model joining onto this base.
#
# Reproducibility: source-of-truth phase CSVs are read-only.
# Run from: package root.
#   Rscript data-raw/build_reference.R
# ═══════════════════════════════════════════════════════════════
suppressPackageStartupMessages({
  library(data.table)
})

# ── Paths ────────────────────────────────────────────────────
PROJ <- "/home/kjb9412/PDAC_biomarker"
TS   <- file.path(PROJ, "analysis/manuscript/tissue_to_serum_biomarker/results")
RNA  <- file.path(PROJ, "analysis/transcriptomics/results/figure1")
PKG  <- rprojroot::find_package_root_file()

# ── 1. RNA layer (phase33 = base) ────────────────────────────
cat("[1/8] phase33 RNA 12-template (atlas surface = 4-Early only) ...\n")
# Canonical 12-template filename (v0.4.0); falls back to legacy alias.
p33_path <- file.path(TS, "phase33_deseq2_coef_12template.csv")
if (!file.exists(p33_path)) {
  p33_path <- file.path(TS, "phase33_deseq2_coef_8template.csv")
}
p33 <- fread(p33_path)
# coef_pat is the canonical 12-template call. We surface only 4-Early in
# rna_pattern; the other 8 (Mid×4, Late×2, Monotonic×2) are flagged via
# excluded_*_pattern columns for transparency.
EARLY <- c("Early_Burst_Up", "Early_Loss_Down", "Early_Peak", "Early_Trough")
MID   <- c("Mid_Peak", "Mid_Trough", "Mid_Plateau_Up", "Mid_Plateau_Down")
LATE  <- c("Late_Burst_Up", "Late_Loss_Down")
MONO  <- c("Monotonic_Up", "Monotonic_Down")

ref <- p33[, .(
  gene_id                     = gene_id,
  gene_symbol                 = gene_symbol,
  rna_lrt_padj                = lrt_padj,
  rna_beta_N                  = beta_N,
  rna_beta_E                  = beta_E,
  rna_beta_M                  = beta_M,
  rna_beta_L                  = beta_L,
  # v0.1.1 detail: per-stage padj + lfcSE (for forest CI)
  rna_padj_E                  = padj_E,
  rna_padj_M                  = padj_M,
  rna_padj_L                  = padj_L,
  rna_lfcSE_E                 = lfcSE_E,
  rna_lfcSE_M                 = lfcSE_M,
  rna_lfcSE_L                 = lfcSE_L,
  rna_pattern_12              = coef_pat,                                 # internal — full 12-template label
  rna_pattern                 = fifelse(coef_pat %in% EARLY, coef_pat, NA_character_),
  rna_pattern_rho             = coef_rho,
  excluded_mid_pattern        = coef_pat %in% MID,
  excluded_late_pattern       = coef_pat %in% LATE,
  excluded_monotonic_pattern  = coef_pat %in% MONO
)]
# Drop rows without gene_symbol (cannot be queried by gene name)
ref <- ref[!is.na(gene_symbol) & gene_symbol != ""]
# Dedup on gene_symbol — keep min lrt_padj
setorder(ref, gene_symbol, rna_lrt_padj)
ref <- ref[, .SD[1], by = gene_symbol]
cat(sprintf("  base genes (LRT-sig with symbol): %d\n", nrow(ref)))

# Add ensembl_id / entrez_id placeholders (TBD — would require biomaRt
# in a follow-up build; left NA for v0.1.0 to keep build self-contained)
ref[, `:=`(ensembl_id = sub("\\..*$", "", gene_id),
            entrez_id  = NA_integer_)]

# NOTE (v0.4.0): per-contrast Wald padj (rna_wald_padj_E/M/L) are
# attached to the final atlas by data-raw/attach_wald_padj_columns.R
# *after* the meta + audit layers (build_meta_analysis.R,
# build_audit_scores_v2.R) have run. This base script intentionally
# does not include the Wald columns so a standalone re-run does not
# discard the downstream layered columns.

# Runner-up rho (second-best template) — phase33 doesn't expose it; we
# approximate by recomputing on z-scored 4-template profile.
load(file.path(PKG, "data", "default_templates.rda"))
tpl_mat <- do.call(rbind, default_templates)
compute_runner_up <- function(beta) {
  z <- c(0, beta[1], beta[2], beta[3])
  if (sd(z) == 0) return(NA_real_)
  z <- (z - mean(z)) / sd(z)
  rhos <- apply(tpl_mat, 1, function(t) cor(z, t))
  sort(rhos, decreasing = TRUE)[2]
}
ref[, rna_pattern_rho_runner_up := mapply(
  function(e, m, l) compute_runner_up(c(e, m, l)),
  rna_beta_E, rna_beta_M, rna_beta_L)]

# ── 2. Multi-cohort Stouffer consistency + per-cohort detail ─
cat("[2/8] multi_cohort_consistency Stouffer + per-cohort ...\n")
mcc_full <- fread(file.path(RNA, "multi_cohort_consistency.csv"))
# Aggregate trend / monotonic per cohort into list-columns
mcc_summary <- mcc_full[, .(
  gene_symbol,
  rna_stouffer_z       = stouffer_z,
  rna_stouffer_p       = stouffer_p,
  rna_stouffer_padj    = stouffer_padj,
  rna_cohort_agreement = n_cohorts_agree / n_cohorts_tested,
  # v0.1.1 detail: per-cohort sign vote + monotonic flag as list-cols
  rna_per_cohort_trend = lapply(seq_len(.N), function(i) {
    setNames(as.list(c(trend_TCGA[i], trend_CPTAC[i],
                        trend_GSE224564[i], trend_GSE79668[i])),
             c("TCGA", "CPTAC", "GSE224564", "GSE79668"))
  }),
  rna_per_cohort_monotonic = lapply(seq_len(.N), function(i) {
    setNames(as.list(as.logical(c(monotonic_TCGA[i], monotonic_CPTAC[i],
                                    monotonic_GSE224564[i], monotonic_GSE79668[i]))),
             c("TCGA", "CPTAC", "GSE224564", "GSE79668"))
  }))]
ref <- merge(ref, mcc_summary, by = "gene_symbol", all.x = TRUE)

# ── 3. Tissue protein layer (phase34) ────────────────────────
cat("[3/8] phase34 tissue protein 4-Early-pattern ...\n")
# Canonical 12-template filename (v0.4.0); falls back to legacy alias.
p34_path <- file.path(TS, "phase34_protein_pooled_12template.csv")
if (!file.exists(p34_path)) {
  p34_path <- file.path(TS, "phase34_protein_pooled_8template.csv")
}
p34 <- fread(p34_path)
p34 <- p34[, .(gene_symbol = gene,
                prot_pattern_8 = prot_pat,
                prot_rho       = prot_rho)]
p34[, prot_pattern := fifelse(prot_pattern_8 %in% EARLY,
                                prot_pattern_8, NA_character_)]
ref <- merge(ref, p34, by = "gene_symbol", all.x = TRUE)
ref[, rnaprot_concordant := !is.na(rna_pattern) & !is.na(prot_pattern) &
                              rna_pattern == prot_pattern]
# Tier1 gold = both 4-Early concordant
ref[, prot_tier := fcase(
  rnaprot_concordant == TRUE,           "Tier1_gold",
  !is.na(prot_pattern) & is.na(rna_pattern), "Tier2_silver_protOnly",
  !is.na(prot_pattern_8) & is.na(prot_pattern), "ProtOnly_Mid_excluded",
  default = NA_character_)]

# ── 4. scRNA cell origin (phase2c) ───────────────────────────
cat("[4/8] phase2c scRNA cell origin + tau ...\n")
p2c <- fread(file.path(TS, "phase2c_celltype_specificity.csv"))
ct_cols <- grep("^mean_", names(p2c), value = TRUE)
# top celltype + distribution list + specificity tau (v0.1.1 detail)
p2c_dt <- p2c[, .(gene_symbol = gene,
                    cell_origin_top = top_celltype,
                    cell_origin_padj = NA_real_,  # phase2c doesn't expose; left NA
                    cell_specificity_tau = tau,
                    cell_origin_distrib = lapply(seq_len(.N), function(i) {
                      v <- as.numeric(p2c[i, .SD, .SDcols = ct_cols])
                      names(v) <- sub("^mean_", "", ct_cols)
                      sort(v, decreasing = TRUE)
                    }))]
ref <- merge(ref, p2c_dt, by = "gene_symbol", all.x = TRUE)

# ── 5. Serum + translation class (phase42 + phase77) ─────────
cat("[5/8] phase42 + phase77 serum / translation_class + raw means ...\n")
p42 <- fread(file.path(TS, "phase42_pancreatitis_check.csv"))
# v0.1.1 detail: keep raw cohort means + t-test p-value
p42 <- p42[, .(gene_symbol = gene,
                serum_log2fc_PDAC_vs_HC = pdac_mean - hc_mean,
                serum_log2fc_Pan_vs_HC  = pan_mean - hc_mean,
                ann_pdac_mean = pdac_mean,
                ann_pan_mean  = pan_mean,
                ann_hc_mean   = hc_mean,
                ann_pdac_vs_pan_pval = t_pval)]

p77 <- fread(file.path(TS, "phase77_strict_RNAprotConvergent_serum.csv"))
# vs_serum schema: "Opposite" (inverse) or "Same" (concordant)
# translation_class assignment is direction-based (NOT padj-gated) per
# manuscript canonical phrasing; LTBP1 (vs_serum=Opposite, padj NA) is
# explicitly Class B exemplar.
p77 <- p77[, .(gene_symbol = gene,
                phase77_strict = TRUE,
                translation_class = fcase(
                  vs_serum == "Opposite", "B",
                  vs_serum == "Same",     "A",
                  default = NA_character_))]

ref <- merge(ref, p42, by = "gene_symbol", all.x = TRUE)
ref <- merge(ref, p77, by = "gene_symbol", all.x = TRUE)
ref[is.na(phase77_strict), phase77_strict := FALSE]
# Class C (decoupled): gene with tissue evidence + serum log2fc data,
# but not in phase77 strict (i.e., direction unstable / decoupled).
ref[is.na(translation_class) & !is.na(serum_log2fc_PDAC_vs_HC),
    translation_class := "C"]

# serum_detected: in phase42 OR phase77
ref[, serum_detected := !is.na(serum_log2fc_PDAC_vs_HC) | phase77_strict]
# n_cohorts_detected — best-effort from phase77 (3 serum cohorts) +
# phase42 (1). Without per-cohort detail in pdactrace build, use 0/1/2/3
# proxy: 1 if only phase42, 1 if only phase77, sum if both. v0.2.0 will
# replace with explicit per-cohort table.
ref[, serum_n_cohorts_detected := fifelse(
  !is.na(serum_log2fc_PDAC_vs_HC), 1L, 0L) +
  fifelse(phase77_strict, 1L, 0L)]

# ── 6. Resectable markers (phase29) ──────────────────────────
cat("[6/8] phase29 resectable markers ...\n")
p29 <- fread(file.path(TS, "phase29_resectable_markers.csv"))
ref[, resectable_marker := gene_symbol %in% p29$gene]
# Add phase29-native pattern (s_pat) — separate from phase33 coef_pat
# because cohort-adjustment can shift profiles. Manuscript canonical
# 76% Early enrichment uses phase29 s_pat, exposed here for transparency.
p29_pat <- p29[, .(gene_symbol = gene, resectable_pattern_phase29 = s_pat)]
ref <- merge(ref, p29_pat, by = "gene_symbol", all.x = TRUE)
n_res     <- ref[resectable_marker == TRUE, .N]
n_p29_e   <- ref[resectable_marker == TRUE &
                  resectable_pattern_phase29 %in% EARLY, .N]
n_p29_m   <- ref[resectable_marker == TRUE &
                  resectable_pattern_phase29 %in% MID, .N]
n_p33_e   <- ref[resectable_marker == TRUE & !is.na(rna_pattern), .N]
n_p33_m   <- ref[resectable_marker == TRUE & excluded_mid_pattern, .N]
cat(sprintf("  resectable: %d total\n", n_res))
cat(sprintf("    phase29 s_pat (manuscript canonical): %d Early / %d Mid\n",
            n_p29_e, n_p29_m))
cat(sprintf("    phase33 coef_pat (cohort-adjusted):  %d Early / %d Mid-excluded\n",
            n_p33_e, n_p33_m))

# ── 6b. Filter status + annotation (phase60_signalP_pipeline) ─
cat("[6b] phase60 filter status (7-step audit trail) ...\n")
p60 <- fread(file.path(TS, "phase60_signalP_pipeline.csv"))
p60_keep <- p60[, .(
  gene_symbol            = gene,
  flt_signal_peptide     = f_sp,
  flt_serum_measurable   = f_serum,
  flt_serum_significant  = f_pool_sig,
  flt_pancreatitis_pdac  = f_pan_pdac,
  flt_pancreatitis_hc    = f_pan_hc,
  flt_direction_match    = f_dir,
  flt_final              = f_final,
  ann_pool_logfc         = pool_logFC,
  ann_pool_padj          = pool_padj,
  # v0.1.1 detail: pancreatitis-vs-HC serum stats (filter underlying numbers)
  ann_pan_vs_hc_logfc    = pan_vs_hc_logFC,
  ann_pan_vs_hc_pval     = pan_vs_hc_pval,
  ann_pan_excluded_phase60 = pan_excluded)]
ref <- merge(ref, p60_keep, by = "gene_symbol", all.x = TRUE)
# Genes outside phase60 universe have NA filter status — fine, message in
# trace_filters() will say "not yet evaluated by 7-step pipeline".
n_p60 <- ref[!is.na(flt_signal_peptide), .N]
n_final <- ref[isTRUE(flt_final) | (!is.na(flt_final) & flt_final == TRUE), .N]
cat(sprintf("  phase60 coverage: %d genes / final pass: %d\n",
            n_p60, n_final))

# ── 7. Panel members (phase80) ───────────────────────────────
cat("[7/8] phase80 panel members ...\n")
p80 <- fread(file.path(TS, "phase80_ltbp1_pancreatitis_predeclared_panels.csv"))
panel_genes <- unique(unlist(strsplit(p80$genes, "[,;+ ]+")))
panel_genes <- panel_genes[panel_genes != "" & !is.na(panel_genes)]
ref[, panel_member := gene_symbol %in% panel_genes]
cat(sprintf("  panel_member genes: %d\n", ref[panel_member == TRUE, .N]))

# ── 8. Provenance + scope + atlas-level metadata ─────────────
cat("[8/8] provenance + T3 NA columns ...\n")
ref[, provenance := {
  parts <- character()
  if (!is.na(rna_lrt_padj))             parts <- c(parts, "phase33")
  if (!is.na(prot_pattern_8))           parts <- c(parts, "phase34")
  if (!is.na(rna_stouffer_z))           parts <- c(parts, "stouffer_consistency")
  if (!is.na(cell_origin_top))          parts <- c(parts, "phase2c")
  if (!is.na(serum_log2fc_PDAC_vs_HC))  parts <- c(parts, "phase42")
  if (!is.na(flt_signal_peptide))       parts <- c(parts, "phase60")
  if (phase77_strict)                    parts <- c(parts, "phase77")
  if (resectable_marker)                 parts <- c(parts, "phase29")
  if (panel_member)                      parts <- c(parts, "phase80")
  paste(parts, collapse = ",")
}, by = gene_symbol]
ref[, last_updated := as.Date("2026-05-04")]
ref[, evidence_scope := fcase(
  panel_member == TRUE,                               "panel_validated",
  serum_detected == TRUE | phase77_strict == TRUE,    "tissue_serum",
  default                                              = "tissue_only")]

# T3 / direction columns: NA placeholder (deferred to follow-up paper)
ref[, `:=`(
  serum_direction_label     = NA_character_,
  direction_model_trainable = NA,
  direction_model_card_ref  = NA_character_)]

# Drop internal helper columns not in schema
helper_cols <- intersect(c("rna_pattern_12", "rna_pattern_8",
                             "prot_pattern_12", "prot_pattern_8",
                             "gene_id", "prot_rho"), names(ref))
if (length(helper_cols) > 0) {
  ref[, (helper_cols) := NULL]
}

# Reorder columns to match schema_spec()
target_cols <- c(
  "gene_symbol", "ensembl_id", "entrez_id",
  "rna_lrt_padj", "rna_beta_N", "rna_beta_E", "rna_beta_M", "rna_beta_L",
  "rna_padj_E", "rna_padj_M", "rna_padj_L",
  # rna_wald_padj_E/M/L are injected post-hoc by attach_wald_padj_columns.R
  "rna_lfcSE_E", "rna_lfcSE_M", "rna_lfcSE_L",
  "rna_pattern", "rna_pattern_rho", "rna_pattern_rho_runner_up",
  "rna_stouffer_z", "rna_stouffer_p", "rna_stouffer_padj",
  "rna_cohort_agreement", "rna_per_cohort_trend",
  "rna_per_cohort_monotonic", "excluded_mid_pattern",
  "excluded_late_pattern", "excluded_monotonic_pattern",
  "prot_pattern", "prot_tier", "rnaprot_concordant",
  "cell_origin_top", "cell_origin_distrib", "cell_origin_padj",
  "cell_specificity_tau",
  "serum_detected", "serum_n_cohorts_detected",
  "serum_log2fc_PDAC_vs_HC", "serum_log2fc_Pan_vs_HC",
  "translation_class", "phase77_strict",
  "resectable_marker", "resectable_pattern_phase29", "panel_member",
  "flt_signal_peptide", "flt_serum_measurable", "flt_serum_significant",
  "flt_pancreatitis_pdac", "flt_pancreatitis_hc", "flt_direction_match",
  "flt_final",
  "ann_pool_logfc", "ann_pool_padj",
  "ann_pan_vs_hc_logfc", "ann_pan_vs_hc_pval", "ann_pan_excluded_phase60",
  "ann_pdac_mean", "ann_pan_mean", "ann_hc_mean",
  "ann_pdac_vs_pan_pval",
  "provenance", "last_updated", "evidence_scope",
  "serum_direction_label", "direction_model_trainable",
  "direction_model_card_ref")
miss_cols <- setdiff(target_cols, names(ref))
if (length(miss_cols) > 0) stop("Missing columns: ",
                                  paste(miss_cols, collapse = ", "))
extra_cols <- setdiff(names(ref), target_cols)
if (length(extra_cols) > 0) {
  cat("  WARN extra cols (will keep at end):",
      paste(extra_cols, collapse = ", "), "\n")
}
setcolorder(ref, c(target_cols, extra_cols))
setkey(ref, gene_symbol)

# ── Save ─────────────────────────────────────────────────────
pdactrace_reference <- ref
out <- file.path(PKG, "data", "pdactrace_reference.rda")
save(pdactrace_reference, file = out, compress = "xz")
cat(sprintf("\nSaved: %s\n", out))
cat(sprintf("Final rows: %d\n", nrow(pdactrace_reference)))
cat(sprintf("Final cols: %d\n", ncol(pdactrace_reference)))

# ── Summary checks (will be re-run as testthat separately) ────
cat("\n═══ Sanity checks ═══\n")
cat(sprintf("4-Early × tissue-prot concordant: %d\n",
            ref[rnaprot_concordant == TRUE &
                rna_pattern %in% EARLY, .N]))
cat(sprintf("Resectable markers: %d in atlas of 25 phase29 candidates\n",
            n_res))
cat(sprintf("  phase29 s_pat:   %d Early / %d Mid (manuscript canonical)\n",
            n_p29_e, n_p29_m))
cat(sprintf("  phase33 coef_pat: %d Early / %d Mid (cohort-adjusted)\n",
            n_p33_e, n_p33_m))
n_outside <- 25 - n_res
cat(sprintf("  Outside atlas (LRT-non-sig in phase33): %d\n", n_outside))
cat(sprintf("LTBP1 row:\n"))
print(ref[gene_symbol == "LTBP1", .(gene_symbol, rna_pattern,
                                     rna_pattern_rho, prot_pattern,
                                     translation_class, panel_member,
                                     evidence_scope, provenance)])
cat(sprintf("\nClass B genes:\n"))
print(ref[translation_class == "B", .(gene_symbol, rna_pattern,
                                        translation_class)])
