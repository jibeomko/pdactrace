#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════
# Top 100 curated candidate table (Manuscript Table 2 prep)
#
# Audit-prioritized nomination table — NOT discovery claim.
# Anchor = NOT training; T1/T2/exploratory tier으로 confirm marker만 분류.
# pdactrace_call: confirmed_anchor / rescued_candidate /
#                  novel_nomination / rejected_artifact
# ═══════════════════════════════════════════════════════════════
suppressPackageStartupMessages({ library(data.table) })
PKG <- rprojroot::find_package_root_file()

# ── Load ──────────────────────────────────────────────────────
mc  <- fread(file.path(PKG, "data-raw", "audit_score_mc_v1.csv"))
load(file.path(PKG, "data", "pdactrace_reference.rda"))
ref <- as.data.table(pdactrace_reference)

anchors <- fread(file.path(PKG, "data-raw",
                  "external_positives_anchors_v2.csv"))
hk_list     <- fread(file.path(PKG, "data-raw",
                  "external_negatives_hpa_housekeeping.csv"))$gene_symbol
plasma_list <- setdiff(fread(file.path(PKG, "data-raw",
                  "external_negatives_hard_plasma.csv"))$gene_symbol,
                       c("APOA2","TTR"))

# Anchor tier lookup
anchor_tier_map <- setNames(anchors$evidence_tier, anchors$gene)
anchor_assay_map <- setNames(anchors$assay_context, anchors$gene)

# ── Join atlas detail to MC results ────────────────────────────
join_cols <- c("gene_symbol","rna_pattern","rna_pattern_rho",
                "rna_lrt_padj","rna_cohort_agreement",
                "prot_pattern","prot_tier",
                "cell_origin_top","cell_origin_distrib","cell_specificity_tau",
                "serum_detected","serum_log2fc_PDAC_vs_HC",
                "serum_log2fc_Pan_vs_HC","translation_class",
                "phase77_strict","resectable_marker","panel_member",
                "flt_signal_peptide","flt_final",
                "max_I2_meta","confidence_tier",
                # 3+2 framework outputs
                "audit_score","audit_class",
                "audit_evidence_strength","audit_biological_coherence",
                "audit_translational_relevance",
                "audit_leakage_gate","audit_heterogeneity_gate")
ref_join <- ref[, ..join_cols]
setkey(ref_join, gene_symbol)
setkey(mc, gene_symbol)
mc <- ref_join[mc]

# ── Top 100 by 3+2 audit_score (canonical ranking) ────────────
top100 <- mc[order(-audit_score)][1:100]

# ── Mechanism class heuristic (controlled vocab) ──────────────
ECM_KEYWORDS <- c("COL","FN","SPARC","TNC","VCAN","MFAP","BGN","LUM",
                   "POSTN","DCN","LOX","FBN","NID","LAMA","LAMB","LAMC",
                   "LRG","THBS","TGFB","COMP","MGP","MIA")
PROTEASE_KEYWORDS <- c("MMP","ADAM","ADAMTS","KLK","CTSB","CTSC","CTSD",
                        "CTSL","CTSS","CTSK","CTSZ","PLAU","PLAT",
                        "CAPN","CAPG","UCH")
MUCIN_KEYWORDS <- c("MUC","TFF","FCGBP","AGR","GP2","REG")
CEACAM_KEYWORDS <- c("CEACAM","DSG","CDH","ITGB","ITGA","CD55","CD44")
KERATIN_KEYWORDS <- c("KRT")
APO_KEYWORDS <- c("APO","SERPIN","ALB","TF","HP","HPX","CP","FGA","FGB",
                   "FGG","ORM","CRP","SAA")
IMMUNE_KEYWORDS <- c("IGH","IGK","IGL","C1Q","C3","C4","C5","C6","C7",
                      "C8","C9","CXCL","CCL","TNF","IL","IFN","CFB","CFH",
                      "CFI","FCN","MASP")

classify_mech <- function(gene, cell, ctd) {
  if (is.na(gene)) return("unclear")
  pre <- substr(gene, 1, 4)
  pre2 <- substr(gene, 1, 2)
  pre3 <- substr(gene, 1, 3)
  # cell-type origin first
  if (!is.na(cell)) {
    if (grepl("CAF|fibrobl|stell|stroma", cell, ignore.case = TRUE))
      return("ECM_stromal")
    if (grepl("ductal", cell, ignore.case = TRUE))
      return("ductal_epithelial")
    if (grepl("acinar", cell, ignore.case = TRUE))
      return("acinar_stress")
    if (grepl("T_cell|T cell|B_cell|macroph|neutroph|monocyte|DC|NK|mast",
               cell, ignore.case = TRUE))
      return("immune_inflammatory")
    if (grepl("endo", cell, ignore.case = TRUE))
      return("ECM_stromal")
  }
  # gene name keywords
  if (any(sapply(ECM_KEYWORDS, function(k) startsWith(gene, k))))
    return("ECM_stromal")
  if (any(sapply(PROTEASE_KEYWORDS, function(k) startsWith(gene, k))))
    return("protease_invasion")
  if (any(sapply(KERATIN_KEYWORDS, function(k) startsWith(gene, k))))
    return("ductal_epithelial")
  if (any(sapply(MUCIN_KEYWORDS, function(k) startsWith(gene, k))))
    return("ductal_epithelial")
  if (any(sapply(CEACAM_KEYWORDS, function(k) startsWith(gene, k))))
    return("cell_surface_adhesion")
  if (any(sapply(IMMUNE_KEYWORDS, function(k) startsWith(gene, k))))
    return("immune_inflammatory")
  if (any(sapply(APO_KEYWORDS, function(k) startsWith(gene, k))))
    return("metabolic_secreted")
  return("unclear")
}

top100[, mechanism_class := mapply(classify_mech, gene_symbol,
                                      cell_origin_top, cell_origin_distrib)]

# ── Anchor tier ────────────────────────────────────────────────
top100[, anchor_tier := anchor_tier_map[gene_symbol]]
top100[, anchor_tier := data.table::fifelse(is.na(anchor_tier), "none", anchor_tier)]

# ── pdactrace_call assignment ──────────────────────────────────
top100[, is_hk := gene_symbol %in% hk_list]
top100[, is_plasma_hi := gene_symbol %in% plasma_list]

top100[, pdactrace_call := data.table::fcase(
  audit_class == "excluded", "rejected_artifact",
  audit_class == "penalized", "rejected_artifact",
  anchor_tier %in% c("T1_validated","T2_literature_db","exploratory"),
    "confirmed_anchor",
  confidence_tier %in% c("ARTIFACT","OTHER") | is.na(confidence_tier),
    "rescued_candidate",
  default = "novel_nomination")]

# ── dominant_evidence_axis ─────────────────────────────────────
classify_axis <- function(rna_p, prot_p, cell, serum, sig_pep) {
  cell_str <- ifelse(is.na(cell), "", cell)
  bits <- c()
  if (!is.na(rna_p) && !is.na(prot_p)) bits <- c(bits, "RNA+Protein")
  if (grepl("CAF|stroma", cell_str, ignore.case = TRUE)) bits <- c(bits, "stromal")
  if (grepl("ductal", cell_str, ignore.case = TRUE)) bits <- c(bits, "ductal")
  if (grepl("acinar", cell_str, ignore.case = TRUE)) bits <- c(bits, "acinar")
  if (!is.na(serum) && serum) bits <- c(bits, "serum-detected")
  if (!is.na(sig_pep) && sig_pep) bits <- c(bits, "secreted")
  if (length(bits) == 0) return("RNA-only")
  paste(bits, collapse = "+")
}

top100[, dominant_evidence_axis := mapply(classify_axis,
        rna_pattern, prot_pattern, cell_origin_top,
        serum_detected, flt_signal_peptide)]

# ── artifact_flags ─────────────────────────────────────────────
flag_artifacts <- function(max_i2, rna_agree, conf_class, rna_pat,
                            prot_pat, serum_det) {
  flags <- c()
  if (!is.na(max_i2) && max_i2 >= 70) flags <- c(flags, "cohort_divergent")
  if (!is.na(rna_agree) && rna_agree < 0.5) flags <- c(flags, "single_cohort_RNA_signal")
  if (conf_class == "high_uncertain") flags <- c(flags, "score_bimodal_due_to_I2")
  if (!is.na(rna_pat) && !is.na(prot_pat)) {
    rna_dir <- ifelse(rna_pat %in% c("Early_Burst_Up","Early_Peak"), "UP", "DOWN")
    prot_dir <- ifelse(prot_pat %in% c("Early_Burst_Up","Early_Peak"), "UP", "DOWN")
    if (rna_dir != prot_dir) flags <- c(flags, "RNA_protein_direction_mismatch")
  }
  if (!is.na(serum_det) && !serum_det) flags <- c(flags, "serum_not_detected")
  if (length(flags) == 0) return("none")
  paste(flags, collapse = ";")
}

top100[, artifact_flags := mapply(flag_artifacts,
        max_I2_meta, rna_cohort_agreement, confidence_class,
        rna_pattern, prot_pattern, serum_detected)]

# ── scrna_support_celltypes simple ─────────────────────────────
top100[, scrna_support_celltypes := data.table::fifelse(
  is.na(cell_origin_top), "none", cell_origin_top)]

# ── serum_support ──────────────────────────────────────────────
top100[, serum_support := data.table::fcase(
  is.na(serum_detected), "not_measured",
  serum_detected & !is.na(translation_class), translation_class,
  serum_detected, "detected_no_class",
  default = "not_detected")]

# ── score_95ci formatted ───────────────────────────────────────
top100[, score_95ci := sprintf("[%.2f, %.2f]",
                                 audit_score_lo95, audit_score_hi95)]

# ── Final Table 2 schema ──────────────────────────────────────
table2_full <- top100[, .(
  rank = seq_len(.N),
  gene = gene_symbol,
  audit_score = round(audit_score, 3),
  audit_class,
  audit_score_median = round(audit_score_median, 3),
  score_95ci,
  confidence_class,
  anchor_tier,
  v0.2_tier = confidence_tier,
  pdactrace_call,
  rna_pattern,
  prot_pattern = prot_pattern,
  scrna_support_celltypes,
  serum_support,
  dominant_evidence_axis,
  mechanism_class,
  artifact_flags,
  note = ""
)]

# Add manual notes for headline genes
table2_full[gene == "LTBP1",
  note := "Headline rescue case — RNA-tier ARTIFACT but multi-layer convergent (myCAF Class B inverse, AUC 0.973 panel)"]
table2_full[gene == "SERPINA1",
  note := "External anchor + LTBP1+SERPINA1 panel partner; high_uncertain reflects cohort heterogeneity"]
table2_full[gene == "THBS2",
  note := "External anchor (Kim 2017 Sci Transl Med); stable_high with tight CI"]
table2_full[gene == "TIMP1",
  note := "External anchor (Capello 2017 + Cohen 2018 CancerSEEK)"]
table2_full[gene == "MUC16",
  note := "External anchor (Cohen 2018) with gene-level mapping caveat (CA125 epitope)"]

# ── Outputs ───────────────────────────────────────────────────
fwrite(table2_full,
        file.path(PKG, "data-raw", "top100_candidate_table_v0_3.csv"))

# Manuscript shortlist 30: high_confidence + supported_uncertain only,
# rejected_artifact 제외 (penalized/excluded는 자동 제외됨)
table2_short <- table2_full[
  pdactrace_call != "rejected_artifact" &
  audit_class %in% c("high_confidence","supported_uncertain")][1:30]
fwrite(table2_short,
        file.path(PKG, "data-raw", "table2_shortlist_30.csv"))

# Mechanism summary
mech_summary <- table2_full[, .(n = .N), by = .(mechanism_class, pdactrace_call)]
mech_summary <- dcast(mech_summary, mechanism_class ~ pdactrace_call,
                       value.var = "n", fill = 0)
fwrite(mech_summary,
        file.path(PKG, "data-raw", "top100_mechanism_summary.csv"))

# ── Reporting ──────────────────────────────────────────────────
cat("========================================\n")
cat("=== pdactrace_call x audit_class (top 100) ===\n")
cat("========================================\n")
print(table(table2_full$pdactrace_call, table2_full$audit_class))

cat("\n========================================\n")
cat("=== Mechanism class breakdown ===\n")
cat("========================================\n")
print(mech_summary)

cat("\n========================================\n")
cat("=== Anchor confirmations in top 100 ===\n")
cat("========================================\n")
print(table2_full[anchor_tier != "none",
        .(rank, gene, audit_score, score_95ci,
          anchor_tier, audit_class, pdactrace_call)])

cat("\n========================================\n")
cat("=== Headline rescued/novel candidates (top 30 shortlist) ===\n")
cat("========================================\n")
print(table2_short[, .(rank, gene, audit_score, score_95ci,
                          audit_class, pdactrace_call,
                          mechanism_class, dominant_evidence_axis)])

cat(sprintf("\nSaved 3 outputs:\n"))
cat(sprintf("  %s (n=100)\n",
    file.path("data-raw", "top100_candidate_table_v0_3.csv")))
cat(sprintf("  %s (n=30, manuscript)\n",
    file.path("data-raw", "table2_shortlist_30.csv")))
cat(sprintf("  %s\n",
    file.path("data-raw", "top100_mechanism_summary.csv")))
