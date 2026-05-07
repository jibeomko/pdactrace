#!/usr/bin/env Rscript
# Attach frozen v0.3.0 audit columns to data/pdactrace_reference.rda.
# (Legacy 7-feature attach kept for reproducibility; the canonical
# 3+2 framework attach is in data-raw/build_audit_scores_v2.R.)
# This script is intentionally mechanical: it joins deterministic
# audit scores and Monte Carlo summaries by gene_symbol without
# changing the scoring rule or regenerating analysis results.

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else
  file.path("data-raw", "attach_audit_columns.R")
pkg <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!file.exists(file.path(pkg, "DESCRIPTION"))) pkg <- normalizePath(getwd())

ref_path <- file.path(pkg, "data", "pdactrace_reference.rda")
det_path <- file.path(pkg, "data-raw", "audit_score_atlas_legacy_7feature.csv")
mc_path <- file.path(pkg, "data-raw", "audit_score_mc_v1.csv")
anchor_path <- file.path(pkg, "data-raw", "external_positives_anchors_v2.csv")

load(ref_path)
ref <- as.data.frame(pdactrace_reference, stringsAsFactors = FALSE)
det <- read.csv(det_path, stringsAsFactors = FALSE, check.names = FALSE)
mc <- read.csv(mc_path, stringsAsFactors = FALSE, check.names = FALSE)

copy_by_gene <- function(target, source, source_col, target_col) {
  idx <- match(target$gene_symbol, source$gene_symbol)
  target[[target_col]] <- source[[source_col]][idx]
  target
}

det_map <- c(
  score_layer = "audit_score_layer",
  score_direction = "audit_score_direction",
  score_early = "audit_score_early",
  score_serum = "audit_score_serum",
  score_rescue = "audit_score_rescue",
  positive_score = "audit_positive_score",
  leakage_mult = "audit_leakage_mult",
  het_mult = "audit_heterogeneity_mult",
  audit_score_raw = "audit_score_raw",
  audit_score = "audit_score",
  is_hk = "audit_is_housekeeping",
  is_plasma_hi = "audit_is_plasma_high_abundance",
  rescue_eligible = "audit_rescue_eligible"
)
for (src in names(det_map)) {
  ref <- copy_by_gene(ref, det, src, unname(det_map[[src]]))
}

mc_map <- c(
  audit_score_median = "audit_score_median",
  audit_score_lo95 = "audit_score_lo95",
  audit_score_hi95 = "audit_score_hi95",
  uncertainty_width = "audit_uncertainty_width",
  rank_median = "audit_rank_median",
  rank_lo95 = "audit_rank_lo95",
  rank_hi95 = "audit_rank_hi95",
  confidence_class = "audit_confidence_class"
)
for (src in names(mc_map)) {
  ref <- copy_by_gene(ref, mc, src, unname(mc_map[[src]]))
}

for (col in c("audit_is_housekeeping", "audit_is_plasma_high_abundance",
              "audit_rescue_eligible")) {
  ref[[col]] <- as.logical(ref[[col]])
}

pdactrace_reference <- ref
save(pdactrace_reference, file = ref_path, compress = "xz")

if (file.exists(anchor_path)) {
  pdactrace_external_anchors <- read.csv(anchor_path, stringsAsFactors = FALSE,
                                         check.names = FALSE)
  save(pdactrace_external_anchors,
       file = file.path(pkg, "data", "pdactrace_external_anchors.rda"),
       compress = "xz")
}

message("Attached audit columns to ", ref_path)
