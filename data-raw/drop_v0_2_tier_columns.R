#!/usr/bin/env Rscript
# ===========================================================
# data-raw/drop_v0_2_tier_columns.R
#
# v0.4.0 cleanup: drop the v0.2.0 confidence_tier system from
# the bundled pdactrace_reference. v0.3.0 audit_class supersedes
# the tier system; dual classification was creating reviewer
# confusion (e.g., LGALS3BP audit_class=high_confidence vs
# confidence_tier=ARTIFACT). See NEWS.md v0.4.0.
#
# Columns removed:
#   - confidence_tier      (5-label v0.2.0 GOLD/SILVER/.../ARTIFACT)
#   - early_onset_score    (composite tier score)
#   - heterogeneity_factor (only used by tier-system internals)
#
# Columns retained (audit-framework dependencies):
#   - max_abs_beta_meta    (used by audit_score)
#   - max_I2_meta          (used by heterogeneity_gate)
#   - meta_*               (per-contrast random-effects estimates)
# ===========================================================
suppressPackageStartupMessages({
  library(data.table)
})

PKG <- rprojroot::find_package_root_file()

cat("[1/3] Loading current atlas ...\n")
load(file.path(PKG, "data", "pdactrace_reference.rda"))
ref <- as.data.table(pdactrace_reference)
cat(sprintf("  before: %d rows x %d cols\n", nrow(ref), ncol(ref)))

drop <- c("confidence_tier", "early_onset_score", "heterogeneity_factor")
present <- intersect(drop, names(ref))
if (length(present) > 0) {
  ref[, (present) := NULL]
  cat(sprintf("  dropped: %s\n", paste(present, collapse = ", ")))
} else {
  cat("  (no v0.2.0 tier columns found; already removed)\n")
}

cat("[2/3] Verifying audit-framework dependencies preserved ...\n")
must_keep <- c("max_abs_beta_meta", "max_I2_meta", "audit_score",
                "audit_class")
missing <- setdiff(must_keep, names(ref))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}
cat("  OK: all audit-framework columns present\n")

cat("[3/3] Saving updated atlas ...\n")
pdactrace_reference <- ref
save(pdactrace_reference,
     file = file.path(PKG, "data", "pdactrace_reference.rda"),
     compress = "xz")
cat(sprintf("  saved: %d rows x %d cols\n", nrow(ref), ncol(ref)))
cat("Done.\n")
