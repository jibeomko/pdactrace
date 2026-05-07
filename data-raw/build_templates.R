#!/usr/bin/env Rscript
# ===========================================================
# Build default_templates.rda  (12-template catalog, v0.4.0)
#
# 12 pre-declared trajectory templates on z-score scale, organized
# by onset-stage family:
#   - Early × 4: Peak, Trough, Burst_Up, Loss_Down  (onset at Stage I)
#   - Mid × 4:   Peak, Trough, Plateau_Up, Plateau_Down (onset Stage II)
#   - Late × 2:  Burst_Up, Loss_Down (onset Stage III/IV; Peak/Trough/
#                Plateau are degenerate at the last stage because they
#                require post-onset return/sustained data points)
#   - Monotonic × 2: Up, Down (linear pan-stage progression, max at L)
#
# Manuscript surface policy (unchanged from v0.3.0):
#   `pdactrace_reference$rna_pattern` exposes only Early × 4. Mid /
#   Late / Monotonic calls are flagged via excluded_*_pattern columns
#   for transparency. The 12-template catalog widens the negative-
#   evidence pool so Early calls compete against 11 alternatives, not
#   7 — strengthening the manuscript's "76% Early-fall in 25 markers"
#   self-validation claim.
#
# Source-of-truth: tissue_to_serum_biomarker/scripts/ncs_theme.R
#   (canonical 4-point profiles for all 12 templates).
# ===========================================================

# Raw 4-point profiles - order: (Normal, Early, Mid, Late)
default_templates_raw <- list(
  # Early × 4 (onset at Stage I)
  Early_Peak       = c(-1.5,  1.5,  0.5, -0.5),
  Early_Trough     = c( 1.5, -1.5, -0.5,  0.5),
  Early_Burst_Up   = c(-1.5,  0.5,  0.5,  0.5),
  Early_Loss_Down  = c( 1.5, -0.5, -0.5, -0.5),
  # Mid × 4 (onset at Stage II)
  Mid_Peak         = c(-0.5,  0.5,  1.5,  0.5),
  Mid_Trough       = c( 0.5, -0.5, -1.5, -0.5),
  Mid_Plateau_Up   = c(-1.0, -0.5,  0.7,  0.8),
  Mid_Plateau_Down = c( 1.0,  0.5, -0.7, -0.8),
  # Late × 2 (onset at Stage III/IV)
  Late_Burst_Up    = c(-0.5, -0.5, -0.5,  1.5),
  Late_Loss_Down   = c( 0.5,  0.5,  0.5, -1.5),
  # Pan-stage monotonic × 2
  Monotonic_Up     = c(-1.5, -0.5,  0.5,  1.5),
  Monotonic_Down   = c( 1.5,  0.5, -0.5, -1.5)
)

# Z-normalize each template (mean 0, sd 1)
default_templates <- lapply(default_templates_raw,
                              function(v) (v - mean(v)) / sd(v))

# Metadata
attr(default_templates, "stage_order")    <- c("Normal", "Early", "Mid", "Late")
attr(default_templates, "onset_class")    <- c(
  Early_Peak       = "Early", Early_Trough     = "Early",
  Early_Burst_Up   = "Early", Early_Loss_Down  = "Early",
  Mid_Peak         = "Mid",   Mid_Trough       = "Mid",
  Mid_Plateau_Up   = "Mid",   Mid_Plateau_Down = "Mid",
  Late_Burst_Up    = "Late",  Late_Loss_Down   = "Late",
  Monotonic_Up     = "Pan",   Monotonic_Down   = "Pan")
attr(default_templates, "duration_class") <- c(
  Early_Peak       = "transient", Early_Trough     = "transient",
  Early_Burst_Up   = "sustained", Early_Loss_Down  = "sustained",
  Mid_Peak         = "transient", Mid_Trough       = "transient",
  Mid_Plateau_Up   = "sustained", Mid_Plateau_Down = "sustained",
  Late_Burst_Up    = "sustained", Late_Loss_Down   = "sustained",
  Monotonic_Up     = "sustained", Monotonic_Down   = "sustained")
attr(default_templates, "direction_class") <- c(
  Early_Peak       = "Up",   Early_Trough     = "Down",
  Early_Burst_Up   = "Up",   Early_Loss_Down  = "Down",
  Mid_Peak         = "Up",   Mid_Trough       = "Down",
  Mid_Plateau_Up   = "Up",   Mid_Plateau_Down = "Down",
  Late_Burst_Up    = "Up",   Late_Loss_Down   = "Down",
  Monotonic_Up     = "Up",   Monotonic_Down   = "Down")
attr(default_templates, "version") <- "v0.4.0"
attr(default_templates, "scope")   <- "12-template catalog (E×4 + M×4 + L×2 + Mono×2)"
attr(default_templates, "surfaced_in_atlas") <- c(
  "Early_Peak", "Early_Trough", "Early_Burst_Up", "Early_Loss_Down")
attr(default_templates, "excluded_by_design") <- c(
  "Mid_Peak", "Mid_Trough", "Mid_Plateau_Up", "Mid_Plateau_Down",
  "Late_Burst_Up", "Late_Loss_Down",
  "Monotonic_Up", "Monotonic_Down")
attr(default_templates, "source") <-
  "tissue_to_serum_biomarker/scripts/ncs_theme.R (canonical 12 templates)"

out <- file.path(rprojroot::find_package_root_file(), "data",
                  "default_templates.rda")
save(default_templates, file = out, compress = "xz")
cat(sprintf("Saved: %s\n", out))
cat(sprintf("12 templates: %s\n",
            paste(names(default_templates), collapse = ", ")))
cat("Atlas surface (rna_pattern): Early × 4 only\n")
cat("Excluded by design (flagged via excluded_*_pattern):\n")
cat("  Mid × 4, Late × 2, Monotonic × 2\n")
