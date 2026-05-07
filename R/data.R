#' Default trajectory templates (12-template catalog, v0.4.0)
#'
#' Built-in 12 hypothesis-driven trajectory templates used by
#' `classify_trajectory()` and `score_trajectory()`:
#' Early × 4 (Burst_Up, Loss_Down, Peak, Trough) + Mid × 4
#' (Plateau_Up, Plateau_Down, Peak, Trough) + Late × 2 (Burst_Up,
#' Loss_Down) + Monotonic × 2 (Up, Down). The atlas surface
#' (`pdactrace_reference$rna_pattern`) restricts visible best-match
#' calls to Early × 4 only — non-Early calls are flagged via
#' `excluded_*_pattern` columns. See `?early_pattern_names`,
#' `?mid_pattern_names_excluded` for membership lists.
#'
#' @format A named list of length 12, each element a numeric vector
#'   of length 4 (z-scored expected values across stages
#'   Normal / Early / Mid / Late).
#' @keywords datasets
"default_templates"

#' Atlas metadata
#'
#' Version, snapshot, source-script, and build metadata for the bundled
#' pdactrace reference atlas.
#'
#' @format A named list.
#' @keywords datasets
"atlas_metadata"

#' Random-effects meta-analysis table
#'
#' Cohort-level and random-effects meta-analysis summaries used to populate
#' v0.2.0 meta-analysis and tier columns in `pdactrace_reference`.
#'
#' @format A data.frame with one row per gene and meta-analysis summary
#'   columns.
#' @keywords datasets
"meta_analysis"
