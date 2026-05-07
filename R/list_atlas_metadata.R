#' Inspect atlas metadata
#'
#' Returns the bundled `atlas_metadata` list with version, snapshot
#' date, scope description, cohort tables, and source CSV provenance.
#'
#' @return Named list. See `?atlas_metadata` for structure.
#' @examples
#'   meta <- list_atlas_metadata()
#'   meta$version
#'   meta$source_csvs
#' @export
list_atlas_metadata <- function() {
  e <- new.env()
  data("atlas_metadata", package = "pdactrace", envir = e)
  e$atlas_metadata
}
