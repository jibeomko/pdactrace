#' List the public datasets that contributed to the pdactrace atlas
#'
#' Returns the bundled `pdactrace_data_sources` table — one row per
#' public dataset (GEO / PRIDE / PDC / ICGC / GTEx / MassIVE / SRA)
#' that was used to populate the v0.99.x reference atlas. Designed so
#' a Bioconductor reviewer or end user can answer "where does this
#' data come from?" programmatically rather than by chasing prose in
#' the manuscript.
#'
#' @param layer Optional character vector to filter by data layer:
#'   one or more of `"RNA"`, `"Protein"`, `"scRNA"`, `"Serum"`,
#'   `"Pancreatitis"`, `"Validation"`. Default `NULL` returns the
#'   full table.
#' @return A `data.table` with columns:
#'   * `accession`   — public identifier (GEO ID, PXD, etc.)
#'   * `layer`       — RNA / Protein / scRNA / Serum / ...
#'   * `source_type` — repository name (GEO / PRIDE / PDC / ...)
#'   * `used_for`    — short role description
#'   * `citation`    — primary publication citation when available
#'   * `url`         — direct repository URL
#' @examples
#' head(list_data_sources())
#' list_data_sources(layer = "Serum")
#' nrow(list_data_sources())
#' @export
list_data_sources <- function(layer = NULL) {
  e <- new.env(parent = emptyenv())
  utils::data("pdactrace_data_sources", package = "pdactrace",
              envir = e)
  out <- e$pdactrace_data_sources
  if (!is.null(layer)) {
    if (!is.character(layer)) {
      stop("`layer` must be NULL or a character vector.",
           call. = FALSE)
    }
    keep <- layer
    out <- out[layer %in% keep]
  }
  out[]
}

#' Atlas provenance (version, build date, source repo, DOIs)
#'
#' Returns a named list summarising where the bundled atlas came
#' from: package version, atlas snapshot date, the manuscript
#' monorepo URL, both Zenodo concept DOIs (the package archive and
#' the manuscript-monorepo archive), and the cohort count. Designed
#' as a one-call provenance dossier for reviewers and citing authors.
#'
#' @return Named list with fields `package_version`, `atlas_version`,
#'   `build_date`, `n_cohorts`, `package_repo`, `manuscript_repo`,
#'   `package_doi`, `manuscript_doi`, and `data_layers`. Missing
#'   metadata fields are returned as `NA_character_`.
#' @examples
#' atlas_provenance()
#' atlas_provenance()$package_doi
#' @export
atlas_provenance <- function() {
  meta <- tryCatch(list_atlas_metadata(),
                   error = function(e) list())
  ds <- tryCatch(list_data_sources(),
                 error = function(e) NULL)
  layer_summary <- if (!is.null(ds)) {
    tab <- table(ds$layer)
    paste(sprintf("%s: %d", names(tab), as.integer(tab)),
          collapse = "; ")
  } else NA_character_

  pkg_v <- tryCatch(as.character(utils::packageVersion("pdactrace")),
                     error = function(e) NA_character_)

  list(
    package_version = pkg_v,
    atlas_version   = if (is.null(meta$version)) NA_character_
                      else meta$version,
    build_date      = if (is.null(meta$build_date)) NA_character_
                      else meta$build_date,
    n_cohorts       = if (is.null(meta$n_cohorts)) NA_integer_
                      else meta$n_cohorts,
    package_repo    = "https://github.com/jibeomko/pdactrace",
    manuscript_repo = "https://github.com/jibeomko/PDAC_biomarker",
    package_doi     = "10.5281/zenodo.20070208",
    manuscript_doi  = "10.5281/zenodo.20067849",
    data_layers     = layer_summary)
}
