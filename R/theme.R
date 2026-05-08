#' pdactrace NCS-grade ggplot2 theme
#'
#' Bundles the publication-grade ggplot2 theme used across the BiB
#' manuscript figures. Source-of-truth: `BiB_framework/scripts/bib_ncs_theme.R`.
#'
#' Two flavors:
#' * [pdactrace_panel_theme()] - `theme_void` base, for schematics
#' * [pdactrace_axes_theme()]  - `theme_minimal` base, for plots with axes
#'
#' Width constants (NCS spec, inches):
#' * `NCS_W_SINGLE = 3.46`  (88 mm, single column)
#' * `NCS_W_15COL  = 4.72`  (120 mm, 1.5-column)
#' * `NCS_W_DOUBLE = 7.08`  (180 mm, double column)
#' * `NCS_W_TRIPLE = 10.50` (full-width composite)
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point() +
#'   pdactrace_axes_theme()
#' identical(NCS_W_DOUBLE, 7.08)
#' @name pdactrace-theme
NULL

#' @rdname pdactrace-theme
#' @export
NCS_W_SINGLE <- 3.46
#' @rdname pdactrace-theme
#' @export
NCS_W_15COL  <- 4.72
#' @rdname pdactrace-theme
#' @export
NCS_W_DOUBLE <- 7.08
#' @rdname pdactrace-theme
#' @export
NCS_W_TRIPLE <- 10.50

.pdactrace_text_common <- function() {
  # Font family: R-standard "sans" alias maps to Helvetica on
  # Linux/macOS and Arial on Windows; both render bold-weight cleanly
  # under cairo_pdf and survive R CMD check's default PostScript
  # device. Hard-coding "Arial" breaks the PostScript device check
  # because Arial isn't in the default PostScript font database.
  ff <- "sans"
  ggplot2::theme(
    text = ggplot2::element_text(family = ff, color = "grey5"),
    plot.title = ggplot2::element_text(family = ff, size = 10,
                                          face = "bold", hjust = 0,
                                          color = "grey5",
                                          margin = ggplot2::margin(0, 0, 2, 0)),
    plot.subtitle = ggplot2::element_text(family = ff, size = 7.5,
                                              color = "grey25",
                                              hjust = 0, face = "italic",
                                              lineheight = 0.95,
                                              margin = ggplot2::margin(0, 0, 2, 0)),
    plot.tag = ggplot2::element_text(family = ff, size = 10,
                                       face = "bold", color = "grey5"),
    axis.title = ggplot2::element_text(family = ff, size = 8.5,
                                         face = "bold", color = "grey5"),
    axis.text  = ggplot2::element_text(family = ff, size = 7.5,
                                         face = "bold", color = "grey15"),
    axis.line  = ggplot2::element_line(linewidth = 0.4, color = "grey15"),
    axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "grey15"),
    axis.ticks.length = ggplot2::unit(2, "pt"),
    legend.text  = ggplot2::element_text(family = ff, size = 7.0,
                                           color = "grey10"),
    legend.title = ggplot2::element_text(family = ff, size = 7.5,
                                           face = "bold", color = "grey5"),
    legend.key.size = ggplot2::unit(8, "pt"),
    legend.margin   = ggplot2::margin(2, 3, 2, 3),
    legend.spacing.y = ggplot2::unit(1, "pt"),
    strip.text = ggplot2::element_text(family = ff, size = 7.5,
                                          face = "bold", color = "grey5",
                                          lineheight = 0.95,
                                          margin = ggplot2::margin(0, 0, 1.5, 0)),
    strip.background = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(linewidth = 0.12, color = "grey93"),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing = ggplot2::unit(0.18, "cm"),
    plot.margin = ggplot2::margin(4, 5, 4, 5),
    plot.background = ggplot2::element_rect(fill = "white", color = NA))
}

#' @rdname pdactrace-theme
#' @export
pdactrace_panel_theme <- function() {
  ggplot2::theme_void(base_size = 6) + .pdactrace_text_common() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 8, face = "bold",
                                            hjust = 0, color = "grey10",
                                            margin = ggplot2::margin(0, 0, 1, 0)))
}

#' @rdname pdactrace-theme
#' @export
pdactrace_axes_theme <- function() {
  ggplot2::theme_minimal(base_size = 6) + .pdactrace_text_common()
}

#' Palettes for pdactrace figures
#'
#' Named character vectors for ggplot2 `scale_*_manual()`.
#' * `pdactrace_pal_class`    - Class A / B / Other
#' * `pdactrace_pal_group`    - HC / Pancreatitis / PDAC
#' * `pdactrace_pal_dir`      - UP / DOWN / NS
#' * `pdactrace_pal_pattern`  - Early × 4 pattern colors (atlas-surfaced subset of the v0.4.0 12-template catalog)
#'
#' @examples
#' pdactrace_pal_class
#' pdactrace_pal_pattern
#' @name pdactrace-palettes
NULL

#' @rdname pdactrace-palettes
#' @export
pdactrace_pal_class <- c(
  "Class A (concordant)"      = "#C62828",
  "Class B (inverse stromal)" = "#0D47A1",
  "Class C (decoupled)"       = "grey60",
  "Other"                     = "grey75")

#' @rdname pdactrace-palettes
#' @export
pdactrace_pal_group <- c(
  HC           = "#90CAF9",
  Pancreatitis = "#FFB74D",
  Pan          = "#FFB74D",
  PDAC         = "#E57373")

#' @rdname pdactrace-palettes
#' @export
pdactrace_pal_dir <- c(UP = "#C62828", DOWN = "#1565C0", NS = "grey82")

#' @rdname pdactrace-palettes
#' @export
pdactrace_pal_pattern <- c(
  Early_Burst_Up   = "#C62828",
  Early_Loss_Down  = "#1565C0",
  Early_Peak       = "#F57C00",
  Early_Trough     = "#FFB74D")

#' Save a ggplot to PDF using cairo_pdf
#'
#' Convenience wrapper used in vignettes; logs filename + size.
#'
#' @param p A ggplot object.
#' @param dir Output directory.
#' @param name File stem (no extension).
#' @param w,h Width/height in inches.
#' @return Invisibly returns the saved file path.
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
#' tmp <- tempdir()
#' fp <- pdactrace_save(p, tmp, "demo", w = 3, h = 2.5)
#' file.exists(fp)
#' unlink(fp)
#' @export
pdactrace_save <- function(p, dir, name, w, h) {
  fp <- file.path(dir, sprintf("%s.pdf", name))
  ggplot2::ggsave(fp, p, width = w, height = h,
                    units = "in", device = grDevices::cairo_pdf)
  message(sprintf("  %-40s  (%g x %g in)", basename(fp), w, h))
  invisible(fp)
}
