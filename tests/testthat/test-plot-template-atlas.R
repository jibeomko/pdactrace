test_that("plot_template_atlas('rna') returns 12 named ggplot panels", {
  panels <- plot_template_atlas("rna")
  expect_length(panels, 12L)
  expect_true(all(c(early_pattern_names(),
                    mid_pattern_names_excluded(),
                    "Late_Burst_Up", "Late_Loss_Down",
                    "Monotonic_Up", "Monotonic_Down")
                  %in% names(panels)))
  for (p in panels) expect_s3_class(p, "ggplot")
})

test_that("plot_template_atlas('protein') returns the same shape", {
  panels <- plot_template_atlas("protein")
  expect_length(panels, 12L)
  for (p in panels) expect_s3_class(p, "ggplot")
})

test_that("every RNA template cohort is non-empty", {
  for (tpl in c(early_pattern_names(),
                 mid_pattern_names_excluded(),
                 "Late_Burst_Up", "Late_Loss_Down",
                 "Monotonic_Up", "Monotonic_Down")) {
    agg <- pdactrace:::.template_aggregate("rna", tpl)
    expect_gt(length(unique(agg$gene_symbol)), 0,
               label = paste0("RNA cohort for template ", tpl))
  }
})

test_that("templates argument subsetting works", {
  out <- plot_template_atlas("rna",
                              templates = c("Early_Burst_Up",
                                              "Late_Loss_Down"))
  expect_length(out, 2L)
  expect_named(out, c("Early_Burst_Up", "Late_Loss_Down"))
})

test_that("invalid template name errors clearly", {
  expect_error(plot_template_atlas("rna",
                                     templates = "Not_A_Pattern"),
               regexp = "Unknown template")
})

test_that("output_dir writes one PDF per template + reports paths", {
  tmp <- tempfile("atlas-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  out <- plot_template_atlas("rna",
                              templates = c("Early_Burst_Up",
                                              "Mid_Trough"),
                              output_dir = tmp)
  files <- attr(out, "files")
  expect_length(files, 2L)
  expect_true(all(file.exists(files)))
  expect_true(all(grepl("\\.pdf$", files)))
})
