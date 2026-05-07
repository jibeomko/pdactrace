# Marker-class regression tests --------------------------------------
# Pin the audit_class assignment for the four canonical case-study
# genes that the BiB manuscript and the audit-framework vignette use
# as illustrative examples. If any of these flips because the audit
# rule, the bundled atlas, or the gate logic changes, this test fires
# loudly — making the v0.99.x audit_class output a stable contract.

test_that("LGALS3BP stays high_confidence (v0.99.x atlas)", {
  ex <- explain_score("LGALS3BP", verbose = FALSE)
  expect_equal(ex$audit_class, "high_confidence")
  expect_gt(ex$audit_score, 0.85)
})

test_that("LTBP1 stays supported_uncertain (Class B exemplar)", {
  ex <- explain_score("LTBP1", verbose = FALSE)
  expect_equal(ex$audit_class, "supported_uncertain")
  expect_gt(ex$audit_score, 0.30)
  expect_lt(ex$audit_score, 0.85)
})

test_that("ALB stays penalized (plasma high-abundance gate)", {
  ex <- explain_score("ALB", verbose = FALSE)
  expect_equal(ex$audit_class, "penalized")
  expect_lt(ex$gates$value[ex$gates$gate == "leakage_gate"], 1)
})

test_that("GAPDH stays excluded (housekeeping gate)", {
  ex <- explain_score("GAPDH", verbose = FALSE)
  expect_equal(ex$audit_class, "excluded")
  expect_equal(ex$audit_score, 0)
  expect_equal(ex$gates$value[ex$gates$gate == "leakage_gate"], 0)
})

test_that("compare_candidates orders the four case studies as expected", {
  cmp <- compare_candidates(c("LGALS3BP", "LTBP1", "ALB", "GAPDH"))
  expect_equal(cmp$gene_symbol[1L], "LGALS3BP")
  expect_equal(tail(cmp$gene_symbol, 1L), "GAPDH")
})

test_that("toy data round-trips through fit_stage_de + classify_trajectory", {
  skip_if_not_installed("DESeq2")
  fit <- fit_stage_de(toy_counts,
                       stage  = toy_coldata$stage,
                       cohort = toy_coldata$cohort)
  expect_s3_class(fit, "data.table")
  pat <- classify_trajectory(fit)
  expect_s3_class(pat, "data.table")
})

test_that("list_data_sources returns the bundled accession table", {
  ds <- list_data_sources()
  expect_s3_class(ds, "data.table")
  expect_gt(nrow(ds), 20L)
  expect_setequal(
    intersect(c("RNA", "Protein", "scRNA", "Serum"), ds$layer),
    c("RNA", "Protein", "scRNA", "Serum"))
  rna_only <- list_data_sources(layer = "RNA")
  expect_true(all(rna_only$layer == "RNA"))
})

test_that("atlas_provenance carries both Zenodo DOIs", {
  ap <- atlas_provenance()
  expect_equal(ap$package_doi,    "10.5281/zenodo.20069896")
  expect_equal(ap$manuscript_doi, "10.5281/zenodo.20067849")
  expect_match(ap$package_repo, "github.com/jibeomko/pdactrace$")
})
