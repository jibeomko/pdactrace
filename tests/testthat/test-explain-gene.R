test_that("explain_gene view='evidence' prints provenance, no math sections", {
  out <- capture.output(
    res <- explain_gene("LGALS3BP", view = "evidence", verbose = TRUE),
    type = "output")
  joined <- paste(out, collapse = "\n")
  expect_match(joined, "Evidence \\(plain-English provenance\\)")
  expect_false(grepl("Trajectory fit", joined))
  expect_false(grepl("Effect magnitude", joined))
  expect_type(res, "list")
  expect_equal(res$view, "evidence")
  expect_null(res$math)
})

test_that("explain_gene view='math' prints all 8 math sections", {
  out <- capture.output(
    res <- explain_gene("LGALS3BP", view = "math", verbose = TRUE),
    type = "output")
  joined <- paste(out, collapse = "\n")
  for (section in c("Trajectory fit", "Effect magnitude",
                    "Cohort consistency", "RNA-protein coupling",
                    "Serum bridge", "Cell specificity",
                    "Filter survival", "Clinical role")) {
    expect_match(joined, section, fixed = TRUE)
  }
  expect_false(grepl("Evidence \\(plain-English provenance\\)", joined))
  expect_equal(res$view, "math")
  expect_type(res$math, "list")
})

test_that("explain_gene view='both' prints both sections separated by a rule", {
  out <- capture.output(
    res <- explain_gene("LGALS3BP", view = "both", verbose = TRUE),
    type = "output")
  joined <- paste(out, collapse = "\n")
  expect_match(joined, "Evidence \\(plain-English provenance\\)")
  expect_match(joined, "Trajectory fit", fixed = TRUE)
  expect_match(joined, "----")  # horizontal rule between sections
  expect_equal(res$view, "both")
})

test_that("explain_gene verbose=FALSE returns the list silently", {
  out <- capture.output(
    res <- explain_gene("LGALS3BP", view = "math", verbose = FALSE),
    type = "output")
  expect_equal(length(out), 0L)
  expect_type(res, "list")
  expect_equal(res$gene, "LGALS3BP")
  expect_type(res$math, "list")
})
