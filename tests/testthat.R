if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  library(pdactrace)
  test_check("pdactrace")
} else {
  message("Package 'testthat' is not available; skipping tests.")
}
