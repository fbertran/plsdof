if (requireNamespace("testthat", quietly = TRUE)) {
  testthat::test_check("plsdof")
} else {
  message("Skipping tests: 'testthat' not installed.")
}
