testthat::test_that("Datasets load", {
  testthat::skip_if_not_installed("plsdof")
  ds <- c()
  if (length(ds) == 0) testthat::skip("No datasets shipped")
  for (d in ds) {
    expect_silent(data(list=d, package="plsdof", envir = environment()))
    expect_true(exists(d, envir = environment()))
  }
})
