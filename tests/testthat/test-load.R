testthat::test_that("Package loads quietly", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("plsdof")
  expect_silent(suppressPackageStartupMessages(library(plsdof)))
  expect_true("plsdof" %in% loadedNamespaces())
  detach("package:plsdof", character.only = TRUE)
})
