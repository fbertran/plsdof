testthat::test_that("Exported symbols exist", {
  testthat::skip_if_not_installed("plsdof")
  suppressPackageStartupMessages(library(plsdof))
  exports <- c("benchmark.pls", "benchmark.regression", "compute.lower.bound", "dA", "dnormalize", "dvvtz", "first.local.minimum", "information.criteria", "kernel.pls.fit", "krylov", "linear.pls.fit", "normalize", "pcr", "pcr.cv", "pls.cv", "pls.dof", "pls.ic", "pls.model", "ridge.cv", "tr", "vvtz")
  for (nm in exports) {
    if (nzchar(nm)) {
      testthat::expect_true(exists(nm, where=asNamespace("plsdof")), info = paste("missing:", nm))
    }
  }
  detach("package:plsdof", character.only = TRUE)
})

testthat::test_that("S3 methods listed in NAMESPACE are findable", {
  testthat::skip_if_not_installed("plsdof")
  suppressPackageStartupMessages(library(plsdof))
  s3pairs <- list(c("coef","plsdof"), c("vcov","plsdof"))
  for (p in s3pairs) {
    if (length(p) == 2) {
      gen <- p[[1]]; cls <- p[[2]]
      meth <- utils::getS3method(gen, cls, optional = TRUE)
      testthat::expect_true(!is.null(meth), info = paste("S3 method not found:", gen, cls))
    }
  }
  detach("package:plsdof", character.only = TRUE)
})

testthat::test_that("exportClasses & exportMethods objects are present", {
  testthat::skip_if_not_installed("plsdof")
  suppressPackageStartupMessages(library(plsdof))
  classes <- c()
  for (cl in classes) {
    if (nzchar(cl)) testthat::expect_true(existsClass(cl), info = paste("Class missing:", cl))
  }
  methods <- c()
  # Can't easily assert methods without package-specific context; basic presence check:
  for (m in methods) {
    if (nzchar(m)) testthat::expect_true(any(getGenerics()@.Data == m), 
        info = paste("Generic not registered for method:", m))
  }
  detach("package:plsdof", character.only = TRUE)
})
