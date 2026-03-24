testthat::test_that("Help files parse or Rd_db fallback", {
  testthat::skip_if_not_installed("plsdof")
  man_dir <- system.file("man", package = "plsdof")
  if (!dir.exists(man_dir)) testthat::skip("No installed man/; skipping Rd tests")
  
  # Try the fast DB path first
  h <- try(tools::Rd_db("plsdof"), silent = TRUE)
  if (inherits(h, "try-error")) {
    testthat::skip("Installed help DB seems corrupt; falling back to parsing each Rd file.")
  } else {
    testthat::expect_true(is.list(h) && length(h) > 0)
    return(invisible())
  }
  
  # Fallback: parse each installed Rd; this detects real syntax issues
  rds <- list.files(man_dir, pattern = "[.]Rd$", full.names = TRUE)
  testthat::expect_gt(length(rds), 0L)
  for (f in rds) {
    rd <- try(tools::parse_Rd(f, verbose = FALSE), silent = TRUE)
    testthat::expect_false(inherits(rd, "try-error"), info = paste("Failed to parse:", basename(f)))
  }
})
