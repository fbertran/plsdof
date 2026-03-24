testthat::test_that("pls.cv returns coherent plsdof object (coef/vcov)", {
  set.seed(42)
  n <- 40; p <- 6
  X <- matrix(rnorm(n*p), n, p)
  beta <- rnorm(p)
  y <- as.vector(scale(X, TRUE, TRUE) %*% beta + rnorm(n, sd = 0.2))

  m <- 4
  # compute.covariance only for use.kernel = FALSE
  res <- pls.cv(X, y, k = 5, m = m, compute.covariance = TRUE)
  testthat::expect_true(inherits(res, "plsdof"))
  testthat::expect_type(res$m.opt, "double")
  testthat::expect_true(res$m.opt >= 0 && res$m.opt <= m)
  # S3 methods should work
  cf <- coef(res)
  testthat::expect_type(cf, "double")
  testthat::expect_length(cf, p)
  V <- vcov(res)
  testthat::expect_true(is.matrix(V))
  testthat::expect_equal(dim(V), c(p, p))
  # covariance should be symmetric (within numerical tolerance)
  testthat::expect_lt(max(abs(V - t(V))), 1e-7)
})

testthat::test_that("pls.ic selects an m within bounds and returns coefficients", {
  set.seed(123)
  n <- 40; p <- 6
  X <- matrix(rnorm(n*p), n, p)
  beta <- rnorm(p)
  y <- as.vector(scale(X, TRUE, TRUE) %*% beta + rnorm(n, sd = 0.25))

  m <- 5
  res <- pls.ic(X, y, m = m, criterion = "bic", compute.jacobian = TRUE, verbose = FALSE)
  testthat::expect_true(inherits(res, "plsdof"))
  testthat::expect_true(res$m.opt >= 0 && res$m.opt <= m)
  testthat::expect_type(res$intercept, "double")
  testthat::expect_type(res$coefficients, "double")
  testthat::expect_length(res$coefficients, p)
  # If we asked for jacobian (linear case), covariance for the chosen m should be a p x p matrix
  V <- res$covariance
  if (!is.null(V)) {
    testthat::expect_equal(dim(V), c(p, p))
  }
})

testthat::test_that("linear.pls.fit returns expected shapes and DoF increases", {
  set.seed(7)
  n <- 30; p <- 5
  X <- matrix(rnorm(n*p), n, p)
  y <- rnorm(n)
  m <- 4

  fit <- linear.pls.fit(X, y, m = m, compute.jacobian = TRUE)

  # Yhat: n x (m+1), coefficients: p x (m+1), intercept: length m+1
  testthat::expect_equal(dim(fit$Yhat), c(n, m + 1))
  testthat::expect_equal(dim(fit$coefficients), c(p, m + 1))
  testthat::expect_length(fit$intercept, m + 1)
  testthat::expect_length(fit$DoF, m + 1)

  # DoF should be nondecreasing and bounded by min(n-1, p+1)
  testthat::expect_true(all(diff(fit$DoF) >= -1e-8))
  testthat::expect_true(max(fit$DoF) <= min(n - 1, p + 1) + 1e-8)
})

testthat::test_that("PCR basics: shapes and jackknife fields from pcr / pcr.cv", {
  set.seed(99)
  n <- 35; p <- 4
  X <- matrix(rnorm(n*p), n, p)
  y <- rnorm(n)
  m <- 3

  pr <- pcr(X, y, m = m)
  testthat::expect_equal(dim(pr$coefficients), c(p, m + 1))
  testthat::expect_length(pr$intercept, m + 1)

  pcv <- pcr.cv(X, y, k = 5, m = m)
  testthat::expect_true(pcv$m.opt >= 0 && pcv$m.opt <= m)
  testthat::expect_length(pcv$coefficients, p)
  testthat::expect_true(is.array(pcv$coefficients.jackknife) || is.list(pcv$coefficients.jackknife))
})

testthat::test_that("first.local.minimum returns an index within vector bounds", {
  idx <- first.local.minimum(c(3, 2, 1, 2, 5, 6))
  testthat::expect_true(is.numeric(idx) || is.integer(idx))
  testthat::expect_true(idx >= 1 && idx <= 6)
})
