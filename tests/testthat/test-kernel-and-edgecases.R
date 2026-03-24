testthat::test_that("kernel.pls.fit returns coherent shapes and DoF (compute.jacobian=TRUE)", {
  set.seed(123)
  n <- 30; p <- 8; m <- 5
  X <- matrix(rnorm(n*p), n, p)
  y <- rnorm(n)

  fit <- kernel.pls.fit(X, y, m = m, compute.jacobian = TRUE)

  testthat::expect_equal(dim(fit$Yhat), c(n, m + 1))
  testthat::expect_equal(dim(fit$coefficients), c(p, m + 1))
  testthat::expect_length(fit$intercept, m + 1)
  testthat::expect_length(fit$DoF, m + 1)
  # DoF monotone nondecreasing and bounded
  testthat::expect_true(all(diff(fit$DoF[1:m]) >= -1e-8))
  testthat::expect_true(max(fit$DoF) <= min(n - 1, p + 1) + 1e-8)
})

testthat::test_that("pls.model with use.kernel=TRUE returns predictions and shapes", {
  set.seed(456)
  n <- 40; p <- 10; m <- 6
  X <- matrix(rnorm(n*p), n, p)
  y <- rnorm(n)
  ntest <- 15
  Xtest <- matrix(rnorm(ntest*p), ntest, p)

  res <- pls.model(X, y, m = m, Xtest = Xtest, compute.DoF = TRUE, use.kernel = TRUE)

  testthat::expect_true(is.matrix(res$coefficients))
  testthat::expect_equal(dim(res$coefficients), c(p, m + 1))
  testthat::expect_length(res$intercept, m + 1)
  testthat::expect_length(res$DoF, m + 1)
  testthat::expect_equal(dim(res$Yhat), c(n, m + 1))
  # predictions on Xtest if given
  testthat::expect_true(is.matrix(res$prediction))
  testthat::expect_equal(dim(res$prediction), c(ntest, m + 1))
  # DoF bounded
  testthat::expect_true(max(res$DoF) <= min(n - 1, p + 1) + 1e-8)
})

testthat::test_that("pls.model keeps covariance when jacobian is requested without DoF", {
  set.seed(654)
  n <- 40; p <- 6; m <- 4
  X <- matrix(rnorm(n*p), n, p)
  y <- rnorm(n)

  res <- pls.model(X, y, m = m, compute.DoF = FALSE, compute.jacobian = TRUE)

  testthat::expect_equal(res$DoF, 1:(m + 1))
  testthat::expect_false(is.null(res$covariance))
  testthat::expect_equal(dim(res$covariance), c(m + 1, p, p))
})

testthat::test_that("High-dimensional case p >> n still returns bounded DoF", {
  set.seed(789)
  n <- 20; p <- 60; m <- 10
  X <- matrix(rnorm(n*p), n, p)
  y <- rnorm(n)

  fit <- linear.pls.fit(X, y, m = m, compute.jacobian = TRUE)

  testthat::expect_equal(dim(fit$coefficients), c(p, m + 1))
  testthat::expect_length(fit$DoF, m + 1)
  testthat::expect_true(max(fit$DoF) <= min(n - 1, p + 1) + 1e-8)
})

testthat::test_that("pls.cv respects 'groups' partitioning (k ignored)", {
  set.seed(2468)
  n <- 50; p <- 6; m <- 4
  X <- matrix(rnorm(n*p), n, p)
  y <- rnorm(n)
  # partition into 5 custom groups of unequal sizes
  groups <- rep(1:5, length.out = n); groups <- sample(groups)

  res <- pls.cv(X, y, k = 10, groups = groups, m = m, compute.covariance = FALSE)

  # When groups are supplied, rows of cv.error.matrix == #unique(groups)
  testthat::expect_true(is.matrix(res$cv.error.matrix))
  testthat::expect_equal(nrow(res$cv.error.matrix), length(unique(groups)))
  testthat::expect_true(res$m.opt >= 0 && res$m.opt <= m)
})

testthat::test_that("pls.cv(use.kernel=TRUE) works without covariance", {
  set.seed(1357)
  n <- 45; p <- 7; m <- 5
  X <- matrix(rnorm(n*p), n, p)
  y <- rnorm(n)

  res <- pls.cv(X, y, m = m, use.kernel = TRUE, compute.covariance = FALSE)

  testthat::expect_true(inherits(res, "plsdof"))
  testthat::expect_type(res$coefficients, "double")
  testthat::expect_length(res$coefficients, p)
  # covariance is only implemented for use.kernel=FALSE, so either NULL or absent
  if (!is.null(res$covariance)) {
    testthat::expect_true(all(dim(res$covariance) != 0))  # not asserting shape here
  }
})
