# Ridge Regression.

This function computes the optimal ridge regression model based on
cross-validation.

## Usage

``` r
ridge.cv(
  X,
  y,
  lambda = NULL,
  scale = TRUE,
  k = 10,
  plot.it = FALSE,
  groups = NULL,
  method.cor = "pearson",
  compute.jackknife = TRUE
)
```

## Arguments

- X:

  matrix of input observations. The rows of `X` contain the samples, the
  columns of `X` contain the observed variables

- y:

  vector of responses. The length of y must equal the number of rows of
  X

- lambda:

  Vector of penalty terms.

- scale:

  Scale the columns of X? Default is scale=TRUE.

- k:

  Number of splits in `k`-fold cross-validation. Default value is
  `k`=10.

- plot.it:

  Plot the cross-validation error as a function of `lambda`? Default is
  FALSE.

- groups:

  an optional vector with the same length as `y`. It encodes a
  partitioning of the data into distinct subgroups. If `groups` is
  provided, `k=10` is ignored and instead, cross-validation is performed
  based on the partioning. Default is `NULL`.

- method.cor:

  How should the correlation to the response be computed? Default is
  ”pearson”.

- compute.jackknife:

  Logical. If `TRUE`, the regression coefficients on each of the
  cross-validation splits is stored. Default is `TRUE`.

## Value

- cv.error.matrix:

  matrix of cross-validated errors based on mean squared error. A row
  corresponds to one cross-validation split.

- cv.error:

  vector of cross-validated errors based on mean squared error

- lambda.opt:

  optimal value of `lambda`, based on mean squared error

- intercept:

  intercept of the optimal model, based on mean squared error

- coefficients:

  vector of regression coefficients of the optimal model, based on mean
  squared error

- cor.error.matrix:

  matrix of cross-validated errors based on correlation. A row
  corresponds to one cross-validation split.

- cor.error:

  vector of cross-validated errors based on correlation

- lambda.opt.cor:

  optimal value of `lambda`, based on correlation

- intercept.cor:

  intercept of the optimal model, based on correlation

- coefficients.cor:

  vector of regression coefficients of the optimal model, based on mean
  squared error

- coefficients.jackknife:

  Array of the regression coefficients on each of the cross-validation
  splits. The dimension is `ncol(X) x length(lambda) x k`.

## Details

Based on the regression coefficients `coefficients.jackknife` computed
on the cross-validation splits, we can estimate their mean and their
variance using the jackknife. We remark that under a fixed design and
the assumption of normally distributed `y`-values, we can also derive
the true distribution of the regression coefficients.

## See also

[`pls.cv`](https://fbertran.github.io/plsdof/reference/pls.cv.md),
[`pcr.cv`](https://fbertran.github.io/plsdof/reference/pcr.cv.md),
[`benchmark.regression`](https://fbertran.github.io/plsdof/reference/benchmark.regression.md)

## Author

Nicole Kraemer

## Examples

``` r
n<-100 # number of observations
p<-60 # number of variables
X<-matrix(rnorm(n*p),ncol=p) 
y<-rnorm(n)
ridge.object<-ridge.cv(X,y)
```
