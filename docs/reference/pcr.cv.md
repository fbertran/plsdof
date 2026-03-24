# Model selection for Princinpal Components regression based on cross-validation

This function computes the optimal model parameter using
cross-validation. Mdel selection is based on mean squared error and
correlation to the response, respectively.

## Usage

``` r
pcr.cv(
  X,
  y,
  k = 10,
  m = min(ncol(X), nrow(X) - 1),
  groups = NULL,
  scale = TRUE,
  eps = 1e-06,
  plot.it = FALSE,
  compute.jackknife = TRUE,
  method.cor = "pearson",
  supervised = FALSE
)
```

## Arguments

- X:

  matrix of predictor observations.

- y:

  vector of response observations. The length of `y` is the same as the
  number of rows of `X`.

- k:

  number of cross-validation splits. Default is 10.

- m:

  maximal number of principal components. Default is
  `m=min(ncol(X),nrow(X)-1)`.

- groups:

  an optional vector with the same length as `y`. It encodes a
  partitioning of the data into distinct subgroups. If `groups` is
  provided, `k=10` is ignored and instead, cross-validation is performed
  based on the partioning. Default is `NULL`.

- scale:

  Should the predictor variables be scaled to unit variance? Default is
  `TRUE`.

- eps:

  precision. Eigenvalues of the correlation matrix of `X` that are
  smaller than `eps` are set to 0. The default value is `eps=10^{-6}.`

- plot.it:

  Logical. If `TRUE`, the function plots the cross-validation-error as a
  function of the number of components. Default is `FALSE`.

- compute.jackknife:

  Logical. If `TRUE`, the regression coefficients on each of the
  cross-validation splits is stored. Default is `TRUE`.

- method.cor:

  How should the correlation to the response be computed? Default is
  ”pearson”.

- supervised:

  Should the principal components be sorted by decreasing squared
  correlation to the response? Default is FALSE.

## Value

- cv.error.matrix:

  matrix of cross-validated errors based on mean squared error. A row
  corresponds to one cross-validation split.

- cv.error:

  vector of cross-validated errors based on mean squared error

- m.opt:

  optimal number of components based on mean squared error

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

- m.opt.cor:

  optimal number of components based on correlation

- intercept.cor:

  intercept of the optimal model, based on correlation

- coefficients.cor:

  vector of regression coefficients of the optimal model, based on
  correlation

- coefficients.jackknife:

  Array of the regression coefficients on each of the cross-validation
  splits, if `compute.jackknife=TRUE`. In this case, the dimension is
  `ncol(X) x (m+1) x k`.

## Details

The function computes the principal components on the scaled predictors.
Based on the regression coefficients `coefficients.jackknife` computed
on the cross-validation splits, we can estimate their mean and their
variance using the jackknife. We remark that under a fixed design and
the assumption of normally distributed `y`-values, we can also derive
the true distribution of the regression coefficients.

## See also

[`pls.model`](https://fbertran.github.io/plsdof/reference/pls.model.md),
[`pls.ic`](https://fbertran.github.io/plsdof/reference/pls.ic.md)

## Author

Nicole Kraemer, Mikio L. Braun

## Examples

``` r
n<-500 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)

# compute PCR 
pcr.object<-pcr.cv(X,y,scale=FALSE,m=3)
pcr.object1<-pcr.cv(X,y,groups=sample(c(1,2,3),n,replace=TRUE),m=3)
```
