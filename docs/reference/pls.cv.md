# Model selection for Partial Least Squares based on cross-validation

This function computes the optimal model parameter using
cross-validation.

## Usage

``` r
pls.cv(
  X,
  y,
  k = 10,
  groups = NULL,
  m = ncol(X),
  use.kernel = FALSE,
  compute.covariance = FALSE,
  method.cor = "pearson"
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

- groups:

  an optional vector with the same length as `y`. It encodes a
  partitioning of the data into distinct subgroups. If `groups` is
  provided, `k=10` is ignored and instead, cross-validation is performed
  based on the partioning. Default is `NULL`.

- m:

  maximal number of Partial Least Squares components. Default is
  `m=ncol(X)`.

- use.kernel:

  Use kernel representation? Default is `use.kernel=FALSE`.

- compute.covariance:

  If `TRUE`, the function computes the covariance for the cv-optimal
  regression coefficients.

- method.cor:

  How should the correlation to the response be computed? Default is
  ”pearson”.

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

  vector of regression coefficients of the optimal model, based on mean
  squared error

- covariance:

  If `TRUE` and `use.kernel=FALSE`, the covariance of the cv-optimal
  regression coefficients (based on mean squared error) is returned.

## Details

The data are centered and scaled to unit variance prior to the PLS
algorithm. It is possible to estimate the covariance matrix of the
cv-optimal regression coefficients (`compute.covariance=TRUE`).
Currently, this is only implemented if `use.kernel=FALSE`.

## References

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial
Least Squares Regression". Journal of the American Statistical
Association 106 (494)
<https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107>

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom,
and Efficient Model Selection", Proceedings of the 24th International
Conference on Machine Learning, Omni Press, 441 - 448

## See also

[`pls.model`](https://fbertran.github.io/plsdof/reference/pls.model.md),
[`pls.ic`](https://fbertran.github.io/plsdof/reference/pls.ic.md)

## Author

Nicole Kraemer, Mikio L. Braun

## Examples

``` r
n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)

# compute linear PLS
pls.object<-pls.cv(X,y,m=ncol(X))

# define random partioning
groups<-sample(c("a","b","c"),n,replace=TRUE)
pls.object1<-pls.cv(X,y,groups=groups)
```
