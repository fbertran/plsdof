# Kernel Partial Least Squares Fit

This function computes the Partial Least Squares fit. This algorithm
scales mainly in the number of observations.

## Usage

``` r
kernel.pls.fit(
  X,
  y,
  m = ncol(X),
  compute.jacobian = FALSE,
  DoF.max = min(ncol(X) + 1, nrow(X) - 1)
)
```

## Arguments

- X:

  matrix of predictor observations.

- y:

  vector of response observations. The length of `y` is the same as the
  number of rows of `X`.

- m:

  maximal number of Partial Least Squares components. Default is
  `m`=ncol(X).

- compute.jacobian:

  Should the first derivative of the regression coefficients be computed
  as well? Default is `FALSE`

- DoF.max:

  upper bound on the Degrees of Freedom. Default is
  `min(ncol(X)+1,nrow(X)-1)`.

## Value

- coefficients:

  matrix of regression coefficients

- intercept:

  vector of regression intercepts

- DoF:

  Degrees of Freedom

- sigmahat:

  vector of estimated model error

- Yhat:

  matrix of fitted values

- yhat:

  vector of squared length of fitted values

- RSS:

  vector of residual sum of error

- covariance:

  `NULL` object.

- TT:

  matrix of normalized PLS components

## Details

We first standardize `X` to zero mean and unit variance.

## References

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial
Least Squares Regression". Journal of the American Statistical
Association 106 (494)
<https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107>

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom,
and Efficient Model Selection", Proceedings of the 24th International
Conference on Machine Learning, Omni Press, 441 - 448

## See also

[`linear.pls.fit`](https://fbertran.github.io/plsdof/reference/linear.pls.fit.md),
[`pls.cv`](https://fbertran.github.io/plsdof/reference/pls.cv.md),[`pls.model`](https://fbertran.github.io/plsdof/reference/pls.model.md),
[`pls.ic`](https://fbertran.github.io/plsdof/reference/pls.ic.md)

## Author

Nicole Kraemer, Mikio L. Braun

## Examples

``` r
n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)


pls.object<-kernel.pls.fit(X,y,m=5,compute.jacobian=TRUE)


```
