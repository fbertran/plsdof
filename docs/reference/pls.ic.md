# Model selection for Partial Least Squares based on information criteria

This function computes the optimal model parameters using one of three
different model selection criteria (aic, bic, gmdl) and based on two
different Degrees of Freedom estimates for PLS.

## Usage

``` r
pls.ic(
  X,
  y,
  m = min(ncol(X), nrow(X) - 1),
  criterion = "bic",
  naive = FALSE,
  use.kernel = FALSE,
  compute.jacobian = FALSE,
  verbose = TRUE
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
  `m=ncol(X)`.

- criterion:

  Choice of the model selection criterion. One of the three options aic,
  bic, gmdl.

- naive:

  Use the naive estimate for the Degrees of Freedom? Default is `FALSE`.

- use.kernel:

  Use kernel representation? Default is `use.kernel=FALSE`.

- compute.jacobian:

  Should the first derivative of the regression coefficients be computed
  as well? Default is `FALSE`

- verbose:

  If `TRUE`, the function prints a warning if the algorithms produce
  negative Degrees of Freedom. Default is `TRUE`.

## Value

The function returns an object of class "plsdof".

- DoF:

  Degrees of Freedom

- m.opt:

  optimal number of components

- sigmahat:

  vector of estimated model errors

- intercept:

  intercept

- coefficients:

  vector of regression coefficients

- covariance:

  if `compute.jacobian=TRUE` and `use.kernel=FALSE`, the function
  returns the covariance matrix of the optimal regression coefficients.

- m.crash:

  the number of components for which the algorithm returns negative
  Degrees of Freedom

## Details

There are two options to estimate the Degrees of Freedom of PLS:
`naive=TRUE` defines the Degrees of Freedom as the number of components
+1, and `naive=FALSE` uses the generalized notion of Degrees of Freedom.
If `compute.jacobian=TRUE`, the function uses the Lanczos decomposition
to derive the Degrees of Freedom, otherwise, it uses the Krylov
representation. (See Kraemer and Sugiyama (2011) for details.) The
latter two methods only differ with respect to the estimation of the
noise level.

## References

Akaikie, H. (1973) "Information Theory and an Extension of the Maximum
Likelihood Principle". Second International Symposium on Information
Theory, 267 - 281.

Hansen, M., Yu, B. (2001). "Model Selection and Minimum Descripion
Length Principle". Journal of the American Statistical Association, 96,
746 - 774

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial
Least Squares Regression". Journal of the American Statistical
Association 106 (494)
<https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107>

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom,
and Efficient Model Selection", Proceedings of the 24th International
Conference on Machine Learning, Omni Press, 441 - 448

Schwartz, G. (1979) "Estimating the Dimension of a Model" Annals of
Statistics 26(5), 1651 - 1686.

## See also

[`pls.model`](https://fbertran.github.io/plsdof/reference/pls.model.md),
[`pls.cv`](https://fbertran.github.io/plsdof/reference/pls.cv.md)

## Author

Nicole Kraemer, Mikio L. Braun

## Examples

``` r
n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)

# compute linear PLS
pls.object<-pls.ic(X,y,m=ncol(X))
```
