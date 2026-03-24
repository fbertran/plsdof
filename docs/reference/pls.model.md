# Partial Least Squares

This function computes the Partial Least Squares fit.

## Usage

``` r
pls.model(
  X,
  y,
  m = ncol(X),
  Xtest = NULL,
  ytest = NULL,
  compute.DoF = FALSE,
  compute.jacobian = FALSE,
  use.kernel = FALSE,
  method.cor = "pearson"
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
  `m=min(ncol(X),nrow(X)-1)`.

- Xtest:

  optional matrix of test observations. Default is `Xtest=NULL`.

- ytest:

  optional vector of test observations. Default is `ytest=NULL`.

- compute.DoF:

  Logical variable. If `compute.DoF=TRUE`, the Degrees of Freedom of
  Partial Least Squares are computed. Default is `compute.DoF=FALSE`.

- compute.jacobian:

  Should the first derivative of the regression coefficients be computed
  as well? This remains independent from `compute.DoF`; callers can
  request covariance estimates while still returning the naive
  `1:(m + 1)` Degrees of Freedom sequence. Default is `FALSE`

- use.kernel:

  Should the kernel representation be used to compute the solution.
  Default is `FALSE`.

- method.cor:

  How should the correlation to the response be computed? Default is
  ”pearson”.

## Value

- coefficients:

  matrix of regression coefficients

- intercept:

  vector of intercepts

- DoF:

  vector of Degrees of Freedom

- RSS:

  vector of residual sum of error

- sigmahat:

  vector of estimated model error

- Yhat:

  matrix of fitted values

- yhat:

  vector of squared length of fitted values

- covariance:

  if `compute.jacobian` is `TRUE`, the function returns the array of
  covariance matrices for the PLS regression coefficients.

`prediction`if `Xtest` is provided, the predicted y-values for `Xtest`.
`mse`if `Xtest` and `ytest` are provided, the mean squared error on the
test data. `cor`if `Xtest` and `ytest` are provided, the correlation to
the response on the test data.

## Details

This function computes the Partial Least Squares fit and its Degrees of
Freedom. Further, it returns the regression coefficients and various
quantities that are needed for model selection in combination with
`information.criteria`.

## References

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial
Least Squares Regression". Journal of the American Statistical
Association 106 (494)
<https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107>

Kraemer, N., Sugiyama, M., Braun, M.L. (2009) "Lanczos Approximations
for the Speedup of Partial Least Squares Regression", Proceedings of the
12th International Conference on Artificial Intelligence and
Stastistics, 272 - 279

## See also

[`pls.ic`](https://fbertran.github.io/plsdof/reference/pls.ic.md),
[`pls.cv`](https://fbertran.github.io/plsdof/reference/pls.cv.md)

## Author

Nicole Kraemer, Mikio L. Braun

## Examples

``` r
n<-50 # number of observations
p<-15 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)

ntest<-200 #
Xtest<-matrix(rnorm(ntest*p),ncol=p) # test data
ytest<-rnorm(ntest) # test data

# compute PLS + degrees of freedom + prediction on Xtest
first.object<-pls.model(X,y,compute.DoF=TRUE,Xtest=Xtest,ytest=NULL)

# compute PLS + test error
second.object=pls.model(X,y,m=10,Xtest=Xtest,ytest=ytest)
```
