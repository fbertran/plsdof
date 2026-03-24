# Comparison of model selection criteria for Partial Least Squares Regression.

This function computes the test error over several runs for different
model selection strategies.

## Usage

``` r
benchmark.pls(
  X,
  y,
  m = ncol(X),
  R = 20,
  ratio = 0.8,
  verbose = TRUE,
  k = 10,
  ratio.samples = 1,
  use.kernel = FALSE,
  criterion = "bic",
  true.coefficients = NULL
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

- R:

  number of runs. Default is 20.

- ratio:

  ratio no of training examples/(no of training examples + no of test
  examples). Default is 0.8

- verbose:

  If `TRUE`, the functions plots the progress of the function. Default
  is `TRUE`.

- k:

  number of cross-validation splits. Default is 10.

- ratio.samples:

  Ratio of (no of training examples + no of test examples)/`nrow(X)`.
  Default is 1.

- use.kernel:

  Use kernel representation? Default is `use.kernel=FALSE`.

- criterion:

  Choice of the model selection criterion. One of the three options aic,
  bic, gmdl. Default is "bic".

- true.coefficients:

  The vector of true regression coefficients (without intercept), if
  available. Default is `NULL`.

## Value

- MSE:

  data frame of size R x 5. It contains the test error for the five
  different methods for each of the R runs.

- M:

  data frame of size R x 5. It contains the optimal number of components
  for the five different methods for each of the R runs.

- DoF:

  data frame of size R x 5. It contains the Degrees of Freedom
  (corresponding to `M`) for the five different methods for each of the
  R runs.

- TIME:

  data frame of size R x 4. It contains the runtime for all methods
  (apart from the zero model) for each of the R runs.

- M.CRASH:

  data frame of size R x 2. It contains the number of components for
  which the Krylov representation and the Lanczos representation return
  negative Degrees of Freedom, hereby indicating numerical problems.

- ME:

  if `true.coefficients` are available, this is a data frame of size R
  x 5. It contains the model error for the five different methods for
  each of the R runs.

- SIGMAHAT:

  data frame of size R x 5. It contains the estimation of the noise
  level provided by the five different methods for each of the R runs.

## Details

The function estimates the optimal number of PLS components based on
four different criteria: (1) cross-validation, (2) information criteria
with the naive Degrees of Freedom DoF(m)=m+1, (3) information criteri
with the Degrees of Freedom computed via a Lanczos represenation of PLS
and (4) information criteri with the Degrees of Freedom computed via a
Krylov represenation of PLS. Note that the latter two options only
differ with respect to the estimation of the model error.

In addition, the function computes the test error of the "zero model",
i.e. `mean(y)` on the training data is used for prediction.

If `true.coefficients` are available, the function also computes the
model error for the different methods, i.e. the sum of squared
differences between the true and the estimated regression coefficients.

## References

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial
Least Squares Regression". Journal of the American Statistical
Association 106 (494)
<https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107>

## See also

[`pls.ic`](https://fbertran.github.io/plsdof/reference/pls.ic.md),
[`pls.cv`](https://fbertran.github.io/plsdof/reference/pls.cv.md)

## Author

Nicole Kraemer

## Examples

``` r
# generate artificial data
n<-50 # number of examples
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
true.coefficients<-runif(p,1,3)
y<-X%*%true.coefficients + rnorm(n,0,5)
my.benchmark<-benchmark.pls(X,y,R=10,true.coefficients=true.coefficients)
#> iteration no  1  
#> iteration no  2  
#> iteration no  3  
#> iteration no  4  
#> iteration no  5  
#> iteration no  6  
#> iteration no  7  
#> iteration no  8  
#> iteration no  9  
#> iteration no  10  
```
