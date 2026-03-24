# Principal Components Regression

This function computes the Principal Components Regression (PCR) fit.

## Usage

``` r
pcr(
  X,
  y,
  scale = TRUE,
  m = min(ncol(X), nrow(X) - 1),
  eps = 1e-06,
  supervised = FALSE
)
```

## Arguments

- X:

  matrix of predictor observations.

- y:

  vector of response observations. The length of `y` is the same as the
  number of rows of `X`.

- scale:

  Should the predictor variables be scaled to unit variance? Default is
  `TRUE`.

- m:

  maximal number of principal components. Default is
  `m=min(ncol(X),nrow(X)-1)`.

- eps:

  precision. Eigenvalues of the correlation matrix of `X` that are
  smaller than `eps` are set to 0. The default value is `eps=10^{-6}.`

- supervised:

  Should the principal components be sorted by decreasing squared
  correlation to the response? Default is FALSE.

## Value

- coefficients:

  matrix of regression coefficients, including the coefficients of the
  null model, i.e. the constant model `mean(y)`.

- intercept:

  vector of intercepts, including the intercept of the null model, i.e.
  the constant model `mean(y)`.

## Details

The function first scales all predictor variables to unit variance, and
then computes the PCR fit for all components. Is `supervised=TRUE`, we
sort the principal correlation according to the squared correlation to
the response.

## See also

[`pcr.cv`](https://fbertran.github.io/plsdof/reference/pcr.cv.md),
[`pls.cv`](https://fbertran.github.io/plsdof/reference/pls.cv.md)

## Author

Nicole Kraemer

## Examples

``` r
n<-50 # number of observations
p<-15 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)

my.pcr<-pcr(X,y,m=10)

```
