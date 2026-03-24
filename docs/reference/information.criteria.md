# Information criteria

This function computes the optimal model parameters using three
different model selection criteria (aic, bic, gmdl).

## Usage

``` r
information.criteria(RSS, DoF, yhat = NULL, sigmahat, n, criterion = "bic")
```

## Arguments

- RSS:

  vector of residual sum of squares.

- DoF:

  vector of Degrees of Freedom. The length of `DoF` is the same as the
  length of `RSS`.

- yhat:

  vector of squared norm of yhat. The length of `yhat` is the same as
  the length of `RSS`. It is only needed for gmdl. Default value is
  `NULL`.

- sigmahat:

  Estimated model error. The length of `sigmahat` is the same as the
  length of `RSS`.

- n:

  number of observations.

- criterion:

  one of the options "aic", "bic" and "gmdl".

## Value

- DoF:

  degrees of freedom

- score:

  vector of the model selection criterion

- par:

  index of the first local minimum of `score`

## Details

The Akaike information criterion (aic) is defined as \$\${aic}=
\frac{{RSS}}{n} + 2\frac{{DoF}}{n} \sigma^ 2\\.\$\$ The Bayesian
information criterion (bic) is defined as \$\${bic}= \frac{{RSS}}{n} +
log(n)\frac{{DoF}}{n} \sigma^ 2\\.\$\$ The generalized minimum
description length (gmdl) is defined as
\$\$gmdl=\frac{n}{2}log(S)+\frac{DoF}{2}log(F)+\frac{1}{2}log(n)\$\$
with \$\$S=\hat \sigma ^2\$\$ Note that it is also possible to use the
function `information.criteria` for other regression methods than
Partial Least Squares.

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

[`pls.ic`](https://fbertran.github.io/plsdof/reference/pls.ic.md)

## Author

Nicole Kraemer, Mikio Braun

## Examples

``` r
## This is an internal function called by pls.ic

```
