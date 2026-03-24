# Computation of the Degrees of Freedom

This function computes the Degrees of Freedom using the Krylov
representation of PLS.

## Usage

``` r
pls.dof(pls.object, n, y, K, m, DoF.max)
```

## Arguments

- pls.object:

  object returned by `linear.pls.fit` or by `kernel.pls.fit`

- n:

  number of observations

- y:

  vector of response observations.

- K:

  kernel matrix X X^t.

- m:

  number of components

- DoF.max:

  upper bound on the Degrees of Freedom.

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

- TT:

  matrix of normalized PLS components

## Details

This computation of the Degrees of Freedom is based on the equivalence
of PLS regression and the projection of the response vector `y` onto the
Krylov space spanned by \$\$Ky,K^2 y,...,K^m y.\$\$ Details can be found
in Kraemer and Sugiyama (2011).

## References

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial
Least Squares Regression". Journal of the American Statistical
Association 106 (494)
<https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107>

Kraemer, N., Sugiyama M., Braun, M.L. (2009) "Lanczos Approximations for
the Speedup of Kernel Partial Least Squares Regression." Proceedings of
the Twelfth International Conference on Artificial Intelligence and
Statistics (AISTATS), p. 272-279

## See also

[`pls.model`](https://fbertran.github.io/plsdof/reference/pls.model.md),
[`pls.ic`](https://fbertran.github.io/plsdof/reference/pls.ic.md)

## Author

Nicole Kraemer, Mikio L. Braun

## Examples

``` r
# this is an internal function
```
