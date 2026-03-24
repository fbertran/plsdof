# Derivative of normalization function

This function computes the derivative of the function \$\$v\mapsto
\frac{v}{\\v\\}\$\$ with respect to y.

## Usage

``` r
dnormalize(v, dv)
```

## Arguments

- v:

  vector of length n.

- dv:

  derivative of v with respect to y. As y is a vector of length n, the
  derivative is a matrix of size nxn.

## Value

the Jacobian matrix of the normalization function. This is a matrix of
size nxn.

## Details

The first derivative of the normalization operator is
\$\$\frac{\partial}{\partial y}\left(v\mapsto
\frac{v}{\\v\\}\right)=\frac{1}{\\v\\}\left(I_n - \frac{v v^
\top}{v^\top v}\right) \frac{\partial v}{\partial y}\$\$

## References

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial
Least Squares Regression". Journal of the American Statistical
Association 106 (494)
<https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107>

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom,
and Efficient Model Selection", Proceedings of the 24th International
Conference on Machine Learning, Omni Press, 441 - 448

## See also

[`normalize`](https://fbertran.github.io/plsdof/reference/normalize.md)

## Author

Nicole Kraemer, Mikio L. Braun

## Examples

``` r
v<-rnorm(15)
dv<-diag(15)
d.object<-dnormalize(v,dv)
```
