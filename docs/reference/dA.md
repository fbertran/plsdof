# Derivative of normalization function

This function computes the derivative of the function \$\$v\mapsto
\frac{w}{\\w\\\_A}\$\$ with respect to y.

## Usage

``` r
dA(w, A, dw)
```

## Arguments

- w:

  vector of length n.

- A:

  square matrix that defines the norm

- dw:

  derivative of w with respect to y. As y is a vector of length n, the
  derivative is a matrix of size nxn.

## Value

the Jacobian matrix of the normalization function. This is a matrix of
size nxn.

## Details

The first derivative of the normalization operator is
\$\$\frac{\partial}{\partial y}\left(w\mapsto
\frac{w}{\\w\\\_A}\right)=\frac{1}{\\w\\}\left(I_n - \frac{w w^ \top
A}{w^\top w}\right) \frac{\partial w}{\partial y}\$\$

## References

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial
Least Squares Regression". Journal of the American Statistical
Association 106 (494)
<https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107>

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom,
and Efficient Model Selection", Proceedings of the 24th International
Conference on Machine Learning, Omni Press, 441 - 448

## See also

[`normalize`](https://fbertran.github.io/plsdof/reference/normalize.md),
[`dnormalize`](https://fbertran.github.io/plsdof/reference/dnormalize.md)

## Author

Nicole Kraemer

## Examples

``` r
w<-rnorm(15)
dw<-diag(15)
A<-diag(1:15)
d.object<-dA(w,A,dw)
```
