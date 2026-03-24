# First derivative of the projection operator

This function computes the first derivative of the projection operator
\$\$P_V z= V V^\top z\$\$

## Usage

``` r
dvvtz(v, z, dv, dz)
```

## Arguments

- v:

  orthonormal basis of the space on which `z` is projected. `v` is
  either a matrix or a vector.

- z:

  vector that is projected onto the columns of `v`

- dv:

  first derivative of the the columns of `v` with respect to a vector y.
  If `v` is a matrix, `dv` is an array of dimension
  `ncol(v)`x`nrow(v)`x`length(y)`. If `v` is a vector, `dv` is a matrix
  of dimension `nrow(v)`x`length(y)`.

- dz:

  first derivative of `z` with respect to a vector y. This is a matrix
  of dimension `nrow(v)`x`length(y)`.

## Value

The first derivative of the projection operator with respect to y. This
is a matrix of dimension `nrow(v)`x`length(y)`.

## Details

For the computation of the first derivative, we assume that the columns
of `v` are normalized and mutually orthogonal. (Note that the function
will not return an error message if these assumptionsa are not
fulfilled. If we denote the columns of `v` by \\v_1,\ldots,v_l\\, the
first derivative of the projection operator is \$\$ \frac{\partial
P}{\partial y}=\sum\_{j=1} ^ l \left\[ \left(v_j z^ \top + v_j^ \top z
I_n \right)\frac{\partial v_j}{\partial y} + v_j v_j ^ \top
\frac{\partial z}{\partial y}\right\] \$\$ Here, n denotes the length of
the vectors \\v_j\\.

## Note

This is an internal function.

## References

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial
Least Squares Regression". Journal of the American Statistical
Association. 106 (494)
<https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107>

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom,
and Efficient Model Selection", Proceedings of the 24th International
Conference on Machine Learning, Omni Press, 441 - 448

## See also

[`vvtz`](https://fbertran.github.io/plsdof/reference/vvtz.md)

## Author

Nicole Kraemer, Mikio L. Braun
