# Projectin operator

This function computes the projection operator \$\$P_V z= V V^\top z\$\$

## Usage

``` r
vvtz(v, z)
```

## Arguments

- v:

  orthonormal basis of the space on which `z` is projected. `v` is
  either a matrix or a vector.

- z:

  vector that is projected onto the columns of `v`

## Value

value of the projection operator

## Details

The above formula is only valid if the columns of `v` are normalized and
mutually orthogonal.

## See also

[`dvvtz`](https://fbertran.github.io/plsdof/reference/dvvtz.md)

## Author

Nicole Kraemer

## Examples

``` r
# generate random orthogonal vectors
X<-matrix(rnorm(10*100),ncol=10)   # random data
S<-cor(X)         # correlation matrix of data
v<-eigen(S)$vectors[,1:3]    # first three eigenvectors of correlation matrix
z<-rnorm(10)        # random vector z
projection.z<-vvtz(v,z)

```
