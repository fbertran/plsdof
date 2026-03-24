# Krylov sequence

This function computes the Krylov sequence of a matrix and a vector.

## Usage

``` r
krylov(A, b, m)
```

## Arguments

- A:

  square matrix of dimension p x p.

- b:

  vector of length p

- m:

  length of the Krylov sequence

## Value

A matrix of size p x m containing the sequence b,Ab,..., A^(m-1)b.

## Author

Nicole Kraemer

## Examples

``` r
A<-matrix(rnorm(8*8),ncol=8)
b<-rnorm(8)
K<-krylov(A,b,4)
```
