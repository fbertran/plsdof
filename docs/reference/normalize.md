# Normalization of vectors

Normalization of vectors.

## Usage

``` r
normalize(v, w = NULL)
```

## Arguments

- v:

  vector

- w:

  optional vector

## Value

- v:

  normalized `v`

- w:

  normalized `w`

## Details

The vector `v` is normalized to length 1. If `w` is given, it is
normalized by the length of `v`.

## Author

Nicole Kraemer, Mikio L. Braun

## Examples

``` r
v<-rnorm(5)
w<-rnorm(10)
dummy<-normalize(v,w)
```
