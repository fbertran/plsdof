# Index of the first local minimum.

This function computes the index of the first local minimum.

## Usage

``` r
first.local.minimum(x)
```

## Arguments

- x:

  vector.

## Value

the index of the first local minimum of `x`.

## References

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial
Least Squares Regression". Journal of the American Statistical
Association. ahead of print 106 (494)
<https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107>

## Author

Nicole Kraemer

## Examples

``` r
v<-rnorm(30)
out<-first.local.minimum(v)
```
