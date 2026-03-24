# Lower bound for the Degrees of Freedom

This function computes the lower bound for the the Degrees of Freedom of
PLS with 1 component.

## Usage

``` r
compute.lower.bound(X)
```

## Arguments

- X:

  matrix of predictor observations.

## Value

- bound:

  logical. bound is `TRUE` if the decay of the eigenvalues is slow
  enough

- lower.bound:

  if bound is TRUE, this is the lower bound, otherwise, it is set to -1

## Details

If the decay of the eigenvalues of `cor(X)` is not too fast, we can
lower-bound the Degrees of Freedom of PLS with 1 component. Note that we
implicitly assume that we use scaled predictor variables to compute the
PLS solution.

## References

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial
Least Squares Regression". Journal of the American Statistical
Association 106 (494)
<https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107>

## See also

[`pls.model`](https://fbertran.github.io/plsdof/reference/pls.model.md)

## Author

Nicole Kraemer

## Examples

``` r
# Boston Housing data
library(MASS)
data(Boston)
X<-Boston[,-14]
my.lower<-compute.lower.bound(X)
```
