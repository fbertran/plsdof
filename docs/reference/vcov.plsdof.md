# Variance-covariance matrix

This function returns the variance-covariance matrix of a plsdof-object.

## Usage

``` r
# S3 method for class 'plsdof'
vcov(object, ...)
```

## Arguments

- object:

  an object of class "plsdof" that is returned by the function
  `linear.pls`

- ...:

  additional parameters

## Value

variance-covariance matrix

## Details

The function returns the variance-covariance matrix for the optimal
number of components. It can be applied to objects returned by `pls.ic`
and `pls.cv`.

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

[`coef.plsdof`](https://fbertran.github.io/plsdof/reference/coef.plsdof.md),
[`pls.ic`](https://fbertran.github.io/plsdof/reference/pls.ic.md),
[`pls.cv`](https://fbertran.github.io/plsdof/reference/pls.cv.md)

## Author

Nicole Kraemer

## Examples

``` r

n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)


pls.object<-pls.ic(X,y,m=5,criterion="bic")
my.vcov<-vcov(pls.object)
#> WARNING: Covariance of regression coefficients is not available.
#> Returning NULL object.
my.sd<-sqrt(diag(my.vcov)) # standard deviation of regression coefficients
```
