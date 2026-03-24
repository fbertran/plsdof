# Regression coefficients

This function returns the regression coefficients of a plsdof-object.

## Usage

``` r
# S3 method for class 'plsdof'
coef(object, ...)
```

## Arguments

- object:

  an object of class "plsdof" that is returned by the functions `pls.ic`
  and `pls.cv`.

- ...:

  additional parameters

## Value

regression coefficients.

## Details

The function returns the regression coefficients (without intercept) for
the optimal number of components.

## References

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial
Least Squares Regression". Journal of the American Statistical
Association 106 (494)
<https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107>

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom,
and Efficient Model Selection", Proceedings of the 24th International
Conference on Machine Learning, Omni Press, 441 - 448

## See also

[`vcov.plsdof`](https://fbertran.github.io/plsdof/reference/vcov.plsdof.md),
[`pls.model`](https://fbertran.github.io/plsdof/reference/pls.model.md),
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


pls.object<-pls.ic(X,y,criterion="bic")
mycoef<-coef(pls.object)
```
