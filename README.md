---
title: "Degrees of Freedom and Statistical Inference for Partial Least Squares Regression"
author: "Frédéric Bertrand"
output: github_document
---

[![CRAN status](https://www.r-pkg.org/badges/version/plsdof)](https://cran.r-project.org/package=plsdof)

<!-- README.md is generated from README.Rmd. Please edit that file -->


# plsdof

The plsdof package provides Degrees of Freedom estimates
        for Partial Least Squares (PLS) Regression. Model selection for
        PLS is based on various information criteria (aic, bic, gmdl)
        or on cross-validation. Estimates for the mean and covariance
        of the PLS regression coefficients are available. They allow
        the construction of approximate confidence intervals and the
        application of test procedures.
        Further, cross-validation procedures for Ridge Regression and 
        Principal Components Regression are available.



The plsdof package was completely coded and developped by Nicole Kraemer and Mikio L. Braun. N. Kraemer, M. Sugiyama (2012). It is mainly based on the article "The Degrees of Freedom of Partial Least Squares Regression", *Journal of the American Statistical Association*, **106**(494):697-705, [doi:10.1198/jasa.2011.tm10107](http://dx.doi.org/doi:10.1198/jasa.2011.tm10107).


Yet due to the regular updates in CRAN policies, it was removed from the CRAN and orphaned since the former maintainer had stopped updating the package. The plsdof package is required by several packages of Frédéric Bertrand who was then selected as the new maintainer since late 2018. 


## Installation

You can install the released version of plsdof from [CRAN](https://CRAN.R-project.org) with:


```r
install.packages("plsdof")
```

You can install the development version of plsdof from [github](https://github.com) with:


```r
devtools::install_github("fbertran/plsdof")
```

## Example

###PLS model example.

The pls.model function computes the Partial Least Squares fit.


```r
n<-50 # number of observations
p<-15 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)

ntest<-200 #
Xtest<-matrix(rnorm(ntest*p),ncol=p) # test data
ytest<-rnorm(ntest) # test data

library(plsdof)
# compute PLS + degrees of freedom + prediction on Xtest
first.object<-pls.model(X,y,compute.DoF=TRUE,Xtest=Xtest,ytest=NULL)

# compute PLS + test error
second.object=pls.model(X,y,m=10,Xtest=Xtest,ytest=ytest)
```

###Model selection for Partial Least Squares based on information criteria

The pls.ic function computes the optimal model parameters using one of three different model selection criteria (aic, bic, gmdl) and based on two different Degrees of Freedom estimates for PLS.


```r
n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)

# compute linear PLS
pls.object<-pls.ic(X,y,m=ncol(X))
```
