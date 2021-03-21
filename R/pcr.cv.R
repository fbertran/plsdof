#' Model selection for Princinpal Components regression based on
#' cross-validation
#' 
#' This function computes the optimal model parameter using cross-validation.
#' Mdel selection is based on mean squared error and correlation to the
#' response, respectively.
#' 
#' The function computes the principal components on the scaled predictors.
#' Based on the regression coefficients \code{coefficients.jackknife} computed
#' on the cross-validation splits, we can estimate their mean and their
#' variance using the jackknife. We remark that under a fixed design and the
#' assumption of normally distributed \code{y}-values, we can also derive the
#' true distribution of the regression coefficients.
#' 
#' @param X matrix of predictor observations.
#' @param y vector of response observations. The length of \code{y} is the same
#' as the number of rows of \code{X}.
#' @param k number of cross-validation splits. Default is 10.
#' @param m maximal number of principal components. Default is
#' \code{m=min(ncol(X),nrow(X)-1)}.
#' @param groups an optional vector with the same length as \code{y}. It
#' encodes a partitioning of the data into distinct subgroups. If \code{groups}
#' is provided, \code{k=10} is ignored and instead, cross-validation is
#' performed based on the partioning. Default is \code{NULL}.
#' @param scale Should the predictor variables be scaled to unit variance?
#' Default is \code{TRUE}.
#' @param eps precision. Eigenvalues of the correlation matrix of \code{X} that
#' are smaller than \code{eps} are set to 0. The default value is
#' \code{eps=10^{-6}.}
#' @param plot.it Logical. If \code{TRUE}, the function plots the
#' cross-validation-error as a function of the number of components. Default is
#' \code{FALSE}.
#' @param compute.jackknife Logical. If \code{TRUE}, the regression
#' coefficients on each of the cross-validation splits is stored. Default is
#' \code{TRUE}.
#' @param method.cor How should the correlation to the response be computed?
#' Default is ''pearson''.
#' @param supervised Should the principal components be sorted by decreasing
#' squared correlation to the response? Default is FALSE.
#' @return \item{cv.error.matrix}{matrix of cross-validated errors based on
#' mean squared error. A row corresponds to one cross-validation split.}
#' \item{cv.error}{vector of cross-validated errors based on mean squared
#' error} \item{m.opt}{optimal number of components based on mean squared
#' error} \item{intercept}{intercept of the optimal model, based on mean
#' squared error} \item{coefficients}{vector of regression coefficients of the
#' optimal model, based on mean squared error} \item{cor.error.matrix}{matrix
#' of cross-validated errors based on correlation. A row corresponds to one
#' cross-validation split.} \item{cor.error}{vector of cross-validated errors
#' based on correlation} \item{m.opt.cor}{optimal number of components based on
#' correlation} \item{intercept.cor}{intercept of the optimal model, based on
#' correlation} \item{coefficients.cor}{vector of regression coefficients of
#' the optimal model, based on correlation}
#' 
#' \item{coefficients.jackknife}{Array of the regression coefficients on each
#' of the cross-validation splits, if \code{compute.jackknife=TRUE}. In this
#' case, the dimension is \code{ncol(X) x (m+1) x k}.}
#' @author Nicole Kraemer, Mikio L. Braun
#' @seealso \code{\link{pls.model}}, \code{\link{pls.ic}}
#' @keywords multivariate
#' @examples
#' 
#' n<-500 # number of observations
#' p<-5 # number of variables
#' X<-matrix(rnorm(n*p),ncol=p)
#' y<-rnorm(n)
#' 
#' # compute PCR 
#' pcr.object<-pcr.cv(X,y,scale=FALSE,m=3)
#' pcr.object1<-pcr.cv(X,y,groups=sample(c(1,2,3),n,replace=TRUE),m=3)
#' 
#' @export pcr.cv
 pcr.cv<-function (X, y, k = 10, m = min(ncol(X), nrow(X) - 1), groups = NULL, 
    scale = TRUE, eps = 1e-06, plot.it = FALSE, compute.jackknife = TRUE,method.cor="pearson",supervised=FALSE) 
{
    n <- nrow(X)
    p <- ncol(X)
    if (is.null(groups) == FALSE) {
        f = as.factor(groups)
        k = length(levels(f))
        my.names = levels(f)
	all.folds <- split(1:n, f)
    }
    if (is.null(groups) == TRUE) {
        f <- rep(1:k, length = n)
        my.names <- 1:k
	all.folds <- split(sample(1:n), f)
    }
    ntrain = vector(length = k)
    for (i in 1:k) {
        ntrain[i] = n - length(all.folds[[i]])
    }
    ntrain.min = min(ntrain)
    m = min(m, ntrain.min - 1)
    cv.error.matrix = matrix(0, k, m + 1)
    rownames(cv.error.matrix) = my.names
    colnames(cv.error.matrix) = 0:m
    cor.error.matrix<-cv.error.matrix
    coefficients.jackknife = NULL
    if (compute.jackknife == TRUE) {
        coefficients.jackknife <- array(dim = c(p, m + 1, k))
    }
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain = X[-omit, , drop = FALSE]
        ytrain = y[-omit]
        Xtest = X[omit, , drop = FALSE]
        ytest = y[omit]
        pcr.object <- pcr(Xtrain, ytrain, scale, m, eps,supervised)
        res <- matrix(, length(ytest), m + 1)
        if (compute.jackknife == TRUE) {
            coefficients.jackknife[, , i] = pcr.object$coefficients
        }
        for (j in 1:(m + 1)) {
            ytest.j<-pcr.object$intercept[j]+ Xtest%*%pcr.object$coefficients[, j]
            res[, j] <- ytest - ytest.j
            if (j>1){
            cor.error.matrix[i,j]<-cor(ytest,ytest.j,method=method.cor)
            }
        }
        cv.error.matrix[i, ] <- apply(res^2, 2, mean)
    }
 #   cat(paste("cv complete \n"))
    cv.error <- apply(cv.error.matrix, 2, mean)
    cor.error <- apply(cor.error.matrix, 2, mean)
    m.opt <- which.min(cv.error) - 1
 #  cat(paste("mse ",m.opt,"\n"))
    m.opt.cor <- which.max(cor.error) - 1
#	cat(paste("mse ",m.opt,"\n"))
    if (plot.it == TRUE) {
        plot(0:m, cv.error, type = "l")
    }
    m.min<-max(2,m.opt,m.opt.cor)
    pcr.object <- pcr(X, y, scale=scale,m=m.min, eps = eps,supervised=supervised)
    coefficients <- pcr.object$coefficients[, m.opt + 1]
    intercept.cor <- pcr.object$intercept[m.opt.cor + 1]
    coefficients.cor <- pcr.object$coefficients[, m.opt.cor + 1]
    intercept <- pcr.object$intercept[m.opt + 1]
    return(list(intercept = intercept, intercept.cor=intercept.cor,coefficients = coefficients, coefficients.cor=coefficients.cor,
        m.opt = m.opt, m.opt.cor=m.opt.cor,cv.error.matrix = cv.error.matrix, cor.error.matrix=cor.error.matrix,cv.error = cv.error, cor.error=cor.error,
        coefficients.jackknife = coefficients.jackknife))
}

