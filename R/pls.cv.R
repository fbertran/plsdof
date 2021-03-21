#' Model selection for Partial Least Squares based on cross-validation
#' 
#' This function computes the optimal model parameter using cross-validation.
#' 
#' The data are centered and scaled to unit variance prior to the PLS
#' algorithm. It is possible to estimate the covariance matrix of the
#' cv-optimal regression coefficients (\code{compute.covariance=TRUE}).
#' Currently, this is only implemented if \code{use.kernel=FALSE}.
#' 
#' @param X matrix of predictor observations.
#' @param y vector of response observations. The length of \code{y} is the same
#' as the number of rows of \code{X}.
#' @param k number of cross-validation splits. Default is 10.
#' @param groups an optional vector with the same length as \code{y}. It
#' encodes a partitioning of the data into distinct subgroups. If \code{groups}
#' is provided, \code{k=10} is ignored and instead, cross-validation is
#' performed based on the partioning. Default is \code{NULL}.
#' @param m maximal number of Partial Least Squares components. Default is
#' \code{m=ncol(X)}.
#' @param use.kernel Use kernel representation? Default is
#' \code{use.kernel=FALSE}.
#' @param compute.covariance If \code{TRUE}, the function computes the
#' covariance for the cv-optimal regression coefficients.
#' @param method.cor How should the correlation to the response be computed?
#' Default is ''pearson''.
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
#' the optimal model, based on mean squared error} \item{covariance}{If
#' \code{TRUE} and \code{use.kernel=FALSE}, the covariance of the cv-optimal
#' regression coefficients (based on mean squared error) is returned.}
#' @author Nicole Kraemer, Mikio L. Braun
#' @seealso \code{\link{pls.model}}, \code{\link{pls.ic}}
#' @references
#' 
#' Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial Least
#' Squares Regression". Journal of the American Statistical Association 106
#' (494) \url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}
#' 
#' Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and
#' Efficient Model Selection", Proceedings of the 24th International Conference
#' on Machine Learning, Omni Press, 441 - 448
#' @keywords multivariate
#' @examples
#' 
#' n<-50 # number of observations
#' p<-5 # number of variables
#' X<-matrix(rnorm(n*p),ncol=p)
#' y<-rnorm(n)
#' 
#' # compute linear PLS
#' pls.object<-pls.cv(X,y,m=ncol(X))
#' 
#' # define random partioning
#' groups<-sample(c("a","b","c"),n,replace=TRUE)
#' pls.object1<-pls.cv(X,y,groups=groups)
#' 
#' @export pls.cv
pls.cv<-function (X, y, k = 10, groups = NULL, m = ncol(X), use.kernel = FALSE, 
    compute.covariance = FALSE,method.cor="pearson") 
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
    m = min(m, ntrain.min - 1, p)
    cv.error.matrix = matrix(0, k, m + 1)
    rownames(cv.error.matrix) = my.names
    colnames(cv.error.matrix) = 0:m
    cor.error.matrix<-cv.error.matrix
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain = X[-omit, , drop = FALSE]
        ytrain = y[-omit]
        Xtest = X[omit, , drop = FALSE]
        ytest = y[omit]
        pls.object <- pls.model(Xtrain, ytrain, m = m, Xtest = Xtest, 
            ytest = ytest, compute.DoF = FALSE, use.kernel = use.kernel,method.cor=method.cor)
        cv.error.matrix[i, ] <- pls.object$mse
        cor.error.matrix[i, ] <- pls.object$cor
    }
    cv.error = apply(cv.error.matrix, 2, mean)
    cor.error<-apply(cor.error.matrix,2,mean)
    m.opt <- which.min(cv.error) - 1
    m.opt.cor<-which.max(cor.error) - 1
    if (compute.covariance == TRUE) {
        use.kernel = FALSE
    }
    pls.object <- pls.model(X, y, m = max(m.opt, m.opt.cor,1), use.kernel = use.kernel, 
        compute.DoF = compute.covariance, compute.jacobian = compute.covariance)
    intercept <- pls.object$intercept[m.opt + 1]
    coefficients <- pls.object$coefficients[, m.opt + 1]
    covariance <- pls.object$covariance
    intercept.cor <- pls.object$intercept[m.opt.cor + 1]
    coefficients.cor <- pls.object$coefficients[, m.opt.cor + 1]
    if (compute.covariance == TRUE) {
        #covariancve.cor<-covariance[m.opt.cor + 1, , ]
        covariance <- covariance[m.opt + 1, , ]
    }
    outlist = list(cv.error.matrix = cv.error.matrix, cor.error.matrix=cor.error.matrix,cv.error = cv.error, cor.error=cor.error,
        m.opt = m.opt, m.opt.cor=m.opt.cor,covariance = covariance, intercept = intercept, intercept.cor=intercept.cor,
        coefficients = coefficients,coefficients.cor=coefficients.cor)
    class(outlist) = "plsdof"
    return(outlist)
}


