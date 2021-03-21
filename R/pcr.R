#' Principal Components Regression
#' 
#' This function computes the Principal Components Regression (PCR) fit.
#' 
#' The function first scales all predictor variables to unit variance, and then
#' computes the PCR fit for all components. Is \code{supervised=TRUE}, we sort
#' the principal correlation according to the squared correlation to the
#' response.
#' 
#' @param X matrix of predictor observations.
#' @param y vector of response observations. The length of \code{y} is the same
#' as the number of rows of \code{X}.
#' @param scale Should the predictor variables be scaled to unit variance?
#' Default is \code{TRUE}.
#' @param m maximal number of principal components. Default is
#' \code{m=min(ncol(X),nrow(X)-1)}.
#' @param eps precision. Eigenvalues of the correlation matrix of \code{X} that
#' are smaller than \code{eps} are set to 0. The default value is
#' \code{eps=10^{-6}.}
#' @param supervised Should the principal components be sorted by decreasing
#' squared correlation to the response? Default is FALSE.
#' @return \item{coefficients}{matrix of regression coefficients, including the
#' coefficients of the null model, i.e. the constant model \code{mean(y)}. }
#' \item{intercept}{vector of intercepts, including the intercept of the null
#' model, i.e. the constant model \code{mean(y)}. }
#' @author Nicole Kraemer
#' @seealso \code{\link{pcr.cv}}, \code{\link{pls.cv}}
#' @keywords multivariate
#' @examples
#' 
#' n<-50 # number of observations
#' p<-15 # number of variables
#' X<-matrix(rnorm(n*p),ncol=p)
#' y<-rnorm(n)
#' 
#' my.pcr<-pcr(X,y,m=10)
#' 
#' 
#' @export pcr
pcr<-function (X, y, scale = TRUE, m = min(ncol(X), nrow(X) - 1), 
    eps = 1e-06,supervised=FALSE) 
{
    p <- ncol(X)
    n <- nrow(X)
    Beta <- matrix(, p, m)
    mean.y <- mean(y)
    y <- scale(y, scale = FALSE)
    mean.X <- apply(X, 2, mean)
    if (scale == FALSE) {
        sd.X <- rep(1, p)
    }
    if (scale == TRUE) {
        sd.X <- apply(X, 2, sd)
        sd.X[sd.X == 0] = 1
    }
    X <- X - rep(1, nrow(X)) %*% t(mean.X)
    X <- X/(rep(1, nrow(X)) %*% t(sd.X))
    my.svd = svd(X)
    sigma<- (my.svd$d)[1:m]
    U <- my.svd$v[,1:m]
    V<-my.svd$u[,1:m]
    
    if (supervised==TRUE){
	cor2<-as.vector(cor(V,y))^2
	my.sort=order(cor2,decreasing=TRUE)
	sigma=sigma[my.sort]
	U=U[,my.sort]
	V=V[,my.sort]
    }
    
    sigmainv <- 1/sigma
    sigmainv[sigma < eps] = 0
    my.d<-as.vector(t(V) %*% y)
    #cat(paste("length of my.d: ",length(my.d),"\n"))
    #cat(paste("length of sinv: ",length(sigmainv),"\n"))
    dummy <- sigmainv*my.d
    for (i in 1:m) {
        Beta[, i] <- U[, 1:i, drop = FALSE] %*% dummy[1:i]
    }
    coefficients <- matrix(0, p, m + 1)
    coefficients[, 2:(m + 1)] = Beta/(sd.X %*% t(rep(1, m)))
    intercept <- rep(mean.y, m + 1) - t(coefficients) %*% mean.X
    return(list(intercept = intercept, coefficients = coefficients))
}
