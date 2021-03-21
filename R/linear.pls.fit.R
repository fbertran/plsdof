#' Linear Partial Least Squares Fit
#' 
#' This function computes the Partial Least Squares solution and the first
#' derivative of the regression coefficients. This implementation scales mostly
#' in the number of variables
#' 
#' We first standardize \code{X} to zero mean and unit variance.
#' 
#' @aliases linear.pls.fit
#' @param X matrix of predictor observations.
#' @param y vector of response observations. The length of \code{y} is the same
#' as the number of rows of \code{X}.
#' @param m maximal number of Partial Least Squares components. Default is
#' \code{m}=ncol(X).
#' @param compute.jacobian Should the first derivative of the regression
#' coefficients be computed as well? Default is \code{FALSE}
#' @param DoF.max upper bound on the Degrees of Freedom. Default is
#' \code{min(ncol(X)+1,nrow(X)-1)}.
#' @return \item{coefficients}{matrix of regression coefficients}
#' \item{intercept}{vector of regression intercepts} \item{DoF}{Degrees of
#' Freedom} \item{sigmahat}{vector of estimated model error} \item{Yhat}{matrix
#' of fitted values} \item{yhat}{vector of squared length of fitted values}
#' \item{RSS}{vector of residual sum of error} \code{covariance}if
#' \code{compute.jacobian} is \code{TRUE}, the function returns the array of
#' covariance matrices for the PLS regression coefficients. \item{TT}{matrix of
#' normalized PLS components}
#' @author Nicole Kraemer
#' @seealso \code{\link{kernel.pls.fit}},
#' \code{\link{pls.cv}},\code{\link{pls.model}}, \code{\link{pls.ic}}
#' @references Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of
#' Partial Least Squares Regression". Journal of the American Statistical
#' Association 106 (494)
#' \url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}
#' @keywords multivariate
#' @examples
#' 
#' n<-50 # number of observations
#' p<-5 # number of variables
#' X<-matrix(rnorm(n*p),ncol=p)
#' y<-rnorm(n)
#' 
#' pls.object<-linear.pls.fit(X,y,m=5,compute.jacobian=TRUE)
#' 
#' @export linear.pls.fit
linear.pls.fit=function (X, y, m = ncol(X),compute.jacobian=FALSE,DoF.max=min(ncol(X)+1,nrow(X)-1)){
    p <- ncol(X)
    n <- nrow(X)
    m<-min(m,DoF.max)
    Beta <- matrix(, p, m) # matrix of regression coefficients
    W <- V <- Beta
    dW<-dBeta<-dV<-NULL
    if (compute.jacobian==TRUE){
    	dW <- dBeta <- dV <- array(dim = c(m, p, n))
    }
    X0<-X
    y0<-y
    # scaling of the data
    mean.y <- mean(y)
    y <- scale(y, scale = FALSE)
    mean.X <- apply(X, 2, mean)
    sd.X <- apply(X, 2, sd)
    sd.X[sd.X == 0] = 1
    X <- X - rep(1, nrow(X)) %*% t(mean.X)
    X <- X/(rep(1, nrow(X)) %*% t(sd.X))
    dcoefficients=NULL
    A <- t(X) %*% X
    b <- t(X) %*% y
    for (i in 1:m) {
        if (i == 1) {
            W[, i] <- b
            if (compute.jacobian==TRUE){
                dW[i, , ] = t(X)
                dW[i, , ] <- dA(W[, i], A, dW[i, , ])
                dV[i, , ] <- dW[i, , ]
            }
            W[, i] <- W[, i]/sqrt((sum((W[, i]) * (A %*% W[,i]))))
            V[, i] <- W[, i]
            Beta[, i] <- sum(V[, i] * b) * V[, i]
            if (compute.jacobian==TRUE){
                dBeta[i, , ] <- dvvtz(V[, i], b, dV[i, , ], t(X))
            }
        }
        if (i > 1) {
            W[, i] <- b - A %*% Beta[, i - 1]
            if (compute.jacobian==TRUE){
                dW[i, , ] <- t(X) - A %*% dBeta[i - 1, , ]
            }
            V[, i] <- W[, i] - vvtz(V[, 1:(i - 1), drop = FALSE],A %*% W[, i])
            if (compute.jacobian==TRUE){
                dV[i, , ] = dW[i, , ] - dvvtz(V[, 1:(i - 1), drop = FALSE],A %*% W[, i], dV[1:(i - 1), , , drop = FALSE],A %*% dW[i, , ])
                dV[i, , ] <- dA(V[, i], A, dV[i, , ])
            }
            V[, i] <- V[, i]/sqrt((sum(t(V[, i]) %*% A %*% V[,i])))
            Beta[, i] = Beta[, i - 1] + sum(V[, i] * b) * V[,i]
            if (compute.jacobian==TRUE){
                dBeta[i, , ] <- dBeta[i - 1, , ] + dvvtz(V[, i], b, dV[i, , ], t(X))
                }
        }
    }
    dcoefficients<-NULL
    if (compute.jacobian==TRUE){
        dcoefficients<-array(0,dim=c(m+1,p,n))
        dcoefficients[2:(m+1),,]=dBeta
    }
    # compute the zero model as well
    sigmahat <- RSS <- yhat <- vector(length=m+1)
    DoF<-1:(m+1)
    Yhat<-matrix(,n,m+1)
    dYhat<-array(dim=c(m+1,n,n))
    coefficients<-matrix(0,p,m+1)
    coefficients[,2:(m+1)] = Beta/(sd.X %*% t(rep(1, m)))
    intercept <- rep(mean.y, m+1) - t(coefficients) %*% mean.X
    covariance<-NULL
    if (compute.jacobian==TRUE){
        covariance<-array(0,dim=c(m+1,p,p))
        DD<-diag(1/sd.X)
    }
    for (i in 1:(m+1)) {
        Yhat[,i]<-X0%*%coefficients[,i] + intercept[i]
        res<-y0-Yhat[,i]
        yhat[i] <- sum((Yhat[,i])^2)
        RSS[i] <- sum(res^2)
        if (compute.jacobian==TRUE){
        dYhat[i,,] <- X %*% dcoefficients[i, , ] + matrix(1,n,n)/n
        DoF[i] <- sum(diag(dYhat[i,,]))
        dummy <- (diag(n) - dYhat[i,,]) %*% (diag(n) - t(dYhat[i,,]))
        sigmahat[i] <- sqrt(RSS[i]/sum(diag(dummy)))
        if (i>1){
        covariance[i,,]<-sigmahat[i]^2*DD%*%dcoefficients[i,,]%*%t(dcoefficients[i,,])%*%DD
        }
        }
    }
    if (compute.jacobian==FALSE){
        sigmahat<-sqrt(RSS/(n-DoF))
    }
    TT=X%*%V
    DoF[DoF>DoF.max]=DoF.max
     intercept<-as.vector(intercept)
    outlist = list(Yhat=Yhat,yhat=yhat,DoF=DoF,coefficients = coefficients, intercept=intercept, RSS=RSS,sigmahat = sigmahat,TT=TT,covariance=covariance)
    return(outlist)
}
