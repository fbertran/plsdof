#' Partial Least Squares
#' 
#' This function computes the Partial Least Squares fit.
#' 
#' This function computes the Partial Least Squares fit and its Degrees of
#' Freedom. Further, it returns the regression coefficients and various
#' quantities that are needed for model selection in combination with
#' \code{information.criteria}.
#' 
#' @param X matrix of predictor observations.
#' @param y vector of response observations. The length of \code{y} is the same
#' as the number of rows of \code{X}.
#' @param m maximal number of Partial Least Squares components. Default is
#' \code{m=min(ncol(X),nrow(X)-1)}.
#' @param Xtest optional matrix of test observations. Default is
#' \code{Xtest=NULL}.
#' @param ytest optional vector of test observations. Default is
#' \code{ytest=NULL}.
#' @param compute.DoF Logical variable. If \code{compute.DoF=TRUE}, the Degrees
#' of Freedom of Partial Least Squares are computed. Default is
#' \code{compute.DoF=FALSE}.
#' @param compute.jacobian Should the first derivative of the regression
#' coefficients be computed as well? Default is \code{FALSE}
#' @param use.kernel Should the kernel representation be used to compute the
#' solution. Default is \code{FALSE}.
#' @param method.cor How should the correlation to the response be computed?
#' Default is ''pearson''.
#' @return \item{coefficients}{matrix of regression coefficients}
#' \item{intercept}{vector of intercepts} \item{DoF}{vector of Degrees of
#' Freedom} \item{RSS}{vector of residual sum of error} \item{sigmahat}{vector
#' of estimated model error} \item{Yhat}{matrix of fitted values}
#' \item{yhat}{vector of squared length of fitted values} \item{covariance}{if
#' \code{compute.jacobian} is \code{TRUE}, the function returns the array of
#' covariance matrices for the PLS regression coefficients.}
#' \code{prediction}if \code{Xtest} is provided, the predicted y-values for
#' \code{Xtest}. \code{mse}if \code{Xtest} and \code{ytest} are provided, the
#' mean squared error on the test data.  \code{cor}if \code{Xtest} and
#' \code{ytest} are provided, the correlation to the response on the test data.
#' @author Nicole Kraemer, Mikio L. Braun
#' @seealso \code{\link{pls.ic}}, \code{\link{pls.cv}}
#' @references Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of
#' Partial Least Squares Regression". Journal of the American Statistical
#' Association 106 (494)
#' \url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}
#' 
#' Kraemer, N., Sugiyama, M., Braun, M.L. (2009) "Lanczos Approximations for
#' the Speedup of Partial Least Squares Regression", Proceedings of the 12th
#' International Conference on Artificial Intelligence and Stastistics, 272 -
#' 279
#' @keywords multivariate
#' @examples
#' 
#' n<-50 # number of observations
#' p<-15 # number of variables
#' X<-matrix(rnorm(n*p),ncol=p)
#' y<-rnorm(n)
#' 
#' ntest<-200 #
#' Xtest<-matrix(rnorm(ntest*p),ncol=p) # test data
#' ytest<-rnorm(ntest) # test data
#' 
#' # compute PLS + degrees of freedom + prediction on Xtest
#' first.object<-pls.model(X,y,compute.DoF=TRUE,Xtest=Xtest,ytest=NULL)
#' 
#' # compute PLS + test error
#' second.object=pls.model(X,y,m=10,Xtest=Xtest,ytest=ytest)
#' 
#' @export pls.model
pls.model<-function (X, y, m = ncol(X), Xtest = NULL, ytest = NULL, compute.DoF = FALSE, 
    compute.jacobian = FALSE, use.kernel = FALSE,method.cor="pearson") 
{
    if (compute.DoF == FALSE) {
        compute.jacobian == FALSE
    }
    n <- nrow(X)
    p <- ncol(X)
    X0 <- X
    DoF.max = min(n - 1, p + 1)
    m <- min(m, DoF.max)
    DoF <- NULL
    mse <- cor<-NULL
    prediction <- NULL
    coefficients <- NULL
    sigmahat <- NULL
    RSS <- NULL
    intercept <- NULL
    if (compute.jacobian == TRUE) {
        if (use.kernel == TRUE) {
            pls.object <- kernel.pls.fit(X, y, m, compute.jacobian = compute.jacobian, 
                DoF.max = DoF.max)
        }
        if (use.kernel == FALSE) {
            pls.object <- linear.pls.fit(X, y, m, compute.jacobian = TRUE, 
                DoF.max = DoF.max)
        }
        coefficients <- pls.object$coefficients
        intercept <- pls.object$intercept
        Yhat <- pls.object$Yhat
        yhat <- pls.object$yhat
        RSS <- pls.object$RSS
        sigmahat <- pls.object$sigmahat
        DoF <- pls.object$DoF
    }
    if (compute.jacobian == FALSE) {
        if (use.kernel == TRUE) {
            pls.object <- kernel.pls.fit(X, y, m, compute.jacobian = FALSE, 
                DoF.max = DoF.max)
            intercept <- pls.object$intercept
            coefficients <- pls.object$coefficients
            Yhat <- pls.object$Yhat
            yhat <- pls.object$yhat
            RSS <- pls.object$RSS
            sigmahat <- pls.object$sigmahat
            DoF <- pls.object$DoF
            if (compute.DoF == TRUE) {
                mean.X <- apply(X, 2, mean)
                sd.X <- apply(X, 2, sd)
                sd.X[sd.X == 0] = 1
                X <- X - rep(1, nrow(X)) %*% t(mean.X)
                X <- X/(rep(1, nrow(X)) %*% t(sd.X))
                K <- X %*% t(X)
                dof.object = pls.dof(pls.object, K = K, y = y, 
                  n = n, m = m, DoF.max = DoF.max - 1)
                DoF = c(0, dof.object$DoF) + 1
                sigmahat = c(sqrt(RSS[1]/(n - 1)), dof.object$sigmahat)
            }
        }
        if (use.kernel == FALSE) {
            pls.object <- linear.pls.fit(X, y, m, compute.jacobian = FALSE, 
                DoF.max = DoF.max)
            sigmahat <- pls.object$sigmahat
            DoF <- pls.object$DoF
            Yhat <- pls.object$Yhat
            yhat <- pls.object$yhat
            RSS <- pls.object$RSS
            if (compute.DoF == TRUE) {
                mean.X <- apply(X, 2, mean)
                sd.X <- apply(X, 2, sd)
                sd.X[sd.X == 0] = 1
                X <- X - rep(1, nrow(X)) %*% t(mean.X)
                X <- X/(rep(1, nrow(X)) %*% t(sd.X))
                K <- X %*% t(X)
                dof.object = pls.dof(pls.object, K = K, y = y, 
                  n = n, m = m, DoF.max = DoF.max - 1)
                DoF = c(0, dof.object$DoF) + 1
                sigmahat = c(sqrt(RSS[1]/(n - 1)), dof.object$sigmahat)
            }
            coefficients <- pls.object$coefficients
            intercept <- pls.object$intercept
        }
    }
    if (is.null(Xtest) == FALSE) {
        prediction = rep(1, nrow(Xtest)) %*% t(intercept) + Xtest %*% 
            coefficients
        if (is.null(ytest) == FALSE) {
            res <- matrix(, nrow(Xtest), m + 1)
            for (l in 1:(m + 1)) {
                res[, l] = ytest - prediction[, l]
                if (l>1){
                cor[l]<-cor(ytest,prediction[,l],method=method.cor)
                }
            }
            mse = apply(res^2, 2, mean)
        }
    }
    if (compute.DoF == FALSE) {
        DoF = 1:(m + 1)
        sigmahat = sqrt(RSS/(n - DoF))
    }
    covariance = pls.object$covariance
    return(list(prediction = prediction, mse = mse, cor=cor,coefficients = coefficients, 
        intercept = intercept, DoF = DoF, RSS = RSS, Yhat = Yhat, 
        sigmahat = sigmahat, yhat = yhat, covariance = covariance))
}
