#' Comparison of Partial Least Squares Regression, Principal Components
#' Regression and Ridge Regression.
#' 
#' This function computes the test error over several runs for (a) PLS, (b) PCR
#' (c) Ridge Regression and (d) the null model, that is the mean of \code{y}.
#' In the first three cases, the optimal model is selected via
#' cross-validation.
#' 
#' The function computes the test error, the cross-validation-optimal model
#' parameters, their corresponding Degrees of Freedom, and the
#' sum-of-squared-residuals (SSR) for PLS and PCR.
#' 
#' @param X matrix of predictor observations.
#' @param y vector of response observations. The length of \code{y} is the same
#' as the number of rows of \code{X}.
#' @param m maximal number of components for PLS. Default is \code{m=ncol(X)}.
#' @param R number of runs. Default is 20.
#' @param ratio ratio no of training examples/(no of training examples + no of
#' test examples). Default is 0.8
#' @param verbose If \code{TRUE}, the functions plots the progress of the
#' function. Default is \code{TRUE}.
#' @param k number of cross-validation splits. Default is 10.
#' @param nsamples number of data points. Default is \code{nrow(X)}.
#' @param use.kernel Use kernel representation for PLS? Default is
#' \code{use.kernel=FALSE}.
#' @param supervised Should the principal components be sorted by decreasing
#' squared correlation to the response? Default is FALSE.
#' @return \item{MSE}{data frame of size R x 4. It contains the test error for
#' the four different methods for each of the R runs.} \item{M}{data frame of
#' size R x 4. It contains the optimal model parameters for the four different
#' methods for each of the R runs.} \item{DoF}{data frame of size R x 4. It
#' contains the Degrees of Freedom (corresponding to \code{M}) for the four
#' different methods for each of the R runs.} \item{res.pls}{matrix of size R x
#' (ncol(X+1)). It contains the SSR for PLS for each of the R runs.}
#' \item{res.pcr}{matrix of size R x (ncol(X+1)). It contains the SSR for PCR
#' for each of the R runs.} \item{DoF.all}{matrix of size R x (ncol(X+1)). It
#' contains the Degrees of Freedom for PLS for all components for each of the R
#' runs.}
#' @author Nicole Kraemer
#' @seealso \code{\link{pls.cv}}, \code{\link{pcr.cv}},
#' \code{\link{benchmark.pls}}
#' @references
#' 
#' Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial Least
#' Squares Regression". Journal of the American Statistical Association 106
#' (494) \url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}
#' @keywords multivariate
#' @examples
#' 
#' \donttest{
#' # Boston Housing data
#' library(MASS)
#' data(Boston)
#' X<-as.matrix(Boston[,1:4]) # select the first 3 columns as predictor variables
#' y<-as.vector(Boston[,14])
#' 
#' my.benchmark<-benchmark.regression(X,y,ratio=0.5,R=10,k=5)
#' 
#' # boxplot of the mean squared error
#' 
#' boxplot(my.benchmark$MSE,outline=FALSE)
#' 
#' # boxplot of the degrees of freedom, without the null model
#' 
#' boxplot(my.benchmark$DoF[,-4])
#' }
#' 
#' @export benchmark.regression
benchmark.regression=function (X, y, m = ncol(X), R = 20, ratio = 0.8, verbose = TRUE,k = 10, nsamples = nrow(X), use.kernel = FALSE,supervised=FALSE) {
    n <- nsamples
    m.pls <- m.pcr<-lambda.ridge<-vector(length = R) # vector of optimal model parameters
    ntrain <- floor(n * ratio) # number of training observations
    DoF.pls <- DoF.pcr <- DoF.ridge<-vector(length = R) # Degrees of Freedom
    mse.pls <- mse.pcr <- mse.ridge <- vector(length = R) # test error
    mse.null <- vector(length = R) # test error of the null model, i.e. mean(ytrain)
    DoF.all<-res.pls<-res.pcr<-matrix(,R,ncol(X)+1) # residuals, do we need this, do we need DoF.all
    for (i in 1:R) {
        if (verbose == TRUE) {
            cat(paste("iteration no ", i, " \n"))
        }
        # subsample
        samples <- sample(nrow(X), n, replace = FALSE)
        XX <- X[samples, ]
        yy <- y[samples]
        # split into training and test
        train <- sample(n, ntrain, replace = FALSE)
        Xtrain <- XX[train, , drop = FALSE]
        Xtest <- XX[-train, , drop = FALSE]
        ytrain <- yy[train]
        ytest <- yy[-train]
        # null model
        mse.null[i] = mean((mean(ytrain) - ytest)^2)
        # pls
        pls.object<- pls.cv(Xtrain, ytrain, use.kernel = use.kernel, m = m,k=k)
        m.pls[i]<- pls.object$m.opt
        pls.object = pls.model(Xtrain, ytrain, Xtest = Xtest,ytest = ytest, m=m,compute.DoF = TRUE, compute.jacobian = FALSE, use.kernel = use.kernel)
        res.train<-ytrain-rep(1,ntrain)%*%t(pls.object$intercept) - Xtrain%*%pls.object$coefficients
        res.pls[i,]<-apply(res.train^2,2,mean)
        DoF.all[i,]<-pls.object$DoF
        mse.pls[i] <- pls.object$mse[m.pls[i] + 1]
        DoF.pls[i]<-pls.object$DoF[m.pls[i]+1]
        #pcr
        pcr.object<-pcr.cv(Xtrain,ytrain,k=k,m=m,supervised=supervised)
        m.pcr[i]<-pcr.object$m.opt
        DoF.pcr[i]<-m.pcr[i]+1
        pcr.final<-pcr(Xtrain,ytrain)
        res.train<-ytrain - rep(1,ntrain)%*%t(pcr.final$intercept) - Xtrain%*%pcr.final$coefficients
        res.pcr[i,]<-apply(res.train^2,2,mean)
        res<-ytest-pcr.object$intercept - Xtest%*%pcr.object$coefficients
        mse.pcr[i]<-apply(res^2,2,mean)
        # ridge
        ridge.object<-ridge.cv(Xtrain,ytrain,k=k)
        lambda.ridge[i]<-ridge.object$lambda.opt
        res<-ytest-ridge.object$intercept - Xtest%*%ridge.object$coefficients
        mse.ridge[i]<-apply(res^2,2,mean)
        XX<-scale(Xtrain)
        SS<-t(XX)%*%XX
        lambda<-eigen(SS)$values
        DoF.ridge[i]=sum(lambda/(lambda+lambda.ridge[i]))
        }
    namen <- c("PLS", "PCR", "RIDGE","NULL")
    MSE <- matrix(, R, 4)
    MSE[, 1] <- mse.pls
    MSE[, 2] <- mse.pcr
    MSE[, 3] <- mse.ridge
    MSE[,4]<-mse.null
    MSE <- data.frame(MSE)
    colnames(MSE) = namen
    M <- matrix(, R, 4)
    M[, 1] <- m.pls
    M[, 2] <- m.pcr
    M[, 3] <- lambda.ridge
    M[, 4] <- rep(0, R)
    M <- data.frame(M)
    colnames(M) = namen
    DoF <- matrix(, R, 4)
    DoF[, 1] <- DoF.pls
    DoF[, 2] <- DoF.pcr
    DoF[, 3] <- DoF.ridge
    DoF[, 4] <- rep(1, R)
    DoF <- data.frame(DoF)
    colnames(DoF) = namen
    return(list(MSE = MSE, M = M, DoF = DoF,res.pls=res.pls,res.pcr=res.pcr,DoF.all=DoF.all))
}
