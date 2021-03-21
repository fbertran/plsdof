#' Model selection for Partial Least Squares based on information criteria
#' 
#' This function computes the optimal model parameters using one of three
#' different model selection criteria (aic, bic, gmdl) and based on two
#' different Degrees of Freedom estimates for PLS.
#' 
#' There are two options to estimate the Degrees of Freedom of PLS:
#' \code{naive=TRUE} defines the Degrees of Freedom as the number of components
#' +1, and \code{naive=FALSE} uses the generalized notion of Degrees of
#' Freedom. If \code{compute.jacobian=TRUE}, the function uses the Lanczos
#' decomposition to derive the Degrees of Freedom, otherwise, it uses the
#' Krylov representation. (See Kraemer and Sugiyama (2011) for details.) The
#' latter two methods only differ with respect to the estimation of the noise
#' level.
#' 
#' @param X matrix of predictor observations.
#' @param y vector of response observations. The length of \code{y} is the same
#' as the number of rows of \code{X}.
#' @param m maximal number of Partial Least Squares components. Default is
#' \code{m=ncol(X)}.
#' @param criterion Choice of the model selection criterion. One of the three
#' options aic, bic, gmdl.
#' @param naive Use the naive estimate for the Degrees of Freedom? Default is
#' \code{FALSE}.
#' @param use.kernel Use kernel representation? Default is
#' \code{use.kernel=FALSE}.
#' @param compute.jacobian Should the first derivative of the regression
#' coefficients be computed as well? Default is \code{FALSE}
#' @param verbose If \code{TRUE}, the function prints a warning if the
#' algorithms produce negative Degrees of Freedom. Default is \code{TRUE}.
#' @return The function returns an object of class "plsdof". \item{DoF}{Degrees
#' of Freedom} \item{m.opt}{optimal number of components}
#' \item{sigmahat}{vector of estimated model errors}
#' \item{intercept}{intercept} \item{coefficients}{vector of regression
#' coefficients} \item{covariance}{if \code{compute.jacobian=TRUE} and
#' \code{use.kernel=FALSE}, the function returns the covariance matrix of the
#' optimal regression coefficients.} \item{m.crash}{the number of components
#' for which the algorithm returns negative Degrees of Freedom}
#' @author Nicole Kraemer, Mikio L. Braun
#' @seealso \code{\link{pls.model}}, \code{\link{pls.cv}}
#' @references Akaikie, H. (1973) "Information Theory and an Extension of the
#' Maximum Likelihood Principle". Second International Symposium on Information
#' Theory, 267 - 281.
#' 
#' Hansen, M., Yu, B. (2001). "Model Selection and Minimum Descripion Length
#' Principle". Journal of the American Statistical Association, 96, 746 - 774
#' 
#' Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial Least
#' Squares Regression". Journal of the American Statistical Association 106
#' (494) \url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}
#' 
#' Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and
#' Efficient Model Selection", Proceedings of the 24th International Conference
#' on Machine Learning, Omni Press, 441 - 448
#' 
#' Schwartz, G. (1979) "Estimating the Dimension of a Model" Annals of
#' Statistics 26(5), 1651 - 1686.
#' @keywords multivariate
#' @examples
#' 
#' n<-50 # number of observations
#' p<-5 # number of variables
#' X<-matrix(rnorm(n*p),ncol=p)
#' y<-rnorm(n)
#' 
#' # compute linear PLS
#' pls.object<-pls.ic(X,y,m=ncol(X))
#' 
#' @export pls.ic
pls.ic=function (X, y, m = min(ncol(X),nrow(X)-1),criterion="bic",naive=FALSE,use.kernel=FALSE,compute.jacobian=FALSE,verbose=TRUE) {
    m.crash=NA
    n <- nrow(X)
    DoF.max = min(n-1, ncol(X)+1)
    compute.DoF=TRUE
    if (naive==TRUE){
        compute.DoF=FALSE
        compute.jacobian=FALSE
    }
    pls.object <- pls.model(X,y,m,compute.DoF=compute.DoF,compute.jacobian=compute.jacobian,use.kernel=use.kernel)
    RSS <- pls.object$RSS
    yhat <- pls.object$yhat
    sigmahat <- (pls.object$sigmahat)
        DoF <- pls.object$DoF
    if (min(DoF)<=0){
        sign.DoF<-sign(DoF)
        sign.DoF[sign.DoF==0]=-1
        dummy<-(1:(m+1))[sign.DoF==-1]
        mini<-min(dummy)-1
        m<-mini-1
        m.crash<-m+1
        if (verbose==TRUE){
        if (compute.jacobian==TRUE){
        cat(paste("Negative DoF for jacobian. Setting maximal number of components to ",m,".\n"))
        }
        if (compute.jacobian==FALSE){
        cat(paste("Negative DoF. Setting maximal number of components to ",m,".\n"))
        }
        }
        #cat(paste("DoF for ",m," components: ",DoF[1,m+1],"\n"))
        DoF[(m:length(DoF))]=Inf
    }
    if (min(DoF)>0){
    ic <- information.criteria(RSS, DoF, yhat = yhat, sigmahat = sigmahat, 
        n, criterion=criterion)
    #DoF = ic$DoF
    m.opt <- ic$par-1
        score<-ic$score
}
coefficients<-pls.object$coefficients[,m.opt]
intercept<-pls.object$intercept[m.opt]
covariance<-pls.object$covariance
if (compute.jacobian==TRUE){
    covariance=covariance[m.opt+1,,]
}
outlist=list(DoF = DoF, m.opt = m.opt,sigmahat=sigmahat,m.crash=m.crash,intercept=intercept,coefficients=coefficients,covariance=covariance)
class(outlist)="plsdof"    
    return(outlist)
}
