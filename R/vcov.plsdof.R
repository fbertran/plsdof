#' Variance-covariance matrix
#' 
#' This function returns the variance-covariance matrix of a plsdof-object.
#' 
#' The function returns the variance-covariance matrix for the optimal number
#' of components. It can be applied to objects returned by \code{pls.ic} and
#' \code{pls.cv}.
#' 
#' @param object an object of class "plsdof" that is returned by the function
#' \code{linear.pls}
#' @param ... additional parameters
#' @return variance-covariance matrix
#' @author Nicole Kraemer
#' @seealso \code{\link{coef.plsdof}}, \code{\link{pls.ic}},
#' \code{\link{pls.cv}}
#' @references Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of
#' Partial Least Squares Regression". Journal of the American Statistical
#' Association 106 (494)
#' \url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}
#' 
#' Kraemer, N., Sugiyama M., Braun, M.L. (2009) "Lanczos Approximations for the
#' Speedup of Kernel Partial Least Squares Regression." Proceedings of the
#' Twelfth International Conference on Artificial Intelligence and Statistics
#' (AISTATS), p. 272-279
#' @keywords models
#' @examples
#' 
#' 
#' n<-50 # number of observations
#' p<-5 # number of variables
#' X<-matrix(rnorm(n*p),ncol=p)
#' y<-rnorm(n)
#' 
#' 
#' pls.object<-pls.ic(X,y,m=5,criterion="bic")
#' my.vcov<-vcov(pls.object)
#' my.sd<-sqrt(diag(my.vcov)) # standard deviation of regression coefficients
#' 
#' @export
vcov.plsdof=function(object,...){
    dummy<-object$covariance
    if (is.null(dummy)==TRUE){
        cat(paste("WARNING: Covariance of regression coefficients is not available.\n"))
        cat(paste("Returning NULL object.\n"))
    }
    return(dummy)
}
