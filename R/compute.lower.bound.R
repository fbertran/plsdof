#' Lower bound for the Degrees of Freedom
#' 
#' This function computes the lower bound for the the Degrees of Freedom of PLS
#' with 1 component.
#' 
#' If the decay of the eigenvalues of \code{cor(X)} is not too fast, we can
#' lower-bound the Degrees of Freedom of PLS with 1 component. Note that we
#' implicitly assume that we use scaled predictor variables to compute the PLS
#' solution.
#' 
#' @param X matrix of predictor observations.
#' @return \item{bound}{logical. bound is \code{TRUE} if the decay of the
#' eigenvalues is slow enough} \item{lower.bound}{if bound is TRUE, this is the
#' lower bound, otherwise, it is set to -1}
#' @author Nicole Kraemer
#' @seealso \code{\link{pls.model}}
#' @references
#' 
#' Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial Least
#' Squares Regression". Journal of the American Statistical Association 106
#' (494) \url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}
#' @keywords math
#' @examples
#' 
#' # Boston Housing data
#' library(MASS)
#' data(Boston)
#' X<-Boston[,-14]
#' my.lower<-compute.lower.bound(X)
#' 
#' @export compute.lower.bound
compute.lower.bound=function(X){
    S=cor(X)
    lower.bound=-1
    lambda<-eigen(S)$values
    bound=FALSE
    if (2*max(lambda)<=sum(lambda)){
        bound=TRUE
        lower.bound<-1+ sum(lambda)/max(lambda)
    } 
    return(list(bound=bound,lower.bound=lower.bound))
}
