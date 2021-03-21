#' Index of the first local minimum.
#' 
#' This function computes the index of the first local minimum.
#' 
#' 
#' @param x vector.
#' @return the index of the first local minimum of \code{x}.
#' @author Nicole Kraemer
#' @references Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of
#' Partial Least Squares Regression". Journal of the American Statistical
#' Association. ahead of print 106 (494)
#' \url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}
#' @keywords math
#' @examples
#' 
#' v<-rnorm(30)
#' out<-first.local.minimum(v)
#' 
#' @export first.local.minimum
first.local.minimum<-function(x){
    if (length(x)<=2){
        dummy<-which.min(x)
    }
    if (length(x)>2){
        m<-length(x)
        ascending.afterwards<-(x[1:(m-1)]<=x[2:m])
        dummy<-m
        if (sum(ascending.afterwards)>0){
        dummy<-(1:(m-1))[ascending.afterwards]
        dummy<-min(dummy)
        }

    }

return(dummy)
}
