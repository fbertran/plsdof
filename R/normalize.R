#' Normalization of vectors
#' 
#' Normalization of vectors.
#' 
#' The vector \code{v} is normalized to length 1. If \code{w} is given, it is
#' normalized by the length of \code{v}.
#' 
#' @param v vector
#' @param w optional vector
#' @return \item{v}{normalized \code{v}} \item{w}{normalized \code{w}}
#' @author Nicole Kraemer, Mikio L. Braun
#' @keywords math
#' @examples
#' 
#' v<-rnorm(5)
#' w<-rnorm(10)
#' dummy<-normalize(v,w)
#' 
#' @export normalize
normalize <-
function(v,w=NULL){
    norm.v<-sqrt(sum(v^2))
    v<-v/norm.v
    if (is.null(w)==TRUE){
        return(v)
    }
    if (is.null(w)==FALSE){
        w<-w/norm.v
        return(list(v=v,w=w))
    }
}

