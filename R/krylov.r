#' Krylov sequence
#' 
#' This function computes the Krylov sequence of a matrix and a vector.
#' 
#' 
#' @param A square matrix of dimension p x p.
#' @param b vector of length p
#' @param m length of the Krylov sequence
#' @return A matrix of size p x m containing the sequence b,Ab,..., A^(m-1)b.
#' @author Nicole Kraemer
#' @keywords math
#' @examples
#' 
#' A<-matrix(rnorm(8*8),ncol=8)
#' b<-rnorm(8)
#' K<-krylov(A,b,4)
#' 
#' @export krylov
krylov<-function(A,b,m){
    K<-matrix(,length(b),m)
    dummy<-b
    for (i in 1:m){
        K[,i]<-dummy
        dummy<-A%*%dummy
    }
    return(K)

}
