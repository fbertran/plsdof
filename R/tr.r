#' Trace of a matrix
#' 
#' This function computes the trace of a matrix.
#' 
#' 
#' @param M square matrix
#' @return The trace of the matrix M.
#' @author Nicole Kraemer
#' @keywords math
#' @examples
#' 
#' M<-matrix(rnorm(8*8),ncol=8)
#' tr.M<-tr(M)
#' 
#' @export tr
tr<-function(M){
return(sum(diag(M)))
}
