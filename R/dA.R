#' Derivative of normalization function
#' 
#' This function computes the derivative of the function \deqn{v\mapsto
#' \frac{w}{\|w\|_A}} with respect to y.
#' 
#' The first derivative of the normalization operator is
#' \deqn{\frac{\partial}{\partial y}\left(w\mapsto
#' \frac{w}{\|w\|_A}\right)=\frac{1}{\|w\|}\left(I_n - \frac{w w^ \top
#' A}{w^\top w}\right) \frac{\partial w}{\partial y}}
#' 
#' @param w vector of length n.
#' @param A square matrix that defines the norm
#' @param dw derivative of w with respect to y. As y is a vector of length n,
#' the derivative is a matrix of size nxn.
#' @return the Jacobian matrix of the normalization function. This is a matrix
#' of size nxn.
#' @author Nicole Kraemer
#' @seealso \code{\link{normalize}}, \code{\link{dnormalize}}
#' @references Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of
#' Partial Least Squares Regression". Journal of the American Statistical
#' Association 106 (494)
#' \url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}
#' 
#' Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and
#' Efficient Model Selection", Proceedings of the 24th International Conference
#' on Machine Learning, Omni Press, 441 - 448
#' @keywords math
#' @examples
#' 
#' w<-rnorm(15)
#' dw<-diag(15)
#' A<-diag(1:15)
#' d.object<-dA(w,A,dw)
#' 
#' @export dA
dA<-function(w,A,dw){
wa<-sqrt(sum((w*(A%*%w))))
dummy<-(1/wa)*(diag(length(w))- w%*%t(w)%*%A/(wa^2))%*%dw
return(dummy)
}
