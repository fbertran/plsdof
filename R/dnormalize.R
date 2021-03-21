#' Derivative of normalization function
#' 
#' This function computes the derivative of the function \deqn{v\mapsto
#' \frac{v}{\|v\|}} with respect to y.
#' 
#' The first derivative of the normalization operator is
#' \deqn{\frac{\partial}{\partial y}\left(v\mapsto
#' \frac{v}{\|v\|}\right)=\frac{1}{\|v\|}\left(I_n - \frac{v v^ \top}{v^\top
#' v}\right) \frac{\partial v}{\partial y}}
#' 
#' @param v vector of length n.
#' @param dv derivative of v with respect to y. As y is a vector of length n,
#' the derivative is a matrix of size nxn.
#' @return the Jacobian matrix of the normalization function. This is a matrix
#' of size nxn.
#' @author Nicole Kraemer, Mikio L. Braun
#' @seealso \code{\link{normalize}}
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
#' v<-rnorm(15)
#' dv<-diag(15)
#' d.object<-dnormalize(v,dv)
#' 
#' @export dnormalize
dnormalize <-
function(v,dv){
    n<-length(v)
    vn<-normalize(v)
    dn=(1/sqrt(sum(v^2)))*(diag(n)- vn%*%t(vn))%*%dv
    return(dn)
}

