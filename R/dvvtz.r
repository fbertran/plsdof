#' First derivative of the projection operator
#' 
#' This function computes the first derivative of the projection operator
#' \deqn{P_V z= V V^\top z}
#' 
#' For the computation of the first derivative, we assume that the columns of
#' \code{v} are normalized and mutually orthogonal. (Note that the function
#' will not return an error message if these assumptionsa are not fulfilled. If
#' we denote the columns of \code{v} by \eqn{v_1,\ldots,v_l}, the first
#' derivative of the projection operator is \deqn{ \frac{\partial P}{\partial
#' y}=\sum_{j=1} ^ l \left[ \left(v_j z^ \top + v_j^ \top z I_n
#' \right)\frac{\partial v_j}{\partial y} + v_j v_j ^ \top \frac{\partial
#' z}{\partial y}\right] } Here, n denotes the length of the vectors \eqn{v_j}.
#' 
#' @param v orthonormal basis of the space on which \code{z} is projected.
#' \code{v} is either a matrix or a vector.
#' @param z vector that is projected onto the columns of \code{v}
#' @param dv first derivative of the the columns of \code{v} with respect to a
#' vector y. If \code{v} is a matrix, \code{dv} is an array of dimension
#' \code{ncol(v)}x\code{nrow(v)}x\code{length(y)}. If \code{v} is a vector,
#' \code{dv} is a matrix of dimension \code{nrow(v)}x\code{length(y)}.
#' @param dz first derivative of \code{z} with respect to a vector y. This is a
#' matrix of dimension \code{nrow(v)}x\code{length(y)}.
#' @return The first derivative of the projection operator with respect to y.
#' This is a matrix of dimension \code{nrow(v)}x\code{length(y)}.
#' @note This is an internal function.
#' @author Nicole Kraemer, Mikio L. Braun
#' @seealso \code{\link{vvtz}}
#' @references Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of
#' Partial Least Squares Regression". Journal of the American Statistical
#' Association. 106 (494)
#' \url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}
#' 
#' Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and
#' Efficient Model Selection", Proceedings of the 24th International Conference
#' on Machine Learning, Omni Press, 441 - 448
#' @keywords math
#' @export dvvtz
dvvtz<-function (v, z, dv, dz) 
{
    if (is.matrix(v) == FALSE) {
        v <- matrix(v, ncol = 1)
        dv <- array(dv, dim = c(1, nrow(dv), ncol(dv)))
    }
    k = ncol(v)
    p <- nrow(v)
    n<-dim(dv)[3]
    dummy <- matrix(0,dim(dv)[2],dim(dv)[3])
    for (i in 1:k) {
        D <- (v[, i] %*% t(z) + sum(v[, i] * z) * diag(p)) %*% 
            dv[i, , ] + v[, i] %*% t(v[, i]) %*% dz
        dummy <- dummy + D
    }
    return(dummy)
}
