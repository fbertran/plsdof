#' Computation of the Degrees of Freedom
#' 
#' This function computes the Degrees of Freedom using the Krylov
#' representation of PLS.
#' 
#' This computation of the Degrees of Freedom is based on the equivalence of
#' PLS regression and the projection of the response vector \code{y} onto the
#' Krylov space spanned by \deqn{Ky,K^2 y,...,K^m y.} Details can be found in
#' Kraemer and Sugiyama (2011).
#' 
#' @param pls.object object returned by \code{linear.pls.fit} or by
#' \code{kernel.pls.fit}
#' @param n number of observations
#' @param y vector of response observations.
#' @param K kernel matrix X X^t.
#' @param m number of components
#' @param DoF.max upper bound on the Degrees of Freedom.
#' @return \item{coefficients}{matrix of regression coefficients}
#' \item{intercept}{vector of regression intercepts} \item{DoF}{Degrees of
#' Freedom} \item{sigmahat}{vector of estimated model error} \item{Yhat}{matrix
#' of fitted values} \item{yhat}{vector of squared length of fitted values}
#' \item{RSS}{vector of residual sum of error} \item{TT}{matrix of normalized
#' PLS components}
#' @author Nicole Kraemer, Mikio L. Braun
#' @seealso \code{\link{pls.model}}, \code{\link{pls.ic}}
#' @references
#' 
#' Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial Least
#' Squares Regression". Journal of the American Statistical Association 106
#' (494) \url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}
#' 
#' Kraemer, N., Sugiyama M., Braun, M.L. (2009) "Lanczos Approximations for the
#' Speedup of Kernel Partial Least Squares Regression." Proceedings of the
#' Twelfth International Conference on Artificial Intelligence and Statistics
#' (AISTATS), p. 272-279
#' @keywords multivariate
#' @examples
#' 
#' # this is an internal function
#' 
#' @export pls.dof
pls.dof=function(pls.object,n,y,K,m,DoF.max){
    TT<-pls.object$TT
    Yhat<-pls.object$Yhat[,2:(m+1)]
    TK=matrix(,m,m)
    KY<-krylov(K,K%*%y,m)
    #KK<-diag(n)
    lambda<-eigen(K)$values
    tr.K<-vector(length=m)
    for (i in 1:m){
        #KK<-K%*%KK
        tr.K[i]<-sum(lambda^i)
        #tr.K[i]<-tr(KK)
    }
    BB=t(TT)%*%KY
    BB[row(BB)>col(BB)]=0
    b<-t(TT)%*%y
    DoF=vector(length=m)
    Binv<-backsolve(BB,diag(m))
    tkt<-rep(0,m)
    ykv<-rep(0,m)
    KjT<-array(dim=c(m,n,m))
    dummy<-TT
    for (i in 1:m){
        dummy<-K%*%dummy
        KjT[i,,]<-dummy
    }
    trace.term=rep(0,m)
    for (i in 1:m){
            Binvi<-Binv[1:i,1:i,drop=FALSE]
            ci<-Binvi%*%b[1:i]
            Vi<-TT[,1:i,drop=FALSE]%*%t(Binvi)
            trace.term[i]<-sum(ci*tr.K[1:i])
            ri<-y-Yhat[,i]
            for (j in 1:i){
                KjTj=KjT[j,,]
                tkt[i]<-tkt[i]+ci[j]*tr(t(TT[,1:i,drop=FALSE])%*%KjTj[,1:i,drop=FALSE])
                ri<-K%*%ri
                ykv[i]<-ykv[i]+ sum(ri*Vi[,j])
            }
            }
    DoF<-trace.term + 1:m - tkt + ykv
    DoF[DoF>DoF.max]=DoF.max
    sigmahat=sqrt(pls.object$RSS[-1]/(n-DoF))
        return(list(DoF=DoF,sigmahat=sigmahat))





}
