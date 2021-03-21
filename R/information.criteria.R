#' Information criteria
#' 
#' This function computes the optimal model parameters using three different
#' model selection criteria (aic, bic, gmdl).
#' 
#' The Akaike information criterion (aic) is defined as \deqn{{aic}=
#' \frac{{RSS}}{n} + 2\frac{{DoF}}{n} \sigma^ 2\,.} The Bayesian information
#' criterion (bic) is defined as \deqn{{bic}= \frac{{RSS}}{n} +
#' log(n)\frac{{DoF}}{n} \sigma^ 2\,.} The generalized minimum description
#' length (gmdl) is defined as
#' \deqn{gmdl=\frac{n}{2}log(S)+\frac{DoF}{2}log(F)+\frac{1}{2}log(n)} with
#' \deqn{S=\hat \sigma ^2} Note that it is also possible to use the function
#' \code{information.criteria} for other regression methods than Partial Least
#' Squares.
#' 
#' @param RSS vector of residual sum of squares.
#' @param DoF vector of Degrees of Freedom. The length of \code{DoF} is the
#' same as the length of \code{RSS}.
#' @param yhat vector of squared norm of yhat. The length of \code{yhat} is the
#' same as the length of \code{RSS}. It is only needed for gmdl. Default value
#' is \code{NULL}.
#' @param sigmahat Estimated model error. The length of \code{sigmahat} is the
#' same as the length of \code{RSS}.
#' @param n number of observations.
#' @param criterion one of the options "aic", "bic" and "gmdl".
#' @return \item{DoF}{degrees of freedom} \item{score}{vector of the model
#' selection criterion} \item{par}{index of the first local minimum of
#' \code{score}}
#' @author Nicole Kraemer, Mikio Braun
#' @seealso \code{\link{pls.ic}}
#' @references Akaikie, H. (1973) "Information Theory and an Extension of the
#' Maximum Likelihood Principle". Second International Symposium on Information
#' Theory, 267 - 281.
#' 
#' Hansen, M., Yu, B. (2001). "Model Selection and Minimum Descripion Length
#' Principle". Journal of the American Statistical Association, 96, 746 - 774
#' 
#' Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial Least
#' Squares Regression". Journal of the American Statistical Association 106
#' (494) \url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}
#' 
#' Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and
#' Efficient Model Selection", Proceedings of the 24th International Conference
#' on Machine Learning, Omni Press, 441 - 448
#' 
#' Schwartz, G. (1979) "Estimating the Dimension of a Model" Annals of
#' Statistics 26(5), 1651 - 1686.
#' @keywords model
#' @examples
#' 
#' ## This is an internal function called by pls.ic
#' 
#' 
#' @export information.criteria
information.criteria=function (RSS, DoF, yhat=NULL, sigmahat, n,criterion="bic"){
    if (criterion=="aic"){
        score <- as.vector(RSS/n + 2 * (DoF/n) * sigmahat^2)
    }
    if (criterion=="bic"){
        score <- as.vector(RSS/n + log(n) * (DoF/n) * sigmahat^2)
}
        
    if (criterion=="gmdl"){
       SS<-sigmahat^2
        denominator<-DoF*SS
        FF <- (yhat)/(DoF * SS)
        FF[1,FF==0]=Inf
        score <- as.vector((n/2) * log(SS) + (DoF/2) * log(FF) + (1/2) * log(n))
}
    par<-first.local.minimum(score)

    return(list(DoF = DoF, par = par, score=score))
}
