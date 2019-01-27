\name{pls.ic}
\Rdversion{1.1}
\alias{pls.ic}
\title{Model selection for Partial Least Squares based on information criteria
}
\description{This function computes the optimal model parameters using one of three different model selection criteria (aic, bic, gmdl) and based on
two different Degrees of Freedom estimates for PLS.
}
\usage{
pls.ic(X, y, m,criterion="bic",naive,use.kernel,compute.jacobian,verbose) 
}
\arguments{
 \item{X}{matrix of predictor observations.
}
  \item{y}{vector of response observations. The length of \code{y} is the same as the number of rows of \code{X}.
}
  \item{m}{maximal number of Partial Least Squares components. Default is \code{m=ncol(X)}.
}
  \item{criterion}{Choice of the model selection criterion. One of the three options aic, bic, gmdl.
}
  \item{naive}{Use the naive estimate for the Degrees of Freedom? Default is \code{FALSE}.
}
\item{use.kernel}{Use kernel representation? Default is \code{use.kernel=FALSE}.
}
  \item{compute.jacobian}{Should the first derivative of the regression coefficients be computed as well? Default is \code{FALSE}}
\item{verbose}{If \code{TRUE}, the function prints a warning if the algorithms produce negative Degrees of Freedom. Default is \code{TRUE}.}

}
\details{There are two options to estimate the Degrees of Freedom of PLS: \code{naive=TRUE} defines the Degrees of Freedom as
the number of components +1, and \code{naive=FALSE} uses the generalized notion of Degrees of Freedom. If \code{compute.jacobian=TRUE},
the function uses the Lanczos decomposition to derive the Degrees of Freedom, otherwise, it uses the Krylov representation. (See 
Kraemer and Sugiyama (2011) for details.) The latter two methods only differ with respect to the estimation of the noise level.

}
\value{
The function returns an object of class "plsdof".
\item{DoF}{Degrees of Freedom}
\item{m.opt}{optimal number of components}
\item{sigmahat}{vector of estimated model errors}
\item{intercept}{intercept}
\item{coefficients}{vector of regression coefficients}
\item{covariance}{if \code{compute.jacobian=TRUE} and \code{use.kernel=FALSE}, the function returns the covariance matrix of the optimal regression coefficients.}
\item{m.crash}{the number of components for which the algorithm returns negative Degrees of Freedom}

}
\references{
Akaikie, H. (1973) "Information Theory and an Extension of the Maximum Likelihood Principle". Second International Symposium on Information Theory, 267 - 281.

Hansen, M., Yu, B. (2001). "Model Selection and Minimum Descripion Length Principle". Journal of
the American Statistical Association, 96, 746 - 774

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial Least Squares Regression". Journal of the American Statistical Association 106 (494)
\url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}


Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 

Schwartz, G. (1979) "Estimating the Dimension of a Model" Annals of Statistics 26(5), 1651 - 1686.


}
\author{Nicole Kraemer, Mikio L. Braun
}

\seealso{
\code{\link{pls.model}}, \code{\link{pls.cv}}
}
\examples{
n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)

# compute linear PLS
pls.object<-pls.ic(X,y,m=ncol(X))
}

\keyword{multivariate}
