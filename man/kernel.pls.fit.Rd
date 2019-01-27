\name{kernel.pls.fit}
\Rdversion{1.1}
\alias{kernel.pls.fit}
\title{Kernel Partial Least Squares Fit
}
\description{
This function computes the Partial Least Squares fit. This algorithm scales mainly in the number of observations.
}
\usage{
kernel.pls.fit(X, y, m, compute.jacobian,DoF.max)
}
\arguments{
 \item{X}{matrix of predictor observations.
}
  \item{y}{vector of response observations. The length of \code{y} is the same as the number of rows of \code{X}.
}
  \item{m}{maximal number of Partial Least Squares components. Default is \code{m}=ncol(X).
}
  \item{compute.jacobian}{Should the first derivative of the regression coefficients be computed as well? Default is \code{FALSE}}.
\item{DoF.max}{upper bound on the Degrees of Freedom. Default is \code{min(ncol(X)+1,nrow(X)-1)}.}
}
\details{We first standardize \code{X} to zero mean and unit variance. 
}
\value{
\item{coefficients}{matrix of regression coefficients}
\item{intercept}{vector of regression intercepts}
\item{DoF}{Degrees of Freedom}
\item{sigmahat}{vector of estimated model error}
\item{Yhat}{matrix of fitted values}
\item{yhat}{vector of squared length of fitted values}
\item{RSS}{vector of residual sum of error}
\item{covariance}{\code{NULL} object.}
\item{TT}{matrix of normalized PLS components}
}

\references{
Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial Least Squares Regression". Journal of the American Statistical Association 106 (494)
\url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}


Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 

}
\author{Nicole Kraemer, Mikio L. Braun
}

\seealso{
\code{\link{linear.pls.fit}}, \code{\link{pls.cv}},\code{\link{pls.model}}, \code{\link{pls.ic}}
}
\examples{
n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)


pls.object<-kernel.pls.fit(X,y,m=5,compute.jacobian=TRUE)


}
\keyword{multivariate}
