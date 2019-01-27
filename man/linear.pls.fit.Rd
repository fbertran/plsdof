\name{linear.pls}
\Rdversion{1.1}
\alias{linear.pls.fit}
\title{Linear Partial Least Squares Fit
}
\description{This function computes the  Partial Least Squares solution and the first derivative of the regression coefficients. This implementation scales mostly in the number of variables
}
\usage{
linear.pls.fit(X, y, m,compute.jacobian,DoF.max)
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
\item{Yhat}{matrix of  fitted values}
\item{yhat}{vector of squared length of fitted values}
\item{RSS}{vector of residual sum of error}
\code{covariance}{if \code{compute.jacobian} is \code{TRUE}, the function returns the array of covariance matrices for the PLS regression coefficients.}
\item{TT}{matrix of normalized PLS components}
}

\references{
Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial Least Squares Regression". Journal of the American Statistical Association 106 (494)
\url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}

}

\author{Nicole Kraemer
}




\seealso{
\code{\link{kernel.pls.fit}}, \code{\link{pls.cv}},\code{\link{pls.model}}, \code{\link{pls.ic}}
}
\examples{
n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)


pls.object<-linear.pls.fit(X,y,m=5,compute.jacobian=TRUE)


}

\keyword{multivariate}
