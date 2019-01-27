\name{coef.plsdof}
\Rdversion{1.1}
\alias{coef.plsdof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Regression coefficients
}
\description{
This function returns the regression coefficients of a plsdof-object.
}
\usage{
\method{coef}{plsdof}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class "plsdof" that is returned by the functions \code{pls.ic} and \code{pls.cv}.
}
\item{...}{additional parameters}
}
\details{
The function returns the regression coefficients (without intercept) for the optimal number of components.
}
\value{
regression coefficients.
}
\references{

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial Least Squares Regression". Journal of the American Statistical Association 106 (494)
\url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 


}
\author{
Nicole Kraemer
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{vcov.plsdof}}, \code{\link{pls.model}}, \code{\link{pls.ic}}, \code{\link{pls.cv}}
}
\examples{

n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)


pls.object<-pls.ic(X,y,criterion="bic")
mycoef<-coef(pls.object)
}

\keyword{models}
