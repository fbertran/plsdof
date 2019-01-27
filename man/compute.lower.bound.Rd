\name{compute.lower.bound}
\Rdversion{1.1}
\alias{compute.lower.bound}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Lower bound for the Degrees of Freedom
}
\description{
This function computes the lower bound for the the Degrees of Freedom of PLS with 1 component. }
\usage{
compute.lower.bound(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{matrix of predictor observations.
}
}
\details{
If the decay of the eigenvalues of \code{cor(X)} is not too fast, we can lower-bound the Degrees of Freedom of PLS with 1 component.
Note that we implicitly assume that we use scaled predictor variables to compute the PLS solution.

}

\value{
\item{bound}{logical. bound is \code{TRUE} if the decay of the eigenvalues is slow enough}
\item{lower.bound}{if bound is TRUE, this is the lower bound, otherwise, it is set to -1}
}

\references{

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial Least Squares Regression". Journal of the American Statistical Association 106 (494)
\url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}
}

\author{Nicole Kraemer
}
\seealso{
 \code{\link{pls.model}}
}
\examples{
# Boston Housing data
library(MASS)
data(Boston)
X<-Boston[,-14]
my.lower<-compute.lower.bound(X)
}

\keyword{math}
