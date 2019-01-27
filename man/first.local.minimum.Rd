\name{first.local.minimum}
\Rdversion{1.1}
\alias{first.local.minimum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Index of the first local minimum.
}
\description{
This function computes the index of the first local minimum.
}
\usage{
first.local.minimum(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector.

}
  
}

\value{
the index of the first local minimum of \code{x}.
}


\references{
Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial Least Squares Regression". Journal of the American Statistical Association. ahead of print 106 (494)
\url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}

}
\author{
Nicole Kraemer
}


\examples{
v<-rnorm(30)
out<-first.local.minimum(v)
}
\keyword{math}
