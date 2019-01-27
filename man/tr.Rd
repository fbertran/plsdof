\name{tr}
\Rdversion{1.1}
\alias{tr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Trace of a matrix
}
\description{
This function computes the trace of a matrix.
}
\usage{
tr(M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{square matrix
}
}

\value{
The trace of the matrix M.
}

\author{
Nicole Kraemer
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
M<-matrix(rnorm(8*8),ncol=8)
tr.M<-tr(M)
}

\keyword{math}
