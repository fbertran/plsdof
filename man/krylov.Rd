\name{krylov}
\Rdversion{1.1}
\alias{krylov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Krylov sequence
}
\description{
This function computes the Krylov sequence of a matrix and a vector.
}
\usage{
krylov(A,b,m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{square matrix of dimension p x p.
}
\item{b}{vector of length p}
\item{m}{length of the Krylov sequence}
}

\value{A matrix of size p x m containing the sequence b,Ab,..., A^(m-1)b. 
}

\author{
Nicole Kraemer
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
A<-matrix(rnorm(8*8),ncol=8)
b<-rnorm(8)
K<-krylov(A,b,4)
}

\keyword{math}
