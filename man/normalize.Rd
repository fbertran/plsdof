\name{normalize}
\Rdversion{1.1}
\alias{normalize}
\title{Normalization of vectors
}
\description{Normalization of vectors.
}
\usage{
normalize(v, w = NULL)
}
\arguments{
  \item{v}{vector
}
  \item{w}{optional vector
}
}
\details{The vector \code{v} is normalized to length 1. If \code{w} is given, it is normalized by the length of \code{v}.
}
\value{
    \item{v}{normalized \code{v}}
    \item{w}{normalized \code{w}}
}

\author{Nicole Kraemer, Mikio L. Braun
}

\examples{
v<-rnorm(5)
w<-rnorm(10)
dummy<-normalize(v,w)
}

\keyword{math}
