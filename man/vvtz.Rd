\name{vvtz}
\Rdversion{1.1}
\alias{vvtz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Projectin operator
}
\description{
This function computes the projection operator
\deqn{P_V z= V V^\top z}
}
\usage{
vvtz(v, z)
}
\arguments{
 \item{v}{orthonormal basis of the space on which \code{z} is projected. \code{v} is either a matrix or a vector.
}
  \item{z}{vector that is projected onto the columns of \code{v}
}
}
\details{The above formula is only valid if the columns of \code{v} are normalized and mutually orthogonal.
}
\value{value of the projection operator

}

\author{Nicole Kraemer
}

\seealso{
\code{\link{dvvtz}}
}
\examples{
# generate random orthogonal vectors
X<-matrix(rnorm(10*100),ncol=10) 	# random data
S<-cor(X) 				# correlation matrix of data
v<-eigen(S)$vectors[,1:3]		# first three eigenvectors of correlation matrix
z<-rnorm(10)				# random vector z
projection.z<-vvtz(v,z)

}

\keyword{math}
