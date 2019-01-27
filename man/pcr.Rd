\name{pcr}
\Rdversion{1.1}
\alias{pcr}
\title{Principal Components Regression
}
\description{This function computes the Principal Components Regression (PCR) fit.
}
\usage{
pcr(X,y,scale,m,eps,supervised)
}
\arguments{
 \item{X}{matrix of predictor observations.
}
  \item{y}{vector of response observations. The length of \code{y} is the same as the number of rows of \code{X}.
}
\item{scale}{Should the predictor variables be scaled to unit variance? Default is \code{TRUE}.}
\item{m}{maximal number of principal components. Default is \code{m=min(ncol(X),nrow(X)-1)}.
}
  \item{eps}{precision. Eigenvalues of the correlation matrix of \code{X} that are smaller than \code{eps} are set to 0. The default value is \code{eps=10^{-6}.}
}
 \item{supervised}{Should the principal components be sorted by decreasing squared correlation to the response? Default is FALSE.} 
}
\details{The function first scales all predictor variables to unit variance, and then computes the PCR fit for all components. Is \code{supervised=TRUE}, we sort the principal correlation according to the squared correlation to the response.
}
\value{
\item{coefficients}{matrix of regression coefficients, including the coefficients of the null model, i.e. the constant model \code{mean(y)}. }
\item{intercept}{vector of intercepts, including the intercept of the null model, i.e. the constant model \code{mean(y)}. }
}


\author{Nicole Kraemer
}




\seealso{
\code{\link{pcr.cv}}, \code{\link{pls.cv}}
}
\examples{
n<-50 # number of observations
p<-15 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)

my.pcr<-pcr(X,y,m=10)

}

\keyword{multivariate}
