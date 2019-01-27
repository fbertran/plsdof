\name{pls.model}
\Rdversion{1.1}
\alias{pls.model}
\title{Partial Least Squares
}
\description{This function computes the Partial Least Squares fit.
}
\usage{
pls.model(X,y,m,Xtest=NULL,ytest=NULL,compute.DoF,compute.jacobian,use.kernel,method.cor)
}
\arguments{
 \item{X}{matrix of predictor observations.
}
  \item{y}{vector of response observations. The length of \code{y} is the same as the number of rows of \code{X}.
}
  \item{m}{maximal number of Partial Least Squares components. Default is \code{m=min(ncol(X),nrow(X)-1)}.
}
  \item{Xtest}{optional matrix of test observations. Default is \code{Xtest=NULL}.
}
  \item{ytest}{optional vector of test observations. Default is \code{ytest=NULL}.
}
  \item{compute.DoF}{Logical variable. If \code{compute.DoF=TRUE}, the Degrees of Freedom of Partial Least Squares are computed. Default is \code{compute.DoF=FALSE}.
}
\item{compute.jacobian}{Should the first derivative of the regression coefficients be computed as well? Default is \code{FALSE}}
  \item{use.kernel}{Should the kernel representation be used to compute the solution. Default is \code{FALSE}.
}
 \item{method.cor}{How should the correlation to the response be computed? Default is ''pearson''.}
}
\details{This function computes the Partial Least Squares fit and its Degrees of Freedom. Further, it returns the 
regression coefficients and various quantities that are needed for model selection in combination with \code{information.criteria}. 
}
\value{
\item{coefficients}{matrix of regression coefficients}
\item{intercept}{vector of intercepts}
\item{DoF}{vector of Degrees of Freedom}
\item{RSS}{vector of residual sum of error}
\item{sigmahat}{vector of estimated model error}
\item{Yhat}{matrix of  fitted values}
\item{yhat}{vector of squared length of fitted values}
\item{covariance}{if \code{compute.jacobian} is \code{TRUE}, the function returns the array of covariance matrices for the PLS regression coefficients.}
\code{prediction}{if \code{Xtest} is provided, the predicted y-values for \code{Xtest}.}
\code{mse}{if \code{Xtest} and \code{ytest} are provided, the mean squared error on the test data.} 
\code{cor}{if \code{Xtest} and \code{ytest} are provided, the correlation to the response on the test data.} 
}

\references{
Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial Least Squares Regression". Journal of the American Statistical Association 106 (494)
\url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}

Kraemer, N., Sugiyama, M., Braun, M.L. (2009) "Lanczos Approximations for the Speedup of Partial Least Squares Regression", Proceedings of the 12th International Conference on Artificial Intelligence and Stastistics, 272 - 279 

}
\author{Nicole Kraemer, Mikio L. Braun
}




\seealso{
\code{\link{pls.ic}}, \code{\link{pls.cv}}
}
\examples{
n<-50 # number of observations
p<-15 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)

ntest<-200 #
Xtest<-matrix(rnorm(ntest*p),ncol=p) # test data
ytest<-rnorm(ntest) # test data

# compute PLS + degrees of freedom + prediction on Xtest
first.object<-pls.model(X,y,compute.DoF=TRUE,Xtest=Xtest,ytest=NULL)

# compute PLS + test error
second.object=pls.model(X,y,m=10,Xtest=Xtest,ytest=ytest)
}

\keyword{multivariate}
