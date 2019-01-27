\name{pls.cv}
\Rdversion{1.1}
\alias{pls.cv}
\title{Model selection for Partial Least Squares based on cross-validation
}
\description{This function computes the optimal model parameter using cross-validation.
}
\usage{
pls.cv(X, y, k, groups,m,use.kernel=FALSE,compute.covariance=FALSE,method.cor)
}
\arguments{
 \item{X}{matrix of predictor observations.
}
  \item{y}{vector of response observations. The length of \code{y} is the same as the number of rows of \code{X}.
}
\item{k}{number of cross-validation splits. Default is 10.}
\item{groups}{an optional vector with the same length as \code{y}. It encodes a partitioning of the data into distinct subgroups. If \code{groups} is provided, \code{k=10} is ignored and instead, cross-validation is performed based on the partioning. Default is \code{NULL}.}
  \item{m}{maximal number of Partial Least Squares components. Default is \code{m=ncol(X)}.
}
  \item{use.kernel}{Use kernel representation? Default is \code{use.kernel=FALSE}.
}
  
\item{compute.covariance}{If \code{TRUE}, the function computes the covariance for the cv-optimal regression coefficients.} 
 \item{method.cor}{How should the correlation to the response be computed? Default is ''pearson''.}
}
\details{The data are centered and scaled to unit variance prior to the PLS algorithm. It is possible to estimate the covariance
matrix of the cv-optimal regression coefficients (\code{compute.covariance=TRUE}). Currently, this is only implemented if
\code{use.kernel=FALSE}. 
}
\value{
\item{cv.error.matrix}{matrix of cross-validated errors based on mean squared error. A row corresponds to one cross-validation split.}
\item{cv.error}{vector of cross-validated errors based on mean squared error}
\item{m.opt}{optimal number of components based on mean squared error}
\item{intercept}{intercept of the optimal model, based on mean squared error}
\item{coefficients}{vector of regression coefficients of the optimal model, based on mean squared error}
\item{cor.error.matrix}{matrix of cross-validated errors based on  correlation. A row corresponds to one cross-validation split.}
\item{cor.error}{vector of cross-validated errors based on correlation}
\item{m.opt.cor}{optimal number of components based on  correlation}
\item{intercept.cor}{intercept of the optimal model, based on  correlation}
\item{coefficients.cor}{vector of regression coefficients of the optimal model, based on mean squared error}
\item{covariance}{If \code{TRUE} and \code{use.kernel=FALSE}, the covariance of the cv-optimal regression coefficients (based on mean squared error) is returned.}


}
\references{

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial Least Squares Regression". Journal of the American Statistical Association 106 (494)
\url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 


}
\author{Nicole Kraemer, Mikio L. Braun
}

\seealso{
\code{\link{pls.model}}, \code{\link{pls.ic}}
}
\examples{
n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)

# compute linear PLS
pls.object<-pls.cv(X,y,m=ncol(X))

# define random partioning
groups<-sample(c("a","b","c"),n,replace=TRUE)
pls.object1<-pls.cv(X,y,groups=groups)
}
\keyword{multivariate}
