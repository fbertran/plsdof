\name{information.criteria}
\Rdversion{1.1}
\alias{information.criteria}
\title{Information criteria
}
\description{
This function computes the optimal model parameters using three different model selection criteria (aic, bic, gmdl).
}
\usage{
information.criteria(RSS, DoF, yhat, sigmahat, n,criterion="bic")
}
\arguments{
    \item{RSS}{vector of residual sum of squares.
}
  
  \item{DoF}{vector of Degrees of Freedom. The length of \code{DoF} is the same as the length of \code{RSS}.
}
\item{yhat}{vector of squared norm of yhat. The length of \code{yhat} is the same as the length of \code{RSS}. It is only needed for gmdl. Default value is \code{NULL}. 
}
\item{sigmahat}{Estimated model error. The length of \code{sigmahat} is the same as the length of \code{RSS}.
}
\item{n}{number of observations.
}
\item{criterion}{one of the options "aic", "bic" and "gmdl".}
}
\details{
The Akaike information criterion (aic) is defined as 
\deqn{{aic}= \frac{{RSS}}{n} + 2\frac{{DoF}}{n} \sigma^ 2\,.}
The Bayesian information criterion (bic)  is defined as 
\deqn{{bic}= \frac{{RSS}}{n} + log(n)\frac{{DoF}}{n} \sigma^ 2\,.}
The generalized minimum description length (gmdl)  is defined as
\deqn{gmdl=\frac{n}{2}log(S)+\frac{DoF}{2}log(F)+\frac{1}{2}log(n)}
with \deqn{S=\hat \sigma ^2}
Note that it is also possible to use the function \code{information.criteria} for other regression methods than Partial Least Squares.

}

\value{
 \item{DoF}{degrees of freedom}
\item{score}{vector of the model selection criterion}
\item{par}{index of the first local minimum of \code{score}}

}
\references{
Akaikie, H. (1973) "Information Theory and an Extension of the Maximum Likelihood Principle". Second International Symposium on Information Theory, 267 - 281.

Hansen, M., Yu, B. (2001). "Model Selection and Minimum Descripion Length Principle". Journal of
the American Statistical Association, 96, 746 - 774

Kraemer, N., Sugiyama M. (2011). "The Degrees of Freedom of Partial Least Squares Regression". Journal of the American Statistical Association 106 (494)
\url{https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10107}

Kraemer, N., Braun, M.L. (2007) "Kernelizing PLS, Degrees of Freedom, and Efficient Model Selection", Proceedings of the 24th International Conference on Machine Learning, Omni Press, 441 - 448 

Schwartz, G. (1979) "Estimating the Dimension of a Model" Annals of Statistics 26(5), 1651 - 1686.

}
\author{
Nicole Kraemer, Mikio Braun
}

\seealso{
\code{\link{pls.ic}}
}
\examples{
## This is an internal function called by pls.ic

}

\keyword{model}
