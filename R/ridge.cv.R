#' Ridge Regression.
#' 
#' This function computes the optimal ridge regression model based on
#' cross-validation.
#' 
#' Based on the regression coefficients \code{coefficients.jackknife} computed
#' on the cross-validation splits, we can estimate their mean and their
#' variance using the jackknife. We remark that under a fixed design and the
#' assumption of normally distributed \code{y}-values, we can also derive the
#' true distribution of the regression coefficients.
#' 
#' @param X matrix of input observations. The rows of \code{X} contain the
#' samples, the columns of \code{X} contain the observed variables
#' @param y vector of responses. The length of y must equal the number of rows
#' of X
#' @param lambda Vector of penalty terms.
#' @param scale Scale the columns of X? Default is scale=TRUE.
#' @param k Number of splits in \code{k}-fold cross-validation. Default value
#' is \code{k}=10.
#' @param plot.it Plot the cross-validation error as a function of
#' \code{lambda}? Default is FALSE.
#' @param groups an optional vector with the same length as \code{y}. It
#' encodes a partitioning of the data into distinct subgroups. If \code{groups}
#' is provided, \code{k=10} is ignored and instead, cross-validation is
#' performed based on the partioning. Default is \code{NULL}.
#' @param method.cor How should the correlation to the response be computed?
#' Default is ''pearson''.
#' @param compute.jackknife Logical. If \code{TRUE}, the regression
#' coefficients on each of the cross-validation splits is stored. Default is
#' \code{TRUE}.
#' @return \item{cv.error.matrix}{matrix of cross-validated errors based on
#' mean squared error. A row corresponds to one cross-validation split.}
#' \item{cv.error}{vector of cross-validated errors based on mean squared
#' error} \item{lambda.opt}{optimal value of \code{lambda}, based on mean
#' squared error} \item{intercept}{intercept of the optimal model, based on
#' mean squared error} \item{coefficients}{vector of regression coefficients of
#' the optimal model, based on mean squared error}
#' \item{cor.error.matrix}{matrix of cross-validated errors based on
#' correlation. A row corresponds to one cross-validation split.}
#' \item{cor.error}{vector of cross-validated errors based on correlation}
#' \item{lambda.opt.cor}{optimal value of \code{lambda}, based on correlation}
#' \item{intercept.cor}{intercept of the optimal model, based on correlation}
#' \item{coefficients.cor}{vector of regression coefficients of the optimal
#' model, based on mean squared error} \item{coefficients.jackknife}{Array of
#' the regression coefficients on each of the cross-validation splits. The
#' dimension is \code{ncol(X) x length(lambda) x k}.}
#' @author Nicole Kraemer
#' @seealso \code{\link{pls.cv}}, \code{\link{pcr.cv}},
#' \code{\link{benchmark.regression}}
#' @keywords multivariate
#' @examples
#' 
#' n<-100 # number of observations
#' p<-60 # number of variables
#' X<-matrix(rnorm(n*p),ncol=p) 
#' y<-rnorm(n)
#' ridge.object<-ridge.cv(X,y)
#' 
#' @export ridge.cv
ridge.cv<-
function(X,y,lambda=NULL,scale=TRUE,k=10,plot.it=FALSE,groups=NULL,method.cor="pearson",compute.jackknife=TRUE){
    if (is.null(lambda)==TRUE){
        ss<-seq(-10,-1,length=1000)
        ss<-10^ss
        n<-nrow(X)
        nn<-n- floor(n/k)
        lambda<-ss*nn*ncol(X)
    }
	p<-ncol(X)
	coefficients.jackknife=NULL
       if (compute.jackknife==TRUE){
   coefficients.jackknife<-array(dim=c(p,length(lambda),k))
	}
    
    n<-nrow(X)
if (is.null(groups) == FALSE) {
        f = as.factor(groups)
        k = length(levels(f))
        my.names = levels(f)
        all.folds <- split(1:n, f)
    }
    if (is.null(groups) == TRUE) {
        f <- rep(1:k, length = n)
        my.names <- 1:k
       all.folds <- split(sample(1:n), f)
    }
    cv.error.matrix<-cor.error.matrix<-matrix(0,k,length(lambda))
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain=X[-omit,,drop=FALSE]
        ytrain=y[-omit]
        Xtest=X[omit,,drop=FALSE]
        ytest=y[omit]
        ll<-lm.ridge(ytrain~Xtrain,scale=scale,lambda=lambda)
        coef.ll<-coef(ll)
	if (compute.jackknife==TRUE){
        	coefficients.jackknife[,,i]<-t(coef.ll[,-1])
	}        
	res<-matrix(,length(ytest),length(lambda))
        pred<-t(matrix(coef.ll[,1],nrow=length(lambda),ncol=length(ytest))) + Xtest%*%t(coef.ll[,-1])
	res<-pred-matrix(ytest,nrow=length(ytest),ncol=length(lambda))
	cv.error.matrix[i,]= apply(res^2,2,mean)  
        for (j in 1:length(lambda)){   
	cor.error.matrix[i,j]<-cor(pred[,j],ytest)
        }
    }
    cv.error<-apply(cv.error.matrix,2,mean)
    cor.error<-apply(cor.error.matrix,2,mean)
    lambda.opt<-lambda[which.min(cv.error)]
    lambda.opt.cor<-lambda[which.max(cor.error)]
    if (plot.it==TRUE){
        plot(lambda,cv.error,type="l",ylim="mean squared error")
    }
    rr<-lm.ridge(y~X,scale=scale,lambda=lambda.opt)
    coefficients<-coef(rr)
    intercept<-coefficients[1]
    coefficients<-coefficients[-1]
    rr.cor<-lm.ridge(y~X,scale=scale,lambda=lambda.opt.cor)
    coefficients.cor<-coef(rr.cor)
    intercept.cor<-coefficients.cor[1]
    coefficients.cor<-coefficients.cor[-1]
    return(list(cv.error=cv.error,cor.error=cor.error,cv.error.matrix=cv.error.matrix,cor.error.matrix=cor.error.matrix,intercept=intercept,coefficients=coefficients,lambda.opt=lambda.opt,intercept.cor=intercept.cor,coefficients.cor=coefficients.cor,lambda.opt.cor=lambda.opt.cor,coefficients.jackknife=coefficients.jackknife))
}
