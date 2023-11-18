
\name{Summarize SSI results}
\alias{summary.SSI}
\alias{coef.SSI}
\alias{fitted.SSI}
\title{SSI methods}
\usage{
\method{coef}{SSI}(object, ...)

\method{fitted}{SSI}(object, ...)

\method{summary}{SSI}(object, ...)

}
\arguments{
\item{object}{An object of the class 'SSI'}

\item{...}{Other arguments to be passed to \code{coef} method: 
\itemize{
 \item \code{nsup}: (numeric) Average (across testing individuals) number of non-zero regression coefficients. Only the coefficients for the lambda associated to \code{nsup} are returned as a 'matrix' with testing individuals in rows
 \item \code{iy}: (integer vector) Index testing elements (stored in \code{object$tst}) to be considered. Only coefficients corresponding to the testing individuals \code{object$tst[i]} are returned
}
For \code{fitted} and \code{summary} methods:
\itemize{
 \item \code{y}: (numeric vector) An optional response matrix
}
}
}
\value{
Method \code{fitted} returns a matrix with the predicted values for each individual in the testing set (in rows) for each value of lambda (in columns).

Method \code{coef} (list of matrices) returns the regression coefficients for each testing set individual (elements of the list). Each matrix contains the coefficients for each value of lambda (in rows) associated to each training set individual (in columns). 

Method \code{summary} returns a list object containing:
\itemize{
  \item \code{lambda}: (vector) sequence of values of lambda used in the coefficients' estimation.
  \item \code{nsup}: (vector) Number of non-zero coefficients (across testing individuals) at each solution associated to each value of lambda.
  \item \code{accuracy}: (vector) correlation between observed and predicted values associated to each value of lambda.
  \item \code{MSE}: (vector) mean squared error associated to each value of lambda.
  \item \code{optCOR}: (vector) summary of the SSI with maximum accuracy.
  \item \code{optMSE}: (vector) summary of the SSI with minimum MSE.
}
 
}
\description{
Useful methods for retrieving and summarizing important results from the 'SSI' function's output
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  index = which(Y$trial \%in\% 1:6)   # Use only a subset of data
  Y = Y[index,]
  M = scale(M[index,])/sqrt(ncol(M))  # Subset and scale markers
  G = tcrossprod(M)                   # Genomic relationship matrix
  y = as.vector(scale(Y[,"E1"]))      # Scale response variable
  
  # Sets (testing=0, training=1)
  trn_tst = ifelse(Y$trial == 2, 0, 1)
  
  fm = SSI(y,K=G,varU=0.4,varE=0.6,b=0,trn_tst=trn_tst)
  
  uHat = fitted(fm)             # Predicted values for each testing element
  out = summary(fm)             # Useful function to get results
  corTST = out$accuracy         # Testing set accuracy (correlation cor(y,yHat))
  out$optCOR                    # SSI with maximum accuracy
  out$optMSE                    # SSI with minimum MSE
  B = coef(fm)                  # Regression coefficients for all tst
  B = coef(fm, iy=1)            # Coefficients for first tst (tst[1])
  B = coef(fm, ilambda=10)      # Coefficients associated to the 10th lambda
  B = coef(fm, nsup=10)         # Coefficients for which nsup=10
  
}

