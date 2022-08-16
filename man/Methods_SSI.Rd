\name{Methods_SSI}
\alias{summary.SSI}
\alias{plot.SSI}
\alias{coef.SSI}
\alias{fitted.SSI}
\title{SSI methods}
\usage{
\method{coef}{SSI}(object, ..., df=NULL, i=NULL)

\method{fitted}{SSI}(object, ...)

\method{summary}{SSI}(object, ...)
  
\method{plot}{SSI}(..., py=c("accuracy","MSE"), nbreaks.x=6)
}
\arguments{
\item{object}{An object of the class 'SSI'}

\item{...}{Arguments to be passed: 
\itemize{
 \item One or more objects of the class 'SSI' (for method \code{plot})
 \item Other arguments for method \code{plot}: 'xlab', 'ylab', 'main', 'lwd', 'xlim', 'ylim'
 \item An optional vector of observations \code{y} (with a similar heritability as the one declared in 'SSI' function) for methods \code{summary}, \code{plot}, and \code{fitted}
}
}

\item{df}{(numeric) Average (across testing individuals) number of non-zero regression coefficients}

\item{i}{(integer vector) Index testing elements (stored in \code{object$tst}) to be considered. Default \code{i=NULL} will consider all elements in \code{object$tst}}

\item{py}{(character) Either 'accuracy' (correlation between observed and predicted values) or 'MSE' (mean squared error) to plot in the y-axis}

\item{nbreaks.x}{(integer) Number of breaks in the x-axis}

}
\value{
Method \code{fitted} returns a matrix with the predicted values for each individual in the testing set (in rows) for each value of lambda (in columns).

Method \code{coef} (list of matrices) returns the regression coefficients for each testing set individual (elements of the list). Each matrix contains the coefficients for each value of lambda (in rows) associated to each training set individual (in columns). If \code{tst} is specified, the elements of the list will correspond only to the testing individuals given in \code{tst}. If \code{df} is specified, only the coefficients for the lambda associated to \code{df} are returned as a 'matrix' with testing individuals in rows. 

Method \code{summary} returns a list object containing:
\itemize{
  \item \code{lambda}: (vector) sequence of values of lambda used in the coefficients' estimation.
  \item \code{df}: (vector) degrees of freedom (across testing individuals) at each solution associated to each value of lambda.
  \item \code{accuracy}: (vector) correlation between observed and predicted values associated to each value of lambda.
  \item \code{MSE}: (vector) mean squared error associated to each value of lambda.
  \item \code{optCOR}: (vector) summary of the SSI with maximum accuracy.
  \item \code{optMSE}: (vector) summary of the SSI with minimum MSE.
}
 
Method \code{plot} creates a plot of either accuracy or MSE versus the (average across testing individuals) number of predictors (with non-zero regression coefficient) and versus lambda. 
}
\description{
Useful methods for retrieving, summarizing and visualizing important results from an object of the class 'SSI'
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  index = which(Y$trial \%in\% 1:6)   # Use only a subset of data
  Y = Y[index,]
  M = scale(M[index,])/sqrt(ncol(M))  # Subset and scale markers
  G = tcrossprod(M)                   # Genomic relationship matrix
  y = as.vector(scale(Y[,"E1"]))      # Scale response variable
  
  # Training and testing sets
  tst = which(Y$trial == 2)
  trn = which(Y$trial != 2)
  
  fm1 = SSI(y,K=G,theta=1,b=0,tst=tst,trn=trn)
  
  uHat = fitted(fm1)                # Predicted values for each testing element
  out = summary(fm1)                # Useful function to get results
  corTST = out$accuracy             # Testing set accuracy (correlation cor(y,yHat))
  out$optCOR                        # SSI with maximum accuracy
  out$optMSE                        # SSI with minimum MSE
  B = coef(fm1)                     # Regression coefficients for all tst
  B = coef(fm1, i=1)                # Regression coefficients for first tst (tst[1])
  B = coef(fm1, df=10)              # Regression coefficients for which df=10
  plot(fm1,main=expression('corr('*y[obs]*','*y[pred]*') vs sparsity'))   
  plot(fm1,py="MSE",ylab='Mean Square Error', xlab='Sparsity')     
  
}
\author{
Marco Lopez-Cruz (\email{maraloc@gmail.com}) and Gustavo de los Campos
}
