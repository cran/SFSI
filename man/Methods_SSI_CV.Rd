\name{Methods_SSI_CV}
\alias{summary.SSI_CV}
\alias{plot.SSI_CV}
\title{SSI_CV methods}
\usage{
\method{summary}{SSI_CV}(object, ...)

\method{plot}{SSI_CV}(..., py=c("accuracy","MSE"), title=NULL, showFolds=FALSE)
}
\arguments{
\item{object}{An object of the class 'SSI_CV'}

\item{...}{Arguments to be passed: 
\itemize{
 \item Not needed for method \code{summary}
 \item \code{object}: One or more objects of the class 'SSI_CV' (for method \code{plot})
}
}

\item{py}{(character) Either 'accuracy' or 'MSE' to plot the correlation between observed and predicted values or the mean squared error, respectively, in the y-axis}

\item{title}{(character/expression) Title of the plot}

\item{showFolds}{\code{TRUE} or \code{FALSE} to whether add results for individuals folds}
}
\value{
Method \code{summary} returns a list object containing:
\itemize{
  \item \code{lambda}: (matrix) sequence of (average across folds) values of lambda (in columns) used in each CV partition (in rows).
  \item \code{df}: (matrix) degrees of freedom (average across folds) at each solution associated to each value of lambda for each CV partition (in rows).
  \item \code{accuracy}: (matrix) correlation between observed and predicted values (average across folds) given by each value of lambda (in columns) in each CV partition (in rows).
  \item \code{MSE}: (matrix) mean squared error (average across folds) given by each value of lambda (in columns) in each CV partition (in rows).
  \item \code{optCOR}: (matrix) summary of the SSI with maximum accuracy within and across CV partitions (in rows).
  \item \code{optMSE}: (matrix) summary of the SSI with minimum MSE within and across CV partitions (in rows).
}
 
Method \code{plot} creates a plot of either accuracy or MSE versus the (average across folds and partitions) number of predictors (with non-zero regression coefficient) and versus lambda. 
}
\description{
Useful methods for retrieving, summarizing and visualizing important results from an object of the class 'SSI_CV'
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  index = which(Y$CV \%in\% 1:2)        # Use only a subset of data
  M = scale(M[index,])/sqrt(ncol(M))    # Subset and scale markers
  G = tcrossprod(M)                     # Genomic relationship matrix
  y = as.vector(scale(Y[index,"E1"]))   # Subset ans scale response variable
 
  trn = seq(1,length(y),by=2)
  fm1 = SSI_CV(y,K=G,trn=trn,nFolds=5,nCV=1)
  
  out = summary(fm1)    # Useful results
  out$accuracy          # Testing set accuracy (cor(y,yHat))
  out$optCOR            # SSI with maximum accuracy
  out$optMSE            # SSI with minimum MSE
  
  plot(fm1,title=expression('corr('*y[obs]*','*y[pred]*') vs sparsity'))   
  plot(fm1,showFolds=TRUE)     
}
\author{
Marco Lopez-Cruz (\email{maraloc@gmail.com}) and Gustavo de los Campos
}
