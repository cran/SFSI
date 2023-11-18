
\name{9. SSI accuracy vs penalization plot}
\alias{plot.SSI}
\title{Accuracy vs penalization plot}
\usage{
\method{plot}{SSI}(..., x.stat = c("nsup","lambda"),
          y.stat = c("accuracy","MSE"),
          nbreaks.x=7)
}
\arguments{
%\item{x}{Object of the class 'SSI'. Additional objects can be passed through arguments \code{y} and '\code{...}'}

%\item{y}{An (optional) object of the class 'SSI'}

\item{...}{Other arguments to be passed: 
\itemize{
 \item One or more objects of the class 'SSI'
 \item Other arguments for method \code{plot}: 'xlab', 'ylab', 'main', 'lwd', 'xlim', 'ylim'
}
}

\item{x.stat}{(character) Either 'nsup' (number of non-zero regression coefficients entering in the prediction of a given testing individual) or 'lambda' (penalization parameter in log scale) to plot in the x-axis}

\item{y.stat}{(character) Either 'accuracy' (correlation between observed and predicted values) or 'MSE' (mean squared error) to plot in the y-axis}

\item{nbreaks.x}{(integer) Number of breaks in the x-axis}

}
\value{
Creates a plot of either accuracy or MSE versus either the support set size (average number of predictors with non-zero regression coefficient) or versus lambda. 
}
\description{
Accuracy as a function of the penalization plot for an object of the class 'SSI'
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
  
  plot(fm)
  plot(fm, main=expression('corr('*y[obs]*','*y[pred]*') vs sparsity'))   
  plot(fm, y.stat="MSE",ylab='Mean Square Error', xlab='Sparsity')     
  
}
