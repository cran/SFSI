
\name{Multitrait SSI accuracy vs penalization plot}
\alias{multitrait.plot}
\title{Accuracy vs penalization from Multitrait SSI}
\usage{
multitrait.plot(object, x.stat = c("nsup","lambda"),
                y.stat = c("accuracy","MSE"),
                line.color = "orange",
                nbreaks.x = 7, ...)

}
\arguments{
\item{object}{An object of the class 'SSI' for a multitrait case}

\item{x.stat}{(character) Either 'nsup' (number of non-zero regression coefficients entering in the prediction of a given testing individual) or 'lambda' (penalization parameter in log scale) to plot in the x-axis}

\item{y.stat}{(character) Either 'accuracy' (correlation between observed and predicted values) or 'MSE' (mean squared error) to plot in the y-axis}

\item{line.color}{(character) Color of the lines}

\item{nbreaks.x}{(integer) Number of breaks in the x-axis}

\item{...}{Other arguments for method \code{plot}: 'xlab', 'ylab', 'main', 'lwd', 'xlim', 'ylim'}
}
\value{
Creates a plot of either accuracy or MSE versus either the support set size (average number of predictors with non-zero regression coefficient) or versus lambda. This is done separately for each trait 
}
\description{
Visualizing results from an object of the class 'SSI'
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  index = which(Y$trial \%in\% 1:9)   # Use only a subset of data
  Y = Y[index,]
  M = scale(M[index,])/sqrt(ncol(M))  # Subset and scale markers
  G = tcrossprod(M)                   # Genomic relationship matrix
  y = scale(Y[,4:5])    # Response variable
  
  # Sets (testing=0, training=1)
  trn_tst = matrix(NA,ncol=ncol(y),nrow=nrow(y))
  trn_tst[,1] = ifelse(Y$trial \%in\% 2, 0, 1)
  trn_tst[,2] = ifelse(Y$trial \%in\% 2, 0, 1)
  \donttest{
  fm = SSI(y, K=G, trn_tst=trn_tst, mc.cores=1)
  multitrait.plot(fm)
  }
}
