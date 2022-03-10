\name{path.plot}
\alias{path.plot}
\title{Coefficients path plot}
\usage{
path.plot(object, Z = NULL, K = NULL, tst = NULL, 
          cor.max = 0.85, nbreaks.x=6, ...)
}
\arguments{
\item{object}{An object of the 'LASSO' or 'SSI' class}

\item{Z}{(numeric matrix) Design matrix for random effects. When \code{Z=NULL} an identity matrix is considered (default) thus \ifelse{html}{\out{<b>G</b> = <b>K</b>}}{\eqn{\textbf{G = K}}{G = K}}; otherwise \ifelse{html}{\out{<b>G</b> = <b>Z K Z'</b>}}{\eqn{\textbf{G} = \textbf{Z K Z'}}{G = Z K Z'}} is used. Only needed for a \code{fm} object of the class 'SSI'}

\item{K}{(numeric matrix) Kinship relationships. This can be a name of a binary file where the matrix is stored. Only needed for a \code{fm} object of the class 'SSI'}

\item{tst}{(integer vector) Which elements from vector \code{y} (stored in \code{fm$y}) are in testing set and to plot. They must be contained in \code{fm$tst}. Default \code{tst=NULL} will consider the whole vector \code{fm$tst} to plot}

\item{cor.max}{(numeric) Maximum correlation allowed for two different coefficients. A group of coeffcients with a correlation greater than \code{cor.max} are likely to overlap in the plot thus only one is kept}

\item{nbreaks.x}{(integer) Number of breaks in the x-axis}

\item{...}{Other arguments for method \code{plot}: 'xlab', 'ylab', 'main', 'lwd'

}
}
\value{
Returns the plot of the coefficients' evolution path along the regularization parameter
}
\description{
Coefficients evolution path plot from an object of the class 'LASSO' or 'SSI'
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  index = which(Y$CV == 1)              # Use only a subset of data
  M = scale(M[index,])/sqrt(ncol(M))    # Subset and scale markers
  G = tcrossprod(M)                     # Genomic relationship matrix
  y = as.vector(scale(Y[index,"E1"]))   # Subset response variable
  X = scale(X_E1[index,])               # Reflectance data
  
  # Sparse phenotypic regression
  fm1 = LARS(var(X),cov(y,X))
  
  # Sparse family index
  fm2 = SSI(y,K=G,tst=1:10,trn=11:50)
  
  \donttest{
  path.plot(fm1)
  path.plot(fm2, cor.max=0.6)
  path.plot(fm2, K=G, cor.max=0.6)
  }
  
  # Path plot for the first individual in testing set for the SSI
  path.plot(fm2, K=G, tst=fm2$tst[1])
}
\author{
Marco Lopez-Cruz (\email{maraloc@gmail.com}) and Gustavo de los Campos
}
