\name{path.plot}
\alias{path.plot}
\title{Coefficients path plot}
\usage{
path.plot(object, Z = NULL, K = NULL, 
          i = NULL, prune = FALSE, cor.max = 0.97, 
          lambda.min = .Machine$double.eps^0.5,
          nbreaks.x=6, ...)
}
\arguments{
\item{object}{An object of the 'LASSO' or 'SSI' class}

\item{Z}{(numeric matrix) Design matrix for random effects. When \code{Z=NULL} an identity matrix is considered (default) thus \ifelse{html}{\out{<b>G</b> = <b>K</b>}}{\eqn{\textbf{G = K}}{G = K}}; otherwise \ifelse{html}{\out{<b>G</b> = <b>Z K Z'</b>}}{\eqn{\textbf{G} = \textbf{Z K Z'}}{G = Z K Z'}} is used. Only needed for a \code{fm} object of the class 'SSI'}

\item{K}{(numeric matrix) Kinship relationships. This can be a name of a binary file where the matrix is stored. Only needed for a \code{fm} object of the class 'SSI'}

\item{i}{(integer vector) Index a response variable (columns of matrix \code{Gamma}) for an object of the class 'LASSO'. Index testing elements (stored in \code{object$tst}) for an object of the class 'SSI'. Default \code{i=NULL} will consider either all columns in matrix \code{Gamma} or all elements in \code{object$tst}, respectively}

\item{prune}{\code{TRUE} or \code{FALSE} to whether prune within groups of correlated coefficients, keeping only one per group. A group of coefficients that are highly correlated are likely to overlap in the plot}

\item{cor.max}{(numeric) Correlation threshold to prune within groups of correlated coefficients}

\item{lambda.min}{(numeric) Minimum value of lambda to show in the plot as \code{-log(lambda)}. This prevents \code{-log(lambda)} going to infinite for near-zero lambda values}

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
  
  index = which(Y$trial \%in\% 1:6)       # Use only a subset of data
  Y = Y[index,]
  X = scale(X_E1[index,])               # Reflectance data
  M = scale(M[index,])/sqrt(ncol(M))    # Subset and scale markers
  G = tcrossprod(M)                     # Genomic relationship matrix
  y = as.vector(scale(Y[,'E1']))        # Subset response variable
  
  # Sparse phenotypic regression
  fm1 = LARS(var(X),cov(X,y))
  
  # Sparse family index
  fm2 = SSI(y,K=G,tst=1:10,trn=11:50)
  \donttest{
  path.plot(fm1)
  path.plot(fm2, prune=TRUE)
  path.plot(fm2, K=G, prune=TRUE, cor.max=0.9)
  }
  # Path plot for the first individual in testing set for the SSI
  path.plot(fm2, K=G, i=1)
}
\author{
Marco Lopez-Cruz (\email{maraloc@gmail.com}) and Gustavo de los Campos
}
