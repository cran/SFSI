\name{net.plot}
\alias{net.plot}
\title{Network plot}
\usage{
net.plot(object, Z = NULL, K = NULL, tst = NULL,
        show.names = FALSE, group = NULL, group.shape = NULL,
        set.color = NULL, set.size = NULL, df = NULL, main,
        axis.labels = TRUE, curve = FALSE, bg.color = "white",
        unified = TRUE, ntst = 36, line.color = "gray80", 
        line.tick = 0.3, legend.pos="right", point.color = "gray20",
        sets = c("Testing","Supporting","Non-active"), 
        eps = .Machine$double.eps)
}
\arguments{
\item{object}{An object of the 'SSI' class or a matrix of coefficients}

\item{Z}{(numeric matrix) Design matrix for random effects. When \code{Z=NULL} an identity matrix is considered (default) thus \ifelse{html}{\out{<b>G</b> = <b>K</b>}}{\eqn{\textbf{G = K}}{G = K}}; otherwise \ifelse{html}{\out{<b>G</b> = <b>Z K Z'</b>}}{\eqn{\textbf{G} = \textbf{Z K Z'}}{G = Z K Z'}} is used}

\item{K}{(numeric matrix) Kinship relationships. This can be a (character) name of a binary file where the matrix is stored}

\item{tst}{(integer vector) Which elements from vector \code{y} (stored in \code{fm$y}) are in testing set and to plot. They must be contained in \code{fm$tst}. Default \code{tst=NULL} will consider the whole vector \code{fm$tst} to plot}

\item{group}{(data.frame) Column grouping for the individuals. The rows must match with the rows in \code{G} matrix}

\item{df}{(numeric) Average number of training individuals contributing to the prediction (active) of testing individuals. Default \code{df=NULL} will use the \code{df} that yielded the optimal accuracy}

\item{main}{(character/expression) Title of the plot}

\item{bg.color}{(character) Plot background color}

\item{line.color}{(character) Color of lines connecting nodes in rows with those in columns}

\item{line.tick}{(numeric) Tick of lines connecting nodes in rows with those in columns}

\item{curve}{\code{TRUE} or \code{FALSE} to whether draw curve lines connecting nodes in rows with those in columns}

\item{show.names}{\code{TRUE} or \code{FALSE} to whether show node names given by the row/column names of either \code{K} or \code{object} (when this is a matrix)}

\item{set.color}{(character vector) Color point of each type of node: row, 'active' column, and 'non-active' column, respectively}

\item{set.size}{(numeric vector) Size of each type of node: row, 'active' column, and 'non-active' column, respectively}

\item{group.shape}{(integer vector) Shape of each level of the grouping column provided as \code{group}}

\item{axis.labels}{\code{TRUE} or \code{FALSE} to whether show labels in both axes}

\item{unified}{\code{TRUE} or \code{FALSE} to whether show an unified plot or separated for each individual in 'testing'}

\item{point.color}{(character) Color of the points in the plot}

\item{ntst}{(integer) Maximum number of row nodes ('testing') that are plotted separated as indicated by \code{unified=FALSE}}

\item{legend.pos}{(character) Either "right", topright","bottomleft","bottomright","topleft", or "none" indicating where the legend is positioned in the plot}

\item{sets}{(character vector) Names of the types of node: row, 'active' column, and 'non-active' column, respectively}

\item{eps}{(numeric) Minumum value to declare nodes to be connected. Default is the machine precision (numerical zero)}
}
\value{
Returns the top-2 PC's plot connecting testing (predicted) individuals with training (predictors) individuals
}
\description{
Network plot of testing and training individuals from an object of the class 'SSI'
}
\details{
Plot edges in the plane xy given by a numerical matrix. When the object is a 'SSI' object the edges are taken from the regression coefficients from the 'SSI'. Edges are plotted according to the Fruchterman-Reingold algorithm. When a matrix \ifelse{html}{\out{<b>K</b>}}{\eqn{\textbf{K}}{K}} is provided, edges are plotted according to the spectral value decomposition of \ifelse{html}{\out{<b>Z K Z'</b> = <b>U D U'</b>}}{\eqn{\textbf{Z K Z'} = \textbf{U D U'}}{Z K Z' = U D U'}}
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  index = which(Y$CV == 1)              # Use only a subset of data
  M = scale(M[index,])/sqrt(ncol(M))    # Subset and scale markers
  G = tcrossprod(M)                     # Genomic relationship matrix
  y = as.vector(scale(Y[index,"E1"]))   # Subset ans scale response variable
  
  # Training and testing sets
  tst = seq(20,length(y),by=10)
  trn = seq_along(y)[-tst]

  fm = SSI(y,K=G,tst=tst,trn=trn)
  
  # Basic setting
  net.plot(fm)
  net.plot(fm, show.names=c(TRUE,TRUE,FALSE), set.size=c(4,2,1), df=10)
  
  \donttest{
  # Passing a matrix of coefficients
  B = as.matrix(coef(fm,df=15))
  net.plot(B, curve=TRUE, set.size=c(3.5,1.5,1))
  }
}
\author{
Marco Lopez-Cruz (\email{maraloc@gmail.com}) and Gustavo de los Campos
}

