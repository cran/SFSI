\name{collect}
\alias{collect}
\title{collect function}
\usage{
collect(prefix = "")
}
\arguments{
\item{prefix}{(character) Prefix that was added to the output files name, this may include a path}
}
\value{
An object of the class 'SSI' for which methods \code{fitted}, \code{plot} and \code{summary} exist
}
\description{
Collects all outputs saved at the provided \code{saveAt} parameter from the SSI analysis when testing data was splited
according to argument \code{subset}.
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  index = which(Y$CV == 1)             # Use only a subset of data
  M = scale(M[index,])/sqrt(ncol(M))   # Subset and scale markers
  G = tcrossprod(M)                    # Genomic relationship matrix
  y = as.vector(scale(Y[index,"E1"]))  # Subset response variable

  prefix <- paste0(tempdir(),"/testSSI")
  
  # Run the analysis into 4 subsets and save them at a given prefix
  fm <- SSI(y,K=G,tst=1:80,trn=81:length(y),subset=c(1,4),save.at=prefix)
  fm <- SSI(y,K=G,tst=1:80,trn=81:length(y),subset=c(2,4),save.at=prefix)
  fm <- SSI(y,K=G,tst=1:80,trn=81:length(y),subset=c(3,4),save.at=prefix)
  fm <- SSI(y,K=G,tst=1:80,trn=81:length(y),subset=c(4,4),save.at=prefix)

  # Collect all results after completion
  fm <- collect(prefix)
}

\author{
Marco Lopez-Cruz (\email{maraloc@gmail.com}) and Gustavo de los Campos
}
