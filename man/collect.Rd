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
  
  index = which(Y$trial \%in\% 1:6)      # Use only a subset of data
  Y = Y[index,]
  M = scale(M[index,])/sqrt(ncol(M))   # Subset and scale markers
  G = tcrossprod(M)                    # Genomic relationship matrix
  y = as.vector(scale(Y[,"E1"]))       # Response variable

  # Training and testing sets
  tst = which(Y$trial == 2)
  trn = which(Y$trial != 2)
  
  prefix <- paste0(tempdir(),"/testSSI")
  
  # Run the analysis into 4 subsets and save them at a given prefix
  fm <- SSI(y,K=G,tst=tst,trn=trn,subset=c(1,4),save.at=prefix)
  fm <- SSI(y,K=G,tst=tst,trn=trn,subset=c(2,4),save.at=prefix)
  fm <- SSI(y,K=G,tst=tst,trn=trn,subset=c(3,4),save.at=prefix)
  fm <- SSI(y,K=G,tst=tst,trn=trn,subset=c(4,4),save.at=prefix)

  # Collect all results after completion
  fm <- collect(prefix)
}

\author{
Marco Lopez-Cruz (\email{maraloc@gmail.com}) and Gustavo de los Campos
}
