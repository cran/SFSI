
\name{Collecting SSI outputs}
\alias{collect}
\title{Collecting SSI outputs}
\usage{
collect(path = "")
}
\arguments{
\item{path}{(character) Path where output files were saved, this may include a prefix}
}
\value{
An object of the class 'SSI' for which methods \code{fitted}, \code{plot} and \code{summary} exist
}
\description{
Collects all outputs saved at the provided \code{save.at} parameter from the SSI analysis when testing data was splited
according to argument \code{subset}.
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  index = which(Y$trial \%in\% 1:8)      # Use only a subset of data
  Y = Y[index,]
  M = scale(M[index,])/sqrt(ncol(M))   # Subset and scale markers
  G = tcrossprod(M)                    # Genomic relationship matrix
  y = as.vector(scale(Y[,"E1"]))       # Response variable

  # Sets (testing=0, training=1)
  trn_tst = ifelse(Y$trial \%in\% 1:2, 0, 1)
  
  path = paste0(tempdir(),"/testSSI_")
  
  # Run the analysis into 4 subsets and save them at a given path
  fm = SSI(y,K=G,trn_tst=trn_tst,subset=c(1,4),save.at=path)
  fm = SSI(y,K=G,trn_tst=trn_tst,subset=c(2,4),save.at=path)
  fm = SSI(y,K=G,trn_tst=trn_tst,subset=c(3,4),save.at=path)
  fm = SSI(y,K=G,trn_tst=trn_tst,subset=c(4,4),save.at=path)

  # Collect all results after completion
  fm = collect(path)
}

