
\name{getGenCov}
\alias{getGenCov}
\title{Genetic covariances}
\usage{
getGenCov(y1, y2, X = NULL, Z = NULL, K = NULL, 
          U = NULL, d = NULL, scale = TRUE, 
          mc.cores = 1, warn = FALSE, ...)
}
\arguments{
\item{y1}{(numeric vector) Response variable 1}

\item{y2}{(numeric matrix) Response variable 2. The number of rows must be equal to length of vector y1}

\item{X}{(numeric matrix) Design matrix for fixed effects. When \code{X=NULL} a vector of ones is constructed only for the intercept (default)}

\item{Z}{(numeric matrix) Design matrix for random effects. When \code{Z=NULL} an identity matrix is considered (default) thus \ifelse{html}{\out{<b>G</b> = <b>K</b>}}{\eqn{\textbf{G = K}}{G = K}}; otherwise \ifelse{html}{\out{<b>G</b> = <b>Z K Z'</b>}}{\eqn{\textbf{G} = \textbf{Z K Z'}}{G = Z K Z'}} is used}

\item{K}{(numeric matrix) Kinship relationships}

\item{U}{(numeric matrix) Eigenvectors from spectral value decomposition of \ifelse{html}{\out{<b>G</b> = <b>U D U'</b>}}{\eqn{\textbf{G} = \textbf{U D U'}}{G = U D U'}}}

\item{d}{(numeric vector) Eigenvalues from spectral value decomposition of \ifelse{html}{\out{<b>G</b> = <b>U D U'</b>}}{\eqn{\textbf{G} = \textbf{U D U'}}{G = U D U'}}}

\item{scale}{\code{TRUE} or \code{FALSE} to scale y1 and y2 by their corresponding standard deviations so the resulting variables will have unit variance}

\item{mc.cores}{(integer) Number of cores used. The analysis is run in parallel when \code{mc.cores} is greater than 1. Default is \code{mc.cores=1}}

\item{warn}{\code{TRUE} or \code{FALSE} to whether show warnings}

\item{...}{Other arguments passed to the function 'fitBLUP'} 
}
\value{
Returns a list object that contains the elements:
\itemize{
  \item \code{varU1}: genetic variance for response variable 1.
  \item \code{varU2}: (vector) genetic variances for response variable 2.
  \item \code{varE1}: error variance for response variable 1.
  \item \code{varE2}: (vector) error variances for response variable 2.
  \item \code{covU}: (vector) genetic covariances between response variables 1 and 2.
  \item \code{covE}: (vector) environmental covariances between response variables 1 and 2.
}
}
\description{
Pairwise genetic covariances for variables with the same experimental design and equal variance
}
\details{
Assumes that both \ifelse{html}{\out{<b>y</b><sub>1</sub>}}{\eqn{\textbf{y}_1}{y_1}} and
\ifelse{html}{\out{<b>y</b><sub>2</sub>}}{\eqn{\textbf{y}_2}{y_2}} follow the basic linear mixed model that relates phenotypes with genetic values of the form

\ifelse{html}{\out{<center><b>y</b><sub>1</sub> = <b>X b</b><sub>1</sub> + <b>Z u</b><sub>1</sub> + <b>e</b><sub>1</sub></center>}}{\deqn{\textbf{y}_1=\textbf{X}\textbf{b}_1+\textbf{Z}\textbf{u}_1+\textbf{e}_1}{y = X b_1 + Z u_1 + e_1}}

\ifelse{html}{\out{<center><b>y</b><sub>2</sub> = <b>X b</b><sub>2</sub> + <b>Z u</b><sub>2</sub> + <b>e</b><sub>2</sub></center>}}{\deqn{\textbf{y}_2=\textbf{X}\textbf{b}_2+\textbf{Z}\textbf{u}_2+\textbf{e}_2}{y = X b_2 + Z u_2 + e_2}}

where
\ifelse{html}{\out{<b>b</b><sub>1</sub>}}{\eqn{\textbf{b}_1}{b_1}} and 
\ifelse{html}{\out{<b>b</b><sub>2</sub>}}{\eqn{\textbf{b}_2}{b_2}} are the specific fixed effects,
\ifelse{html}{\out{<b>u</b><sub>1</sub>}}{\eqn{\textbf{u}_1}{u_1}} and 
\ifelse{html}{\out{<b>u</b><sub>2</sub>}}{\eqn{\textbf{u}_2}{u_2}} are the specific genetic values of the genotypes,
\ifelse{html}{\out{<b>e</b><sub>1</sub>}}{\eqn{\textbf{e}_1}{e_1}} and 
\ifelse{html}{\out{<b>e</b><sub>2</sub>}}{\eqn{\textbf{e}_2}{e_2}} are the vectors of specific environmental residuals, and
\ifelse{html}{\out{<b>X</b>}}{\eqn{\textbf{X}}{X}} and \ifelse{html}{\out{<b>Z</b>}}{\eqn{\textbf{Z}}{Z}} are common design matrices conecting the fixed and genetic effects with replicates. Genetic values are assumed to follow a Normal distribution as
\ifelse{html}{\out{<b>u</b><sub>1</sub> ~ N(<b>0</b>,&sigma;<sup>2</sup><sub>u1</sub><b>K</b>)}}{\eqn{\textbf{u}_1\sim N(\textbf{0},\sigma^2_{u_1}\textbf{K})}{u_1 ~ N(0,sigma^2_u1 K)}} and 
\ifelse{html}{\out{<b>u</b><sub>2</sub> ~ N(<b>0</b>,&sigma;<sup>2</sup><sub>u2</sub><b>K</b>)}}{\eqn{\textbf{u}_2\sim N(\textbf{0},\sigma^2_{u_2}\textbf{K})}{u_2 ~ N(0,sigma^2_u2 K)}}, and environmental terms are assumed
\ifelse{html}{\out{<b>e</b><sub>1</sub> ~ N(<b>0</b>,&sigma;<sup>2</sup><sub>e1</sub><b>I</b>)}}{\eqn{\textbf{e}_1\sim N(\textbf{0},\sigma^2_{e_1}\textbf{I})}{e_1 ~ N(0,sigma^2_e1 I)}} and
\ifelse{html}{\out{<b>e</b><sub>2</sub> ~ N(<b>0</b>,&sigma;<sup>2</sup><sub>e2</sub><b>I</b>)}}{\eqn{\textbf{e}_2\sim N(\textbf{0},\sigma^2_{e_2}\textbf{I})}{u_2 ~ N(0,sigma^2_e2 I)}}.

The genetic covariance 
\ifelse{html}{\out{&sigma;<sup>2</sup><sub>u1,u2</sub>}}{\eqn{\sigma^2_{u_1 u_2}}{sigma^2_u1u2)}}
is estimated from the formula for the variance for the sum of two variables as 

\ifelse{html}{\out{<center>&sigma;<sup>2</sup><sub>u1,u2</sub> = 1/2(&sigma;<sup>2</sup><sub>u3</sub> - &sigma;<sup>2</sup><sub>u1</sub> - &sigma;<sup>2</sup><sub>u2</sub>)</center>}}{\deqn{\sigma^2_{u_1 u_2}=\frac{1}{2}(\sigma^2_{u_3}-\sigma^2_{u_1}-\sigma^2_{u_2})}{sigma^2_u1u2 = 1/2 (sigma^2_u3 - sigma^2_u1 - sigma^2_u2)}}

where \ifelse{html}{\out{&sigma;<sup>2</sup><sub>u3</sub>}}{\eqn{\sigma^2_{u_3}}{sigma^2_u3)}}
is the genetic variance of the variable 
\ifelse{html}{\out{<b>y</b><sub>3</sub> = <b>y</b><sub>1</sub> + <b>y</b><sub>2</sub>}}{\eqn{\textbf{y}_3=\textbf{y}_1+\textbf{y}_2}{y_3=y_1+y_2}} that also follows the same model as for 
\ifelse{html}{\out{<b>y</b><sub>1</sub>}}{\eqn{\textbf{y}_1}{y_1}} and
\ifelse{html}{\out{<b>y</b><sub>2</sub>}}{\eqn{\textbf{y}_2}{y_2}}.

Likewise, the environmental covariance 
\ifelse{html}{\out{&sigma;<sup>2</sup><sub>e1,e2</sub>}}{\eqn{\sigma^2_{e_1 e_2}}{sigma^2_e1e2)}}
is estimated as 

\ifelse{html}{\out{<center>&sigma;<sup>2</sup><sub>e1,e2</sub> = 1/2(&sigma;<sup>2</sup><sub>e3</sub> - &sigma;<sup>2</sup><sub>e1</sub> - &sigma;<sup>2</sup><sub>e2</sub>)</center>}}{\deqn{\sigma^2_{e_1 e_2}=\frac{1}{2}(\sigma^2_{e_3}-\sigma^2_{e_1}-\sigma^2_{e_2})}{sigma^2_e1e2 = 1/2 (sigma^2_e3 - sigma^2_e1 - sigma^2_e2)}}

where \ifelse{html}{\out{&sigma;<sup>2</sup><sub>e3</sub>}}{\eqn{\sigma^2_{e_3}}{sigma^2_e3)}}
is the error variance of the variable 
\ifelse{html}{\out{<b>y</b><sub>3</sub>}}{\eqn{\textbf{y}_3}{y_3}}.

Solutions are found using the function 'fitBLUP' (see \code{help(fitBLUP)}) to sequentialy fit mixed models for all the variables \ifelse{html}{\out{<b>y</b><sub>1</sub>}}{\eqn{\textbf{y}_1}{y_1}}, \ifelse{html}{\out{<b>y</b><sub>2</sub>}}{\eqn{\textbf{y}_2}{y_2}} and
\ifelse{html}{\out{<b>y</b><sub>3</sub>}}{\eqn{\textbf{y}_3}{y_3}}.
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  X = scale(X[1:200,])/sqrt(ncol(X))    # Subset and scale markers
  G = tcrossprod(X)                     # Genomic relationship matrix
  y = as.vector(scale(Y[1:200,"YLD"]))  # Subset response variable
  WL = scale(WL[1:200,20:40])           # Subset reflectance data
  
  fm = getGenCov(y,WL,K=G)
  
  covU = fm$covU                     # Genetic covariance
  covP_corrected = fm$covU+fm$covE   # Phenotypic covariance
  covP_uncorrected = cov(y,WL)       # Genetic covariance
  
  plot(covP_corrected,covP_uncorrected)
  plot(covU,covP_uncorrected)
  plot(covU,covP_corrected)
}
\author{
Marco Lopez-Cruz (\email{maraloc@gmail.com}) and Gustavo de los Campos
}
\keyword{fitBLUP}
