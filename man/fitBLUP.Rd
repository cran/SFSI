\name{fitBLUP}
\alias{fitBLUP}
\title{Function fitBLUP}
\usage{
fitBLUP(y, X = NULL, Z = NULL, K = NULL, U = NULL, 
        d = NULL, h2 = NULL, BLUP = TRUE, 
        method = c("REML","ML"), return.Hinv = FALSE, 
        tol = 1E-5, maxIter = 1000, interval = c(1E-9,1E9),
        warn = TRUE)  
}
\arguments{
\item{y}{(numeric vector) Response variable}

\item{X}{(numeric matrix) Design matrix for fixed effects. When \code{X=NULL} a vector of ones is constructed only for the intercept (default)}

\item{Z}{(numeric matrix) Design matrix for random effects. When \code{Z=NULL} an identity matrix is considered (default) thus \ifelse{html}{\out{<b>G</b> = <b>K</b>}}{\eqn{\textbf{G = K}}{G = K}}; otherwise \ifelse{html}{\out{<b>G</b> = <b>Z K Z'</b>}}{\eqn{\textbf{G} = \textbf{Z K Z'}}{G = Z K Z'}} is used}

\item{K}{(numeric matrix) Kinship relationships. This can be of the "float32" type as per the 'float' R-package, or a (character) name of a binary file where the matrix is stored}

\item{U}{(numeric matrix) Eigenvectors from spectral value decomposition of \ifelse{html}{\out{<b>G</b> = <b>U D U'</b>}}{\eqn{\textbf{G} = \textbf{U D U'}}{G = U D U'}}}

\item{d}{(numeric vector) Eigenvalues from spectral value decomposition of \ifelse{html}{\out{<b>G</b> = <b>U D U'</b>}}{\eqn{\textbf{G} = \textbf{U D U'}}{G = U D U'}}}

\item{h2}{(numeric) Heritability of the response variable. When it is not \code{NULL}, the variance components are not estimated and only fixed and random effects are computed}

\item{BLUP}{\code{TRUE} or \code{FALSE} to whether return the random effects estimates}

\item{method}{(character) Either 'REML' (Restricted Maximum Likelihood) or 'ML' (Maximum Likelihood)}

\item{return.Hinv}{\code{TRUE} or \code{FALSE} to whether return the inverse of the matrix H}

\item{tol}{(numeric) Maximum error between two consecutive solutions (convergence tolerance) when finding the root of the log-likelihood's first derivative}

\item{maxIter}{(integer) Maximum number of iterations to run before convergence is reached}

\item{interval}{(numeric vector) Range of values in which the root is searched}

\item{warn}{\code{TRUE} or \code{FALSE} to whether show warnings}
}
\value{
Returns a list object that contains the elements:
\itemize{
  \item \code{b}: (vector) fixed effects solutions (including the intercept).
  \item \code{u}: (vector) random effects solutions.
  \item \code{varU}: random effect variance.
  \item \code{varE}: error variance.
  \item \code{h2}: heritability.
  \item \code{convergence}: (logical) whether Brent's method converged.
  \item \code{method}: either 'REML' or 'ML' method used.
}
}
\description{
Solves the Linear Mixed Model and calculates the Best Linear Unbiased Predictor (BLUP)
}
\details{
The basic linear mixed model that relates phenotypes with genetic values is of the form

\ifelse{html}{\out{<center><b>y</b> = <b>X b</b> + <b>Z u</b> + <b>e</b></center>}}{\deqn{\textbf{y}=\textbf{X}\textbf{b}+\textbf{Z}\textbf{u}+\textbf{e}}{y = X b + Z u + e}}

where
\ifelse{html}{\out{<b>y</b>}}{\eqn{\textbf{y}}{y}} is a vector with the response,
\ifelse{html}{\out{<b>b</b>}}{\eqn{\textbf{b}}{b}} is the vector of fixed effects,
\ifelse{html}{\out{<b>u</b>}}{\eqn{\textbf{u}}{u}} is the vector of the (random) genetic values of the genotypes,
\ifelse{html}{\out{<b>e</b>}}{\eqn{\textbf{e}}{e}} is the vector of environmental residuals (random error), and
\ifelse{html}{\out{<b>X</b>}}{\eqn{\textbf{X}}{X}} and \ifelse{html}{\out{<b>Z</b>}}{\eqn{\textbf{Z}}{Z}} are design matrices conecting the fixed and genetic effects with replicates. Genetic values are assumed to follow a Normal distribution as
\ifelse{html}{\out{<b>u</b> ~ N(<b>0</b>,&sigma;<sup>2</sup><sub>u</sub><b>K</b>)}}{\eqn{\textbf{u}\sim N(\textbf{0},\sigma^2_u\textbf{K})}{u ~ N(0,sigma^2_u K)}}, and the error terms are assumed
\ifelse{html}{\out{<b>e</b> ~ N(<b>0</b>,&sigma;<sup>2</sup><sub>e</sub><b>D</b>)}}{\eqn{\textbf{e}\sim N(\textbf{0},\sigma^2_e\textbf{D})}{e ~ N(0,sigma^2_e D)}}, with \ifelse{html}{\out{<b>D</b>=<b>I</b>}}{\eqn{\textbf{D}=\textbf{I}}{D = I}} being an identity matrix.

The resulting vector of genetic values
\ifelse{html}{\out{<b>g</b> = <b>Z u</b>}}{\eqn{\textbf{g}=\textbf{Z}\textbf{u}}{g = Z u}} will therefore follow 
\ifelse{html}{\out{<b>g</b> ~ N(<b>0</b>,&sigma;<sup>2</sup><sub>u</sub><b>G</b>)}}{\eqn{\textbf{g}\sim N(\textbf{0},\sigma^2_u\textbf{G})}{g ~ N(0,sigma^2_u G)}} where
\ifelse{html}{\out{<b>G</b> = <b>Z K Z'</b>}}{\eqn{\textbf{G}=\textbf{Z}\textbf{K}\textbf{Z}'}{G = Z K Z'}}.
In the un-replicated case, \ifelse{html}{\out{<b>Z</b> = <b>I</b>}}{\eqn{\textbf{Z}=\textbf{I}}{Z = I}} is an identity matrix, and hence 
\ifelse{html}{\out{<b>g</b> = <b>u</b>}}{\eqn{\textbf{g}=\textbf{u}}{g = u}} and
\ifelse{html}{\out{<b>G</b> = <b>K</b>}}{\eqn{\textbf{G}=\textbf{K}}{G = K}}.

Only information from observed data (training set) is used to derive predictions of the values \ifelse{html}{\out{<b>u</b><sub>trn</sub> = (u<sub>i</sub>)}}{\eqn{\textbf{u}_{trn}=(u_i)}{u_trn = (u_i)}},
\ifelse{html}{\out{i = 1,2,...,n<sub>trn</sub>}}{\eqn{i=1,2,...,n_{trn}}{i = 1,2,...,n_trn}}, as

\ifelse{html}{\out{<center><b>u</b><sub>trn</sub> = <b>H</b> (<b>y</b><sub>trn</sub> - <b>X</b><sub>trn</sub><b>b</b>)</center>}}{\deqn{\textbf{u}_{tst}=\textbf{H}(\textbf{y}_{trn}-\textbf{X}_{trn}\textbf{b})}{u_tst = H (y_trn - X_trn*b)}}

where \ifelse{html}{\out{<b>H</b>}}{\eqn{\textbf{H}}{H}}
is a matrix of weights given by
\ifelse{html}{\out{<center><b>H</b> = <b>G</b><sub>trn</sub> (<b>G</b><sub>trn</sub> + &theta;<b>D</b>)<sup>-1</sup></center>}}{\deqn{\textbf{H}=\textbf{G}_{trn}(\textbf{G}_{trn} + \theta\textbf{D})^{-1}}{H = G[trn](G[trn] + theta D)^-1}}

where \ifelse{html}{\out{<b>G</b><sub>trn</sub>}}{\eqn{\textbf{G}_{trn}}{G[trn]}}
is the sub-matrix corresponding to the training set, and \ifelse{html}{\out{&theta; = &sigma;<sup>2</sup><sub>e</sub>/&sigma;<sup>2</sup><sub>u</sub>}}{\eqn{\theta=\sigma_e^2/\sigma_u^2}{theta = sigma^2_e/sigma^2_u}} is the variance components ratio representing a shrinkage parameter. This parameter is expressed in terms of the heritability, \ifelse{html}{\out{h<sup>2</sup> = &sigma;<sup>2</sup><sub>u</sub>/(&sigma;<sup>2</sup><sub>u</sub> + &sigma;<sup>2</sup><sub>e</sub>)}}{\eqn{h^2=\sigma^2_u/(\sigma^2_u+\sigma^2_e)}{h^2 = sigma^2_u/(sigma^2_u + sigma^2_e)}}, as \ifelse{html}{\out{&theta; = (1 - h<sup>2</sup>)/h<sup>2</sup>}}{\eqn{\theta=(1-h^2)/h^2}{theta = (1 - h^2)/h^2}}.

The prediction of values \ifelse{html}{\out{<b>u</b><sub>tst</sub>}}{\eqn{\textbf{u}_{tst}}{u_tst}} corresponding to un-observed data (testing set) can be done by using 
\ifelse{html}{\out{<center><b>H</b> = <b>G</b><sub>tst,trn</sub> (<b>G</b><sub>trn</sub> + &lambda;<sub>0</sub><b>D</b>)<sup>-1</sup></center>}}{\deqn{\textbf{H}=\textbf{G}_{tst,trn}(\textbf{G}_{trn} + \lambda_0\textbf{D})^{-1}}{H = G[tst,trn](G[trn] + lambda_0 D)^-1}}

where
\ifelse{html}{\out{<b>G</b><sub>tst,trn</sub>}}{\eqn{\textbf{G}_{tst,trn}}{G[tst,trn]}}
is the sub-matrix of \ifelse{html}{\out{<b>G</b>}}{\eqn{\textbf{G}}{G}} corresponding to the testing set (in rows) and training set (in columns).  

Solutions are found using the GEMMA (Genome-wide Efficient Mixed Model Analysis) approach (Zhou & Stephens, 2012). GEMMA implements the Brent's method to efficiently optimize the first derivative of the log-likelihood to solve for the variance components ratio.
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  index = which(Y$CV \%in\% 1:2)       # Use only a subset of data
  M = scale(M[index,])/sqrt(ncol(M))   # Subset and scale markers
  G = tcrossprod(M)                    # Genomic relationship matrix
  y = as.vector(scale(Y[index,"E1"]))  # Subset response variable

  # Training and testing sets
  tst = seq(1,length(y),by=3)
  trn = seq_along(y)[-tst]

  yNA <- y
  yNA[tst] <- NA
  fm1 = fitBLUP(yNA, K=G)
  plot(y[tst],fm1$u[tst])      # Predicted vs observed values in testing set
  cor(y[tst],fm1$u[tst])       # Prediction accuracy in testing set
  cor(y[trn],fm1$u[trn])       # Prediction accuracy in training set
  fm1$h2                       # Heritability
  
  \donttest{
  if(requireNamespace("float")){
   # Using a 'float' type variable
   G2 = float::fl(G)
   fm2 = fitBLUP(yNA, K=G2)
   max(abs(fm1$u-fm2$u))  # Check for discrepances
  }
  }

}
\references{
VanRaden PM (2008). Efficient methods to compute genomic predictions. \emph{Journal of Dairy Science}, \bold{91}(11), 4414–4423.

Zhou X, Stephens M (2012). Genome-wide efficient mixed-model analysis for association studies. \emph{Nature Genetics}, \bold{44}(7), 821-824
}
\author{
Paulino Perez-Rodriguez, Marco Lopez-Cruz (\email{maraloc@gmail.com}) and Gustavo de los Campos
}
\keyword{fitBLUP}
