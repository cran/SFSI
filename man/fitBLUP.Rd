
\name{7. BLUP estimation from Linear Mixed Model}
\alias{fitBLUP}
\title{Fitting a Linear Mixed model to calculate BLUP}
\usage{
fitBLUP(y, X = NULL, Z = NULL, K = NULL, 
        U = NULL, d = NULL, varU = NULL, 
        varE = NULL, intercept = TRUE, BLUP = TRUE, 
        method = c("REML","ML"), 
        interval = c(1E-9,1E9), tol = 1E-8, 
        maxiter = 1000, n.regions = 10,
        verbose = TRUE)
            
}
\arguments{
\item{y}{(numeric matrix) Response variable. It can contain >1 columns, each of them will be fitted separately}

\item{X}{(numeric matrix) Design matrix for fixed effects. When \code{X=NULL} a vector of ones is constructed only for the intercept (default)}

\item{Z}{(numeric matrix) Design matrix for random effects. When \code{Z=NULL} an identity matrix is considered (default) thus \ifelse{html}{\out{<b>G</b> = <b>K</b>}}{\eqn{\textbf{G = K}}{G = K}}; otherwise \ifelse{html}{\out{<b>G</b> = <b>Z K Z'</b>}}{\eqn{\textbf{G} = \textbf{Z K Z'}}{G = Z K Z'}} is used}

\item{K}{(numeric matrix) Kinship relationship matrix}

\item{U, d}{(numeric matrix/vector) Eigenvectors and eigenvalues from eigen value decomposition of \ifelse{html}{\out{<b>G</b> = <b>U D U'</b>}}{\eqn{\textbf{G} = \textbf{U D U'}}{G = U D U'}}}

\item{varU, varE}{(numeric) Genetic and residual variances. When both \code{varU} and \code{varE} are not \code{NULL} they are not calculated; otherwise, the likelihood function (REML or ML) is optimized to search for the genetic/residual variances ratio}

\item{intercept}{\code{TRUE} or \code{FALSE} to whether fit an intercept. When \code{FALSE}, the model assumes a null intercept}

\item{BLUP}{\code{TRUE} or \code{FALSE} to whether return the random effects estimates}

\item{method}{(character) Either 'REML' (Restricted Maximum Likelihood) or 'ML' (Maximum Likelihood)}

\item{tol}{(numeric) Maximum error between two consecutive solutions (convergence tolerance) when finding the root of the log-likelihood's first derivative}

\item{maxiter}{(integer) Maximum number of iterations to run before convergence is reached}

\item{interval}{(numeric vector) Range of values in which the root is searched}

\item{n.regions}{(numeric) Number of regions in which the searched 'interval' is divided for local optimization}

\item{verbose}{\code{TRUE} or \code{FALSE} to whether show progress}
}
\value{
Returns a list object that contains the elements:
\itemize{
  \item \code{b}: (vector) fixed effects solutions (including the intercept).
  \item \code{u}: (vector) random effects solutions.
  \item \code{varU}: random effect variance.
  \item \code{varE}: residual variance.
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

\ifelse{html}{\out{<p style='text-align:center'><b>y</b> = <b>X b</b> + <b>Z u</b> + <b>e</b></p>}}{\deqn{\textbf{y}=\textbf{X}\textbf{b}+\textbf{Z}\textbf{u}+\textbf{e}}{y = X b + Z u + e}}

where
\ifelse{html}{\out{<b>y</b>}}{\eqn{\textbf{y}}{y}} is a vector with the response,
\ifelse{html}{\out{<b>b</b>}}{\eqn{\textbf{b}}{b}} is the vector of fixed effects,
\ifelse{html}{\out{<b>u</b>}}{\eqn{\textbf{u}}{u}} is the vector of the (random) genetic values of the genotypes,
\ifelse{html}{\out{<b>e</b>}}{\eqn{\textbf{e}}{e}} is the vector of environmental residuals (random error), and
\ifelse{html}{\out{<b>X</b>}}{\eqn{\textbf{X}}{X}} and \ifelse{html}{\out{<b>Z</b>}}{\eqn{\textbf{Z}}{Z}} are design matrices conecting the fixed and genetic effects with replicates. Genetic values are assumed to follow a Normal distribution as
\ifelse{html}{\out{<b>u</b> ~ N(<b>0</b>,&sigma;<sup>2</sup><sub>u</sub><b>K</b>)}}{\eqn{\textbf{u}\sim N(\textbf{0},\sigma^2_u\textbf{K})}{u ~ N(0,sigma^2_u K)}}, and the error terms are assumed
\ifelse{html}{\out{<b>e</b> ~ N(<b>0</b>,&sigma;<sup>2</sup><sub>e</sub><b>I</b>)}}{\eqn{\textbf{e}\sim N(\textbf{0},\sigma^2_e\textbf{I})}{e ~ N(0,sigma^2_e I)}}.

The vector of genetic values
\ifelse{html}{\out{<b>g</b> = <b>Z u</b>}}{\eqn{\textbf{g}=\textbf{Z}\textbf{u}}{g = Z u}} will therefore follow 
\ifelse{html}{\out{<b>g</b> ~ N(<b>0</b>,&sigma;<sup>2</sup><sub>u</sub><b>G</b>)}}{\eqn{\textbf{g}\sim N(\textbf{0},\sigma^2_u\textbf{G})}{g ~ N(0,sigma^2_u G)}} where
\ifelse{html}{\out{<b>G</b> = <b>Z K Z'</b>}}{\eqn{\textbf{G}=\textbf{Z}\textbf{K}\textbf{Z}'}{G = Z K Z'}}.
In the un-replicated case, \ifelse{html}{\out{<b>Z</b> = <b>I</b>}}{\eqn{\textbf{Z}=\textbf{I}}{Z = I}} is an identity matrix, and hence 
\ifelse{html}{\out{<b>g</b> = <b>u</b>}}{\eqn{\textbf{g}=\textbf{u}}{g = u}} and
\ifelse{html}{\out{<b>G</b> = <b>K</b>}}{\eqn{\textbf{G}=\textbf{K}}{G = K}}.

The predicted values \ifelse{html}{\out{<b>u</b><sub>trn</sub> = (u<sub>i</sub>)}}{\eqn{\textbf{u}_{trn}=(u_i)}{u_trn = (u_i)}},
\ifelse{html}{\out{i = 1,2,...,n<sub>trn</sub>}}{\eqn{i=1,2,...,n_{trn}}{i = 1,2,...,n_trn}}, corresponding to observed data (training set) are derived as

\ifelse{html}{\out{<p style='text-align:center'><b>u</b><sub>trn</sub> = <b>B</b> (<b>y</b><sub>trn</sub> - <b>X</b><sub>trn</sub><b>b</b>)</p>}}{\deqn{\textbf{u}_{tst}=\textbf{B}(\textbf{y}_{trn}-\textbf{X}_{trn}\textbf{b})}{u_tst = B (y_trn - X_trn*b)}}

where \ifelse{html}{\out{<b>B</b>}}{\eqn{\textbf{B}}{B}}
is a matrix of weights given by

\ifelse{html}{\out{<p style='text-align:center'><b>B</b> = &sigma;<sup>2</sup><sub>u</sub><b>G</b><sub>trn</sub> (&sigma;<sup>2</sup><sub>u</sub><b>G</b><sub>trn</sub> + &sigma;<sup>2</sup><sub>e</sub><b>I</b>)<sup>-1</sup></p>}}{\deqn{\textbf{B}=\sigma_u^2\textbf{G}_{trn}(\sigma_u^2\textbf{G}_{trn} + \sigma_e^2\textbf{I})^{-1}}{B = sigma^2_u G[trn](sigma^2_u G[trn] + sigma^2_e I)^-1}}

where \ifelse{html}{\out{<b>G</b><sub>trn</sub>}}{\eqn{\textbf{G}_{trn}}{G[trn]}}
is the sub-matrix corresponding to the training set. This matrix can be rewritten as

\ifelse{html}{\out{<p style='text-align:center'><b>B</b> = <b>G</b><sub>trn</sub> (<b>G</b><sub>trn</sub> + &theta;<b>I</b>)<sup>-1</sup></p>}}{\deqn{\textbf{B}=\textbf{G}_{trn}(\textbf{G}_{trn} + \theta\textbf{I})^{-1}}{B = G[trn](G[trn] + theta I)^-1}}

where \ifelse{html}{\out{&theta; = &sigma;<sup>2</sup><sub>e</sub>/&sigma;<sup>2</sup><sub>u</sub>}}{\eqn{\theta=\sigma_e^2/\sigma_u^2}{theta = sigma^2_e/sigma^2_u}} is the residual/genetic variances ratio representing a ridge-like shrinkage parameter.

The matrix \ifelse{html}{\out{<b>H</b> = <b>G</b><sub>trn</sub> + &theta;<b>I</b>}}{\eqn{\textbf{H}=\textbf{G}_{trn} + \theta\textbf{I}}{H = G[trn] + theta I}}
in the above equation can be used to obtain predictions corresponding to un-observed data (testing set), \ifelse{html}{\out{<b>u</b><sub>tst</sub>}}{\eqn{\textbf{u}_{tst}}{u_tst}}, by  

\ifelse{html}{\out{<p style='text-align:center'><b>B</b> = <b>G</b><sub>tst,trn</sub> (<b>G</b><sub>trn</sub> + &theta;<b>I</b>)<sup>-1</sup></p>}}{\deqn{\textbf{B}=\textbf{G}_{tst,trn}(\textbf{G}_{trn} + \theta\textbf{I})^{-1}}{B = G[tst,trn](G[trn] + theta I)^-1}}

where
\ifelse{html}{\out{<b>G</b><sub>tst,trn</sub>}}{\eqn{\textbf{G}_{tst,trn}}{G[tst,trn]}}
is the sub-matrix of \ifelse{html}{\out{<b>G</b>}}{\eqn{\textbf{G}}{G}} corresponding to the testing set (in rows) and training set (in columns).  

Solutions are found using the GEMMA (Genome-wide Efficient Mixed Model Analysis) approach (Zhou & Stephens, 2012). First, the Brent's method is implemented to solve for the genetic/residual variances ratio (i.e., \ifelse{html}{\out{1/&theta;}}{\eqn{1/\theta}{1/theta}}) from the first derivative of the log-likelihood (either REML or ML). Then, variances \ifelse{html}{\out{&sigma;<sup>2</sup><sub>u</sub>}}{\eqn{\sigma_u^2}{sigma^2_u}} and \ifelse{html}{\out{&sigma;<sup>2</sup><sub>e</sub>}}{\eqn{\sigma_e^2}{sigma^2_e}} are calculated. Finally, \ifelse{html}{\out{<b>b</b>}}{\eqn{\textbf{b}}{b}} is obtained using Generalized Least Squares.
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  index = which(Y$trial \%in\% 1:10)     # Use only a subset of data
  Y = Y[index,]
  M = scale(M[index,])/sqrt(ncol(M))  # Subset and scale markers
  G = tcrossprod(M)                   # Genomic relationship matrix
  y = as.vector(scale(Y[,"E1"]))      # Scale response variable
  
  # Training and testing sets
  tst = which(Y$trial \%in\% 1:2)
  trn = which(!Y$trial \%in\% 1:2)

  yNA <- y
  yNA[tst] <- NA
  fm = fitBLUP(yNA, K=G)
  plot(y[tst],fm$u[tst])    # Predicted vs observed values in testing set
  cor(y[tst],fm$u[tst])     # Prediction accuracy in testing set
  cor(y[trn],fm$u[trn])     # Prediction accuracy in training set
  fm$varU                   # Genetic variance
  fm$varE                   # Residual variance
  fm$h2                     # Heritability
  fm$b                      # Intercept
  
}
\references{
VanRaden PM (2008). Efficient methods to compute genomic predictions. \emph{Journal of Dairy Science}, \bold{91}(11), 4414–4423.

Zhou X, Stephens M (2012). Genome-wide efficient mixed-model analysis for association studies. \emph{Nature Genetics}, \bold{44}(7), 821-824
}
\keyword{fitBLUP}
