
\name{6. Sparse Selection Index (SSI)}
\alias{SSI}
\alias{SSI.CV}
\title{Sparse Selection Index}
\usage{
SSI(y, X = NULL, b = NULL, Z = NULL, K, trn_tst = NULL, 
    varU = NULL, varE = NULL, intercept = TRUE, 
    alpha = 1, lambda = NULL, nlambda = 100,
    lambda.min = .Machine$double.eps^0.5,
    common.lambda = TRUE, subset = NULL, tol = 1E-4, 
    maxiter = 500, method = c("REML","ML"), name = NULL, 
    save.at = NULL, mc.cores = 1L,
    precision.format = c("single","double"),
    verbose = TRUE)
    
SSI.CV(y, X = NULL, b = NULL, Z = NULL, K, trn_tst = NULL,
       varU = NULL, varE = NULL, intercept = TRUE, 
       alpha = 1, lambda = NULL, nlambda = 100, 
       lambda.min = .Machine$double.eps^0.5,
       common.lambda = TRUE, nCV = 1L, nfolds = 5, seed = NULL,
       tol = 1E-4, maxiter = 500, method = c("REML","ML"),
       name = NULL, mc.cores = 1L, verbose = TRUE)
}
\arguments{
\item{y}{(numeric matrix) Response variable. It can contain >1 columns for multi-trait analysis}

\item{X}{(numeric matrix) Design matrix for fixed effects. When \code{X=NULL}, a vector of ones is constructed only for the intercept (default)}

\item{b}{(numeric vector) Fixed effects. When \code{b=NULL}, only the intercept is estimated from training data using generalized least squares (default)}

\item{K}{(numeric matrix) Kinship relationship matrix}

\item{Z}{(numeric matrix) Design matrix for random effects. When \code{Z=NULL} an identity matrix is considered (default) thus \ifelse{html}{\out{<b>G</b> = <b>K</b>}}{\eqn{\textbf{G = K}}{G = K}}; otherwise \ifelse{html}{\out{<b>G</b> = <b>Z K Z'</b>}}{\eqn{\textbf{G} = \textbf{Z K Z'}}{G = Z K Z'}} is used}

\item{varU, varE}{(numeric) Genetic and residual variances. When either \code{varU=NULL} or \code{varE=NULL} (default), they are calculated from training data using the function 'fitBLUP' (see \code{help(fitBLUP)}). When \code{y} has >1 columns both \code{varU} and \code{varE} are matrices pairwise estimated (see \code{help(getGenCov)})}

\item{intercept}{\code{TRUE} or \code{FALSE} to whether fit an intercept. When \code{FALSE}, the model assumes a null intercept}

\item{trn_tst}{(integer vector) Which elements from matrix \code{y} are in training/testing set (0: testing, 1: training, NA are ignored). Default \code{trn_tst=NULL} will consider all individuals as training}

\item{subset}{(integer vector) \eqn{c(m,M)} to fit the model only for the \ifelse{html}{\out{m<sup>th</sup>}}{\eqn{m^{th}}{m^{th}}}
subset out of \eqn{M} subsets that the testing set will be divided into. Results can be automatically saved when \code{saveAt} argument
is provided and can be retrieved later using function 'collect' (see \code{help(collect)}). Default is \code{subset=NULL} for no subsetting, then the model is fitted for all testing data}

\item{alpha}{(numeric) Value between 0 and 1 for the weights given to the L1 and L2-penalties}

\item{lambda}{(numeric vector) Penalization parameter sequence. Default is \code{lambda=NULL}, in this case a decreasing grid of \code{nlambda} lambdas will be generated
starting from a maximum equal to

\ifelse{html}{\out{<p style='text-align:center;font-family:courier'>max(abs(G[trn,tst])/alpha)</p>}}{\deqn{\code{max(abs(G[trn,tst])/alpha)}}{max(abs(G[trn,tst])/alpha)}}

to a minimum equal to zero. If \code{alpha=0} the grid is generated starting from a maximum equal to 5}

\item{nlambda}{(integer) Number of lambdas generated when \code{lambda=NULL}}

\item{lambda.min}{(numeric) Minimum value of lambda in the generated grid when \code{lambda=NULL}}

\item{nfolds}{(integer/character) Either 2,3,5,10 or 'n' indicating the number of non-overlaping folds in which the data is split into to do cross-validation. When \code{nfolds='n'} leave-one-out CV is performed}

\item{seed}{(numeric vector) Seed to fix randomization when creating folds for cross-validation. If it is a vector, a number equal to its length of CV repetitions are performed}

\item{nCV}{(integer) Number of CV repetitions to be performed. Default is \code{nCV=1}}

\item{common.lambda}{\code{TRUE} or \code{FALSE} to whether computing the coefficients for a grid of lambdas common to all individuals in testing set
 or for a grid of lambdas specific to each individual in testing set. Default is \code{common.lambda=TRUE}}

\item{mc.cores}{(integer) Number of cores used. When \code{mc.cores} > 1, the analysis is run in parallel for each testing set individual. Default is \code{mc.cores=1}}

\item{tol}{(numeric) Maximum error between two consecutive solutions of the CD algorithm to declare convergence}

\item{maxiter}{(integer) Maximum number of iterations to run the CD algorithm at each lambda step before convergence is reached}

\item{save.at}{(character) Path where files (regression coefficients and output object) are to be saved (this may include a prefix added to the files). Default \code{save.at=NULL} will no save any results and they are returned in the output object}

\item{precision.format}{(character) Either 'single' or 'double' for numeric precision and memory occupancy (4 or 8 bytes, respectively) of the regression coefficients. This is only used when \code{save.at} is not \code{NULL}}

\item{method}{(character) Either 'REML' (Restricted Maximum Likelihood) or 'ML' (Maximum Likelihood) to calculate variance components as per the function 'fitBLUP'}

\item{name}{(character) Name given to the output for tagging purposes. Default \code{name=NULL} will give the name of the method used}

\item{verbose}{\code{TRUE} or \code{FALSE} to whether printing each step}
}
\value{
Function 'SSI' returns a list object of the class 'SSI' for which methods \code{coef}, \code{fitted}, \code{plot}, and \code{summary} exist. Functions 'net.plot' and 'path.plot' can be also used. It contains the elements:
\itemize{
  \item \code{b}: (vector) fixed effects solutions (including the intercept).
  \item \code{Xb}: (vector) product of the design matrix 'X' times the fixed effects solutions.
  \item \code{u}: (matrix) genetic values for testing individuals (in rows) associated to each value of lambda (in columns).
  \item \code{varU}, \code{varE}, \code{h2}: variance components solutions.
  \item \code{alpha}: value for the elastic-net weights used.
  \item \code{lambda}: (matrix) sequence of values of lambda used (in columns) for each testing individual (in rows).
  \item \code{nsup}: (matrix) number of non-zero predictors at each solution given by lambda for each testing individual (in rows).
  \item \code{file_beta}: path where regression coefficients are saved.
}

Function 'SSI.CV' returns a list object of length \code{nCV} of the class 'SSI.CV' for which methods \code{plot} and Optimal cross-validated penalization values can be obtained using the\code{summary} method. Method \code{plot} is also available.

}
\description{
Computes the entire Elastic-Net solution for the regression coefficients of a Selection Index for a grid of values of the penalization parameter.

An optimal penalization can be chosen using cross-validation (CV) within a specific training set.
}
\details{
The basic linear mixed model that relates phenotypes with genetic values is of the form

\ifelse{html}{\out{<p style='text-align:center'><b>y</b> = <b>X b</b> + <b>Z g</b> + <b>e</b></p>}}{\deqn{\textbf{y}=\textbf{X}\textbf{b}+\textbf{Z}\textbf{g}+\textbf{e}}{y = X b + Z g + e}}

where
\ifelse{html}{\out{<b>y</b>}}{\eqn{\textbf{y}}{y}} is a vector with the response,
\ifelse{html}{\out{<b>b</b>}}{\eqn{\textbf{b}}{b}} is the vector of fixed effects,
\ifelse{html}{\out{<b>g</b>}}{\eqn{\textbf{g}}{g}} is the vector of the genetic values of the genotypes,
\ifelse{html}{\out{<b>e</b>}}{\eqn{\textbf{e}}{e}} is the vector of environmental residuals, and
\ifelse{html}{\out{<b>X</b>}}{\eqn{\textbf{X}}{X}} and \ifelse{html}{\out{<b>Z</b>}}{\eqn{\textbf{Z}}{Z}} are design matrices conecting the fixed and genetic effects with replicates. Genetic values are assumed to follow a Normal distribution as
\ifelse{html}{\out{<b>g</b> ~ N(<b>0</b>,&sigma;<sup>2</sup><sub>u</sub><b>K</b>)}}{\eqn{\textbf{g}\sim N(\textbf{0},\sigma^2_u\textbf{K})}{g ~ N(0,sigma^2_u K)}}, and environmental terms are assumed
\ifelse{html}{\out{<b>e</b> ~ N(<b>0</b>,&sigma;<sup>2</sup><sub>e</sub><b>I</b>)}}{\eqn{\textbf{e}\sim N(\textbf{0},\sigma^2_e\textbf{I})}{e ~ N(0,sigma^2_e I)}}.

The resulting vector of genetic values
\ifelse{html}{\out{<b>u</b> = <b>Z g</b>}}{\eqn{\textbf{u}=\textbf{Z}\textbf{g}}{u = Z g}} will therefore follow 
\ifelse{html}{\out{<b>u</b> ~ N(<b>0</b>,&sigma;<sup>2</sup><sub>u</sub><b>G</b>)}}{\eqn{\textbf{u}\sim N(\textbf{0},\sigma^2_u\textbf{G})}{u ~ N(0,sigma^2_u G)}} where
\ifelse{html}{\out{<b>G</b> = <b>Z K Z'</b>}}{\eqn{\textbf{G}=\textbf{Z}\textbf{K}\textbf{Z}'}{G = Z K Z'}}.
In the un-replicated case, \ifelse{html}{\out{<b>Z</b> = <b>I</b>}}{\eqn{\textbf{Z}=\textbf{I}}{Z = I}} is an identity matrix, and hence 
\ifelse{html}{\out{<b>u</b> = <b>g</b>}}{\eqn{\textbf{u}=\textbf{g}}{u = g}} and
\ifelse{html}{\out{<b>G</b> = <b>K</b>}}{\eqn{\textbf{G}=\textbf{K}}{G = K}}.

The values \ifelse{html}{\out{<b>u</b><sub>tst</sub> = (u<sub>i</sub>)}}{\eqn{\textbf{u}_{tst}=(u_i)}{u_tst = (u_i)}},
\ifelse{html}{\out{i = 1,2,...,n<sub>tst</sub>}}{\eqn{i=1,2,...,n_{tst}}{i = 1,2,...,n_tst}}, for a testing set are estimated individual-wise using (as predictors) all available observations in a training set as

\ifelse{html}{\out{<p style='text-align:center'>u<sub>i</sub> = <b>&beta;</b>'<sub>i</sub> (<b>y</b><sub>trn</sub> - <b>X</b><sub>trn</sub><b>b</b>)</p>}}{\deqn{u_i=\boldsymbol{\beta}'_i(\textbf{y}_{trn}-\textbf{X}_{trn}\textbf{b})}{u_i = beta'_i (y_trn - X_trn*b)}}

where \ifelse{html}{\out{<b>&beta;</b><sub>i</sub>}}{\eqn{\boldsymbol{\beta}_i}{beta_i}}
is a vector of weights that are found separately for each individual in the testing set, by minimizing the penalized mean squared error function

\ifelse{html}{\out{<p style='text-align:center'>-&sigma;<sup>2</sup><sub>u</sub><b>G</b>'<sub>trn,tst(i)</sub> <b>&beta;</b><sub>i</sub> + 1/2 <b>&beta;</b>'<sub>i</sub>(&sigma;<sup>2</sup><sub>u</sub><b>G</b><sub>trn</sub> + &sigma;<sup>2</sup><sub>e</sub><b>I</b>)<b>&beta;</b><sub>i</sub> + &lambda; J(<b>&beta;</b><sub>i</sub>)</p>}}{\deqn{-\sigma^2_u\textbf{G}_{trn,tst(i)}' \boldsymbol{\beta}_i + 1/2 \boldsymbol{\beta}_i'(\sigma^2_u\textbf{G}_{trn} + \sigma^2_e\textbf{I}) \boldsymbol{\beta}_i + \lambda J(\boldsymbol{\beta}_i)}{-sigma^2_u G'[trn,tst(i)] beta_i + 1/2 beta_i'(sigma^2_u G[trn] + sigma^2_e I)beta_i + lambda J(beta_i)}}

where
\ifelse{html}{\out{<b>G</b><sub>trn,tst(i)</sub>}}{\eqn{\textbf{G}_{trn,tst(i)}}{G[trn,tst(i)}}
is the \ifelse{html}{\out{i<sup>th</sup>}}{\eqn{i^{th}}{i^th}} column of the
sub-matrix of \ifelse{html}{\out{<b>G</b>}}{\eqn{\textbf{G}}{G}} whose rows correspond to the training set and columns to the testing set;
\ifelse{html}{\out{<b>G</b><sub>trn</sub>}}{\eqn{\textbf{G}_{trn}}{G[trn]}}
is the sub-matrix corresponding to the training set; \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}} is the penalization parameter; and
\ifelse{html}{\out{J(<b>&beta;</b><sub>i</sub>)}}{\eqn{J(\boldsymbol{\beta}_i)}{J(beta_i)}}
is a penalty function given by

\ifelse{html}{\out{<p style='text-align:center'>1/2(1-&alpha;)||<b>&beta;</b><sub>i</sub>||<sub>2</sub><sup>2</sup> + &alpha;||<b>&beta;</b><sub>i</sub>||<sub>1</sub></p>}}{\deqn{1/2(1-\alpha)||\boldsymbol{\beta}_i||_2^2 + \alpha||\boldsymbol{\beta}_i||_1}{1/2(1-alpha)||beta_i||_2^2 + alpha||beta_i||_1}}

where \ifelse{html}{\out{0 &le; &alpha; &le; 1}}{\eqn{0\leq\alpha\leq 1}{0 <= alpha <= 1}}, and
\ifelse{html}{\out{||<b>&beta;</b><sub>i</sub>||<sub>1</sub> = &sum;<sub>j=1</sub>|&beta;<sub>ij</sub>|}}{\eqn{||\boldsymbol{\beta}_i||_1=\sum_{j=1}^{n_{trn}}|\beta_{ij}|}{||beta_i||_1 = sum(|beta_ij|)}} and
\ifelse{html}{\out{||<b>&beta;</b><sub>i</sub>||<sub>2</sub><sup>2</sup> = &sum;<sub>j=1</sub>&beta;<sub>ij</sub><sup>2</sup>}}{\eqn{||\boldsymbol{\beta}_i||_2^2=\sum_{j=1}^{n_{trn}}\beta_{ij}^2}{||beta_i||_2^2 = sum(beta_ij^2)}} are the L1 and (squared) L2-norms, respectively.

Function 'SSI' calculates each individual solution using the function 'solveEN' (via the Coordinate Descent algorithm, see \code{help(solveEN)}) by setting the argument \code{Sigma} equal to
\ifelse{html}{\out{&sigma;<sup>2</sup><sub>u</sub><b>G</b><sub>trn</sub> + &sigma;<sup>2</sup><sub>e</sub><b>I</b>}}{\eqn{\sigma^2_u\textbf{G}_{trn}+\sigma^2_e\textbf{I}}{sigma^2_u G[trn] + sigma^2_e I}}
and \code{Gamma} equal to
\ifelse{html}{\out{&sigma;<sup>2</sup><sub>u</sub><b>G</b><sub>trn,tst(i)</sub>}}{\eqn{\sigma^2_u\textbf{G}_{trn,tst(i)}}{sigma^2_u G[trn,tst(i)]}}.

Function 'SSI.CV' performs cross-validation within the training data specified in argument \code{trn}. Training data is divided into \eqn{k}{k} folds and the SSI is sequentially calculated for (all individuals in) one fold (as testing set) using information from the remaining folds (as training set).
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  index = which(Y$trial \%in\% 1:8)     # Use only a subset of data
  Y = Y[index,]
  M = scale(M[index,])/sqrt(ncol(M))  # Subset and scale markers
  G = tcrossprod(M)                   # Genomic relationship matrix
  y = as.vector(scale(Y[,"E1"]))      # Scale response variable
  
  # Sets (testing=0, training=1)
  trn_tst = ifelse(Y$trial \%in\% 1:2, 0, 1)
  
  # Sparse selection index
  # Fixed effects b and variances varU and varE are 
  # computed internally from training data
  fm1 = SSI(y,K=G,trn_tst=trn_tst)
  summary(fm1)$optCOR
  varU = fm1$varU
  varE = fm1$varE
  b = fm1$b
  \donttest{
  #---------------------------------------------------
  # Predicting a testing set using a value of lambda
  # obtained from cross-validation in a traning set
  #---------------------------------------------------
  # Run a cross validation in training set
  fm2 = SSI.CV(y,K=G,varU=varU,varE=varE,b=b,trn_tst=trn_tst,nfolds=5,name="1 5CV")
  lambda = summary(fm2)$optCOR["lambda"]

  # Fit the index with the obtained lambda
  fm3 = SSI(y,K=G,varU=varU,varE=varE,b=b,trn_tst=trn_tst,lambda=lambda)
  summary(fm3)$accuracy        # Testing set accuracy

  # Compare the accuracy with that of the non-sparse index (G-BLUP)
  summary(fm1)$accuracy[fm1$nlambda,1] # we take the last one
  
  # Obtain an 'optimal' lambda by repeating the CV several times
  fm22 = SSI.CV(y,K=G,varU=varU,varE=varE,b=b,trn_tst=trn_tst,nCV=5,name="5 5CV")
  plot(fm22,fm2)
  
  #---------------------------------------------------
  # Multi-trait SSI
  #---------------------------------------------------
  y = scale(Y[,4:6])    # Response variable
  
  # Sets (testing=0, training=1)
  trn_tst = matrix(NA,ncol=ncol(y),nrow=nrow(y))
  trn_tst[,1] = ifelse(Y$trial \%in\% 1:2, 0, 1)
  trn_tst[,2] = ifelse(Y$trial \%in\% 1, 0, 1)
  trn_tst[,3] = ifelse(Y$trial \%in\% 2, 0, 1)
  
  fm = SSI(y, K=G, trn_tst=trn_tst, mc.cores=1)
  }
}
\references{
Efron B, Hastie T, Johnstone I, Tibshirani R (2004). Least angle regression. \emph{The Annals of Statistics}, \bold{32}(2), 407–499.

Friedman J, Hastie T, Höfling H, Tibshirani R (2007). Pathwise coordinate optimization. \emph{The Annals of Applied Statistics}, \bold{1}(2), 302–332.

Hoerl AE, Kennard RW (1970). Ridge regression: biased estimation for nonorthogonal problems. \emph{Technometrics}, \bold{12}(1), 55–67.

Lush JL (1947). Family merit an individual merit as bases for selection. Part I. \emph{The American Naturalist}, \bold{81}(799), 241–261.

Tibshirani R (1996). Regression shrinkage and selection via the LASSO. \emph{Journal of the Royal Statistical Society B}, \bold{58}(1), 267–288.

VanRaden PM (2008). Efficient methods to compute genomic predictions. \emph{Journal of Dairy Science}, \bold{91}(11), 4414–4423.

Zou H, Hastie T (2005). Regularization and variable selection via the elastic net. \emph{Journal of the Royal Statistical Society B}, \bold{67}(2), 301–320
}

\keyword{SSI}
