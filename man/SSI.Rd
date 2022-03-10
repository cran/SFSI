
\name{SSI}
\alias{SSI}
\alias{SSI.CV}
\title{Sparse Selection Index}
\usage{
SSI(y, X = NULL, b = NULL, Z = NULL, K, D = NULL, 
    theta = NULL, h2 = NULL, trn = seq_along(y), tst = seq_along(y), 
    subset = NULL, alpha = 1, lambda = NULL, nlambda = 100,
    lambda.min = .Machine$double.eps^0.5, common.lambda = TRUE,
    tol = 1E-4, maxiter = 500, method = c("REML","ML"),
    save.at = NULL, name = NULL, mc.cores = 1, verbose = TRUE)
    
SSI.CV(y, X = NULL, b = NULL, Z = NULL, K, D = NULL, 
      theta = NULL, h2 = NULL, trn = seq_along(y), alpha = 1, 
      lambda = NULL, nlambda = 100, lambda.min = .Machine$double.eps^0.5,
      nCV = 1, nfolds = 5, seed = NULL, common.lambda = TRUE,
      tol = 1E-4, maxiter = 500, method = c("REML","ML"), 
      name = NULL, mc.cores = 1, verbose = TRUE)
}
\arguments{
\item{y}{(numeric vector) Response variable}

\item{X}{(numeric matrix) Design matrix for fixed effects. When \code{X=NULL}, a vector of ones is constructed only for the intercept (default)}

\item{b}{(numeric vector) Fixed effects. When \code{b=NULL}, only the intercept is estimated from training data using generalized least squares (default)}

\item{Z}{(numeric matrix) Design matrix for random effects. When \code{Z=NULL} an identity matrix is considered (default) thus \ifelse{html}{\out{<b>G</b> = <b>K</b>}}{\eqn{\textbf{G = K}}{G = K}}; otherwise \ifelse{html}{\out{<b>G</b> = <b>Z K Z'</b>}}{\eqn{\textbf{G} = \textbf{Z K Z'}}{G = Z K Z'}} is used}

\item{K}{(numeric matrix) Kinship relationships. This can be of the "float32" type as per the 'float' R-package, or a (character) name of a binary file where the matrix is stored}

\item{D}{(numeric matrix) Relationships among residuals (usually an identity matrix). When \code{D=NULL} an identity matrix is considered (default)}

\item{theta}{(numeric) Residual/genetic variances ratio. When \code{theta=NULL} (default), it is calculated from training data using the function 'fitBLUP' (see \code{help(fitBLUP)})}

\item{h2}{(numeric) Heritability of the response variable. When \code{h2=NULL} (default), it is calculated from training data using the function 'fitBLUP' (see \code{help(fitBLUP)}). It is used to set the residual/genetic variances ratio \code{theta}}

\item{trn}{(integer vector) Which elements from vector \code{y} are in training set. Default \code{trn=seq_along(y)} will consider all individuals as training}

\item{tst}{(integer vector) Which elements from vector \code{y} are in testing set (prediction set). Default \code{tst=seq_along(y)} will consider all individuals as testing}

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

\item{mc.cores}{(integer) Number of cores used. The analysis is run in parallel when \code{mc.cores} is greater than 1. Default is \code{mc.cores=1}}

\item{tol}{(numeric) Maximum error between two consecutive solutions of the CD algorithm to declare convergence}

\item{maxiter}{(integer) Maximum number of iterations to run the CD algorithm at each lambda step before convergence is reached}

\item{save.at}{(character) Prefix name that will be added to the output files name to be saved, this may include a path. Regression coefficients
are saved as a binary file (single-precision: 32 bits, 7 significant digits). Default \code{save.at=NULL} will no save any output}

\item{method}{(character) Either 'REML' (Restricted Maximum Likelihood) or 'ML' (Maximum Likelihood) to calculate variance components as per the function 'fitBLUP'}

\item{name}{(character) Name given to the output for tagging purposes. Default \code{name=NULL} will give the name of the method used}

\item{verbose}{\code{TRUE} or \code{FALSE} to whether printing each step}
}
\value{
Function 'SSI' returns a list object of the class 'SSI' for which methods \code{coef}, \code{fitted}, \code{plot}, and \code{summary} exist. Functions 'net.plot' and 'path.plot' can be also used. It contains the elements:
\itemize{
  \item \code{b}: (vector) fixed effects solutions (including the intercept).
  \item \code{Xb}: (vector) product of the design matrix 'X' times the fixed effects solutions.
  \item \code{varU}, \code{varE}, \code{h2}: variance components solutions.
  \item \code{alpha}: value for the elastic-net weights used.
  \item \code{lambda}: (matrix) sequence of values of lambda used (in columns) for each testing individual (in rows).
  \item \code{df}: (matrix) degrees of freedom (number of non-zero predictors) at each solution given by lambda for each testing individual (in rows).
  \item \code{file_beta}: path where regression coefficients are saved.
}

Function 'SSI.CV' returns a list object of length \code{nCV} of the class 'SSI.CV' for which methods \code{plot} and \code{summary} exist. Each element is also a list containing the elements:
\itemize{
  \item \code{b}: (vector) solutions for the fixed effects (including the intercept) for each fold.
  \item \code{varU}, \code{varE}, \code{h2}: variance components estimated within each fold.
  \item \code{folds}: (matrix) assignation of training individuals to folds used for the cross-validation.
  \item \code{accuracy}: (matrix) correlation between observed and predicted values (in testing set) within each fold (in rows).
  \item \code{MSE}: (matrix) mean squared error of prediction (in testing set) within each fold (in rows).
  \item \code{lambda}: (matrix) with the sequence of values of lambda used (averaged across individuals) within each fold (in rows).
  \item \code{df}: (matrix) with the degrees of freedom (averaged across individuals) within each fold (in rows).
}
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
\ifelse{html}{\out{<b>e</b> ~ N(<b>0</b>,&sigma;<sup>2</sup><sub>e</sub><b>D</b>)}}{\eqn{\textbf{e}\sim N(\textbf{0},\sigma^2_e\textbf{D})}{e ~ N(0,sigma^2_e D)}}, usually \ifelse{html}{\out{<b>D</b>=<b>I</b>}}{\eqn{\textbf{D}=\textbf{I}}{D = I}}.

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

\ifelse{html}{\out{<p style='text-align:center'>-<b>G</b>'<sub>trn,tst(i)</sub> <b>&beta;</b><sub>i</sub> + 1/2 <b>&beta;</b>'<sub>i</sub>(<b>G</b><sub>trn</sub> + &theta;<b>D</b>)<b>&beta;</b><sub>i</sub> + &lambda; J(<b>&beta;</b><sub>i</sub>)</p>}}{\deqn{-\textbf{G}_{trn,tst(i)}' \boldsymbol{\beta}_i + 1/2 \boldsymbol{\beta}_i'(\textbf{G}_{trn} + \theta\textbf{D}) \boldsymbol{\beta}_i + \lambda J(\boldsymbol{\beta}_i)}{-G'[trn,tst(i)] beta_i + 1/2 beta_i'(G[trn] + theta D)beta_i + lambda J(beta_i)}}

where
\ifelse{html}{\out{<b>G</b><sub>trn,tst(i)</sub>}}{\eqn{\textbf{G}_{trn,tst(i)}}{G[trn,tst(i)}}
is the \ifelse{html}{\out{i<sup>th</sup>}}{\eqn{i^{th}}{i^th}} column of the
sub-matrix of \ifelse{html}{\out{<b>G</b>}}{\eqn{\textbf{G}}{G}} whose rows correspond to the training set and columns to the testing set;
\ifelse{html}{\out{<b>G</b><sub>trn</sub>}}{\eqn{\textbf{G}_{trn}}{G[trn]}}
is the sub-matrix corresponding to the training set; \ifelse{html}{\out{&theta; = &sigma;<sup>2</sup><sub>e</sub>/&sigma;<sup>2</sup><sub>u</sub>}}{\eqn{\theta=\sigma^2_e/\sigma^2_u}{theta = sigma^e/sigma^2_u}} is the residual to genetic variance ratio that can be expressed in terms of the heritability, \ifelse{html}{\out{h<sup>2</sup> = &sigma;<sup>2</sup><sub>u</sub>/(&sigma;<sup>2</sup><sub>u</sub> + &sigma;<sup>2</sup><sub>e</sub>)}}{\eqn{h^2=\sigma^2_u/(\sigma^2_u+\sigma^2_e)}{h^2 = sigma^2_u/(sigma^2_u + sigma^2_e)}}, as \ifelse{html}{\out{&theta; = (1 - h<sup>2</sup>)/h<sup>2</sup>}}{\eqn{\theta=(1-h^2)/h^2}{theta = (1 - h^2)/h^2}}; \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}} is the penalization parameter; and
\ifelse{html}{\out{J(<b>&beta;</b><sub>i</sub>)}}{\eqn{J(\boldsymbol{\beta}_i)}{J(beta_i)}}
is a penalty function given by

\ifelse{html}{\out{<p style='text-align:center'>1/2(1-&alpha;)||<b>&beta;</b><sub>i</sub>||<sub>2</sub><sup>2</sup> + &alpha;||<b>&beta;</b><sub>i</sub>||<sub>1</sub></p>}}{\deqn{1/2(1-\alpha)||\boldsymbol{\beta}_i||_2^2 + \alpha||\boldsymbol{\beta}_i||_1}{1/2(1-alpha)||beta_i||_2^2 + alpha||beta_i||_1}}

where \ifelse{html}{\out{0 &le; &alpha; &le; 1}}{\eqn{0\leq\alpha\leq 1}{0 <= alpha <= 1}}, and
\ifelse{html}{\out{||<b>&beta;</b><sub>i</sub>||<sub>1</sub> = &sum;<sub>j=1</sub>|&beta;<sub>ij</sub>|}}{\eqn{||\boldsymbol{\beta}_i||_1=\sum_{j=1}^{n_{trn}}|\beta_{ij}|}{||beta_i||_1 = sum(|beta_ij|)}} and
\ifelse{html}{\out{||<b>&beta;</b><sub>i</sub>||<sub>2</sub><sup>2</sup> = &sum;<sub>j=1</sub>&beta;<sub>ij</sub><sup>2</sup>}}{\eqn{||\boldsymbol{\beta}_i||_2^2=\sum_{j=1}^{n_{trn}}\beta_{ij}^2}{||beta_i||_2^2 = sum(beta_ij^2)}} are the L1 and (squared) L2-norms, respectively.

Function 'SSI' calculates each individual solution using the function 'solveEN' (via the Coordinate Descent algorithm, see \code{help(solveEN)}) by setting the argument \code{Sigma} equal to
\ifelse{html}{\out{<b>G</b><sub>trn</sub> + &theta;<b>D</b>}}{\eqn{\textbf{G}_{trn}+\theta\textbf{D}}{G[trn] + theta D}}
and \code{Gamma} equal to
\ifelse{html}{\out{<b>G</b><sub>trn,tst(i)</sub>}}{\eqn{\textbf{G}_{trn,tst(i)}}{G[trn,tst(i)]}}.

Function 'SSI.CV' performs cross-validation within the training data specified in argument \code{trn}. Training data is divided into \eqn{k}{k} folds and the SSI is sequentially calculated for (all individuals in) one fold (as testing set) using information from the remaining folds (as training set).
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  index = which(Y$CV == 1)              # Use only a subset of data
  M = scale(M[index,])/sqrt(ncol(M))    # Subset and scale markers
  G = tcrossprod(M)                     # Genomic relationship matrix
  y = as.vector(scale(Y[index,"E1"]))   # Subset ans scale response variable
  
  # Training and testing sets
  tst = seq(1,length(y),by=3)
  trn = (seq_along(y))[-tst]

  # Calculate variance components ratio using training data
  yNA <- y
  yNA[tst] <- NA
  fm0 = fitBLUP(yNA,K=G)
  theta = fm0$varE/fm0$varU
  h2 = fm0$varU/(fm0$varU + fm0$varE)
  
  # Sparse selection index
  fm1 = SSI(y,K=G,theta=theta,trn=trn,tst=tst)
  summary(fm1)$optCOR
  
  \donttest{
  if(requireNamespace("float")){
  # Using a 'float' type variable for K
  G2 = float::fl(G)
  fm2 = SSI(y,K=G2,theta=theta,trn=trn,tst=tst)
  summary(fm2)$optCOR  # compare with above results
  }
  }
  
  #---------------------------------------------------
  # Predicting a testing set using a value of lambda
  # obtained from cross-validation in a traning set
  #---------------------------------------------------
  # Run a cross validation in training set
  fm2 = SSI.CV(y,K=G,theta=theta,trn=trn,nfolds=5,name="1 5CV")
  lambda = summary(fm2)$optCOR["lambda"]
  
  # Fit the index with the obtained lambda
  fm3 = SSI(y,K=G,theta=theta,trn=trn,tst=tst,lambda=lambda)
  summary(fm3)$accuracy        # Testing set accuracy

  # Compare the accuracy with that of the non-sparse index (G-BLUP)
  cor(fm0$u[tst],y[tst])
  
  \donttest{
  # Obtain an 'optimal' lambda by repeating the CV several times
  fm22 = SSI.CV(y,K=G,theta=theta,trn=trn,nCV=5,name="5 5CV")
  plot(fm22,fm2)
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
\author{
Marco Lopez-Cruz (\email{maraloc@gmail.com}) and Gustavo de los Campos
}
\keyword{SSI}
