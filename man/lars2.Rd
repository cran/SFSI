
\name{lars2}
\alias{lars2}
\title{Least Angle Regression to solve LASSO-type problems}
\usage{
lars2(Sigma, Gamma, method = c("LAR", "LAR-LASSO"), maxDF = NULL,
  eps = .Machine$double.eps, scale = TRUE, verbose = FALSE)
}
\arguments{
\item{Sigma}{(numeric matrix) Variance-covariance matrix of predictors. It can be of the "float32" type as per the 'float' R-package}

\item{Gamma}{(numeric vector) Covariance between response variable and predictors}

\item{method}{(character) Either:
\itemize{
 \item \code{'LAR'}: Computes the entire sequence of all coefficients. Values of lambdas are calculated at each step.
 \item \code{'LAR-LASSO'}: Similar to \code{'LAR'} but solutions when a predictor leaves the solution are also returned.
}
Default is \code{method='LAR'}}

\item{maxDF}{(integer) Maximum number of non-zero coefficients in the last LARS solution.
Default \code{maxDF=NULL} will calculate solutions for the entire lambda sequence}

\item{eps}{(numeric) An effective zero. Default is the machine precision}

\item{scale}{\code{TRUE} or \code{FALSE} to scale matrix \code{Sigma} for variables with unit variance and scale \code{Gamma} by the standard deviation of the corresponding predictor
taken from the diagonal of \code{Sigma}}

\item{verbose}{\code{TRUE} or \code{FALSE} to whether printing each LARS step}
}
\value{
Returns a list object with the following elements:
\itemize{
  \item \code{lambda}: (vector) all the sequence of values of the LASSO penalty.
  \item \code{beta}: (matrix) regression coefficients for each predictor (in rows) associated to each value of the penalization parameter lambda (in columns).
  \item \code{df}: (vector) degrees of freedom, number of non-zero predictors associated to each value of lambda.
}
The returned object is of the class 'LASSO' for which methods \code{fitted} exist. Function \code{plotPath} can be also used
}
\description{
Computes the entire LASSO solution for the regression coefficients, starting from zero, to the
least-squares estimates, via the Least Angle Regression (LARS) algorithm (Efron, 2004). It uses as inputs
a variance matrix among predictors and a covariance vector between response and predictors.
}
\details{
Finds solutions for the regression coefficients in a linear model
\ifelse{html}{\out{<center>y<sub>i</sub> = <b>x</b>'<sub>i</sub> <b>&beta;</b> + e<sub>i</sub></center>}}{\deqn{y_i=\textbf{x}_i'\boldsymbol{\beta}+e_i}{y_i = x'_i beta + e_i}}

where
\ifelse{html}{\out{y<sub>i</sub>}}{\eqn{y_i}{y_i}} is the response for the \ifelse{html}{\out{i<sup>th</sup>}}{\eqn{i^{th}}{i^th}} observation,
\ifelse{html}{\out{<b>x</b><sub>i</sub>=(x<sub>i1</sub>,...,x<sub>ip</sub>)'}}{\eqn{\textbf{x}_i=(x_{i1},...,x_{ip})'}{x_i=(x_i1,...,x_ip)'}}
is a vector of \eqn{p}{p} predictors assumed to have unit variance,
\ifelse{html}{\out{<b>&beta;</b>=(&beta;<sub>1</sub>,...,&beta;<sub>p</sub>)'}}{\eqn{\boldsymbol{\beta}=(\beta_1,...,\beta_p)'}{beta=(beta_1,...,beta_p)'}}
is a vector of regression coefficients, and
\ifelse{html}{\out{e<sub>i</sub>}}{\eqn{e_i}{e[i]}}
is a residual.

The regression coefficients
\ifelse{html}{\out{<b>&beta;</b>}}{\eqn{\boldsymbol{\beta}}{beta}}
are estimated as function of the variance matrix among
predictors (\ifelse{html}{\out{<b>&Sigma;</b>}}{\eqn{\boldsymbol{\Sigma}}{Sigma}}) and the covariance vector between response and predictors (\ifelse{html}{\out{<b>&Gamma;</b>}}{\eqn{\boldsymbol{\Gamma}}{Gamma}}) by minimizing the penalized mean squared error function

\ifelse{html}{\out{<center>-<b>&Gamma;</b>' <b>&beta;</b> + 1/2 <b>&beta;</b>'<b>&Sigma;</b><b>&beta;</b> + 1/2 &lambda; ||<b>&beta;</b>||<sub>1</sub></center>}}{\deqn{-\boldsymbol{\Gamma}' \boldsymbol{\beta} + 1/2 \boldsymbol{\beta}'\boldsymbol{\Sigma} \boldsymbol{\beta} + 1/2\lambda||\boldsymbol{\beta}||_1}{-Gamma' beta + 1/2 beta' Sigma beta + 1/2 lambda ||beta||_1}}

where \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}}
is the penalization parameter and
\ifelse{html}{\out{||<b>&beta;</b>||<sub>1</sub> = &sum;<sub>j=1</sub>|&beta;<sub>j</sub>|}}{\eqn{||\boldsymbol{\beta}||_1=\sum_{j=1}^p{|\beta_j|}}{||beta||_1=sum(|beta_j|)}}
is the L1-norm.

The algorithm to find solutions for each \ifelse{html}{\out{&beta;<sub>j</sub>}}{\eqn{\beta_j}{beta_j}} is fully described in Efron (2004) in which the "current correlation" between the predictor
\ifelse{html}{\out{x<sub>ij</sub>}}{\eqn{x_{ij}}{xij}}
and the residual
\ifelse{html}{\out{e<sub>i</sub> = y<sub>i</sub> - <b>x</b>'<sub>i</sub> <b>&beta;</b>}}{\eqn{e_i=y_i-\textbf{x}_i'\boldsymbol{\beta}}{e_i = y_i - x'_i beta}}
is expressed (up-to a constant) as

\ifelse{html}{\out{<center>r<sub>j</sub> = &Gamma;<sub>j</sub> - <b>&Sigma;</b>'<sub>j</sub> <b>&beta;</b></center>}}{\deqn{r_j=\Gamma_j-\boldsymbol{\Sigma}_j'\boldsymbol{\beta}}{r_j = Gamma_j - Sigma'_j beta}}

where
\ifelse{html}{\out{&Gamma;<sub>j</sub>}}{\eqn{\Gamma_j}{Gamma_j}}
is the \ifelse{html}{\out{j<sup>th</sup>}}{\eqn{j^{th}}{j^th}} element of
\ifelse{html}{\out{<b>&Gamma;</b>}}{\eqn{\boldsymbol{\Gamma}}{Gamma}} and
\ifelse{html}{\out{<b>&Sigma;</b><sub>j</sub>}}{\eqn{\boldsymbol{\Sigma}_j}{Sigma_j}}
is the \ifelse{html}{\out{j<sup>th</sup>}}{\eqn{j^{th}}{j^th}} column of the matrix
\ifelse{html}{\out{<b>&Sigma;</b>}}{\eqn{\boldsymbol{\Sigma}}{Sigma}}
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  y = as.vector(Y[,"E1"])  # Response variable
  X = scale(X_E1)           # Predictors

  # Training and testing sets
  tst = 1:ceiling(0.3*length(y))
  trn = seq_along(y)[-tst]

  # Calculate covariances in training set
  XtX = var(X[trn,])
  Xty = cov(y[trn],X[trn,])
  
  # Run the penalized regression
  fm1 = lars2(XtX,Xty,method="LAR-LASSO")  
  
  # Predicted values
  yHat1 = fitted(fm1, X=X[trn,])  # training data
  yHat2 = fitted(fm1, X=X[tst,])  # testing data
  
  # Penalization vs correlation
  plot(-log(fm1$lambda[-1]),cor(y[trn],yHat1[,-1]), main="training")
  plot(-log(fm1$lambda[-1]),cor(y[tst],yHat2[,-1]), main="testing")
  
  \donttest{
  if(requireNamespace("float")){
   # Using a 'float' type variable
   XtX2 = float::fl(XtX)
   fm2 = lars2(XtX2,Xty,method="LAR-LASSO")  
   max(abs(fm1$beta-fm2$beta))      # Check for discrepances in beta
   max(abs(fm1$lambda-fm2$lambda))  # Check for discrepances in lambda
  }
  }
}
\references{
Efron B, Hastie T, Johnstone I, Tibshirani R (2004). Least angle regression. \emph{The Annals of Statistics}, \bold{32}(2), 407–499.

Friedman J, Hastie T, Tibshirani R(2010). Regularization paths for generalized linear models via coordinate descent. \emph{Journal of Statistical Software}, \bold{33}(1), 1–22.

Hastie T, Efron B (2013). lars: least angle regression, Lasso and forward stagewise. \url{https://cran.r-project.org/package=lars}.

Tibshirani R (1996). Regression shrinkage and selection via the LASSO. \emph{Journal of the Royal Statistical Society B}, \bold{58}(1), 267–288.
}
\author{
Marco Lopez-Cruz (\email{maraloc@gmail.com}) and Gustavo de los Campos. Adapted from the 'lars' function in package 'lars' (Hastie & Efron, 2013)
}
\keyword{lars2}
