
\name{solveEN}
\alias{solveEN}
\title{Coordinate Descent algorithm to solve Elastic-Net-type problems}
\usage{
solveEN(P, v, alpha = 1, lambda = NULL, nLambda = 100, 
   minLambda = .Machine$double.eps^0.5, scale = TRUE, 
   tol = 1e-05, maxIter = 1000, verbose = FALSE)
}
\arguments{
\item{P}{(numeric matrix) Variance-covariance matrix of predictors. It can be of the "float32" type as per the 'float' R-package}

\item{v}{(numeric vector) Covariance between response variable and predictors}

\item{lambda}{(numeric vector) Penalization parameter sequence. Default is \code{lambda=NULL}, in this case a decreasing grid of
\code{'nLambda'} lambdas will be generated starting from a maximum equal to
\ifelse{html}{\out{<center><font face="Courier">max(abs(v)/alpha)</font></center>}}{\deqn{\code{max(abs(v)/alpha)}}{max(abs(v)/alpha)}}
to a minimum equal to zero. If \code{alpha=0} the grid is generated starting from a maximum equal to 5}

\item{nLambda}{(integer) Number of lambdas generated when \code{lambda=NULL}}

\item{minLambda}{(numeric) Minimum value of lambda that are generated when \code{lambda=NULL}}

\item{alpha}{(numeric) Value between 0 and 1 for the weights given to the L1 and L2-penalties}

\item{scale}{\code{TRUE} or \code{FALSE} to scale matrix \code{P} for variables with unit variance and scale \code{v} by the standard deviation of the corresponding predictor
taken from the diagonal of \code{P}}

\item{tol}{(numeric) Maximum error between two consecutive solutions of the CD algorithm to declare convergence}

\item{maxIter}{(integer) Maximum number of iterations to run the CD algorithm at each lambda step before convergence is reached}

\item{verbose}{\code{TRUE} or \code{FALSE} to whether printing each CD step}
}
\value{
Returns a list object containing the elements:
\itemize{
  \item \code{lambda}: (vector) all the sequence of values of the penalty.
  \item \code{beta}: (matrix) regression coefficients for each predictor (in columns) associated to each value of the penalization parameter lambda (in rows).
  \item \code{df}: (vector) degrees of freedom, number of non-zero predictors associated to each value of lambda.
}
The returned object is of the class 'LASSO' for which methods \code{fitted} exist. Function 'plotPath' can be also used
}
\description{
Computes the entire Elastic-Net solution for the regression coefficients for all
values of the penalization parameter, via the Coordinate Descent (CD) algorithm (Friedman, 2007). It uses as inputs a variance matrix among predictors and a covariance vector between response and predictors
}
\details{
Finds solutions for the regression coefficients in a linear model
\ifelse{html}{\out{<center>y<sub>i</sub> = <b>x</b>'<sub>i</sub> <b>&beta;</b> + e<sub>i</sub></center>}}{\deqn{y_i=\textbf{x}_i'\boldsymbol{\beta}+e_i}{y_i = x'_i beta + e_i}}

where
\ifelse{html}{\out{y<sub>i</sub>}}{\eqn{y_i}{y_i}} is the response for the \ifelse{html}{\out{i<sup>th</sup>}}{\eqn{i^{th}}{i^th}} observation,
\ifelse{html}{\out{<b>x</b><sub>i</sub>=(x<sub>i1</sub>,...,x<sub>ip</sub>)'}}{\eqn{\textbf{x}_i=(x_{i1},...,x_{ip})'}{x_i=(x_i1,...,x_ip)'}}
is a vector of \eqn{p}{p} predictors assumed to have unit variance,
\ifelse{html}{\out{<b>&beta;</b>=(&beta;<sub>1</sub>,...,&beta;<sub>p</sub>)'}}{\eqn{\boldsymbol{\beta}=(\beta_1,...,\beta_p)'}{beta=(beta[1],...,beta[p])'}}
is a vector of regression coefficients, and
\ifelse{html}{\out{e<sub>i</sub>}}{\eqn{e_i}{e_i}}
is a residual.

The regression coefficients
\ifelse{html}{\out{<b>&beta;</b>}}{\eqn{\boldsymbol{\beta}}{beta}}
are estimated as function of the variance matrix among
predictors (\ifelse{html}{\out{<b>P</b>}}{\eqn{\textbf{P}}{P}}) and the covariance vector between response and predictors (\ifelse{html}{\out{<b>v</b>}}{\eqn{\textbf{v}}{v}}) by minimizing the penalized mean squared error function

\ifelse{html}{\out{<center>-<b>v</b>' <b>&beta;</b> + 1/2 <b>&beta;</b>' <b>P</b> <b>&beta;</b> + &lambda; J(<b>&beta;</b>)</center>}}{\deqn{-\textbf{v}' \boldsymbol{\beta} + 1/2 \boldsymbol{\beta}'\textbf{P} \boldsymbol{\beta} + \lambda J(\boldsymbol{\beta})}{-v' beta + 1/2 beta' P beta + lambda J(beta)}}

where \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}}
is the penalization parameter and
\ifelse{html}{\out{J(<b>&beta;</b>)}}{\eqn{J(\boldsymbol{\beta})}{J(beta)}}
is a penalty function given by

\ifelse{html}{\out{<center>1/2(1-&alpha;)||<b>&beta;</b>||<sub>2</sub><sup>2</sup> + &alpha;||<b>&beta;</b>||<sub>1</sub></center>}}{\deqn{1/2(1-\alpha)||\boldsymbol{\beta}||_2^2 + \alpha||\boldsymbol{\beta}||_1}{1/2(1-alpha)||beta||_2^2 + alpha||beta||_1}}

where \ifelse{html}{\out{0 &le; &alpha; &le; 1}}{\eqn{0\leq\alpha\leq 1}{0 <= alpha <= 1}}, and
\ifelse{html}{\out{||<b>&beta;</b>||<sub>1</sub> = &sum;|&beta;<sub>j</sub>|}}{\eqn{||\boldsymbol{\beta}||_1=\sum|\beta_j|}{||beta||_1 = sum(|beta_j|)}} and
\ifelse{html}{\out{||<b>&beta;</b>||<sub>2</sub><sup>2</sup> = &sum;&beta;<sub>j</sub><sup>2</sup>}}{\eqn{||\boldsymbol{\beta}||_2^2=\sum\beta_j^2}{||beta||_2^2 = sum(beta_j^2)}} are the L1 and (squared) L2-norms, respectively.

The "partial residual" excluding the contribution of the predictor
\ifelse{html}{\out{x<sub>ij</sub>}}{\eqn{x_{ij}}{x_ij}} is
\ifelse{html}{\out{<center>e<sub>i</sub><sup>(j)</sup> = y<sub>i</sub> - <b>x</b>'<sub>i</sub> <b>&beta;</b> + x<sub>ij</sub>&beta;<sub>j</sub></center>}}{\deqn{e_i^{(j)}=y_i-\textbf{x}_i'\boldsymbol{\beta}+x_{ij}\beta_j}{r_i = y_i - x'_i beta + x_ij beta_j}}

then the ordinary least-squares (OLS) coefficient of
\ifelse{html}{\out{x<sub>ij</sub>}}{\eqn{x_{ij}}{x_ij}}
on this residual is (up-to a constant)

\ifelse{html}{\out{<center>&beta;<sub>j</sub><sup>(ols)</sup> = v<sub>j</sub> - <b>P</b>'<sub>j</sub> <b>&beta;</b> + &beta;<sub>j</sub></center>}}{\deqn{\beta_j^{(ols)}=v_j-\textbf{P}_j'\boldsymbol{\beta}+\beta_j}{beta_j^ols = v_j - P'_j beta}}

where
\ifelse{html}{\out{v<sub>j</sub>}}{\eqn{v_j}{v_j}}
is the \ifelse{html}{\out{j<sup>th</sup>}}{\eqn{j^{th}}{j^th}} element of
\ifelse{html}{\out{<b>v</b>}}{\eqn{\textbf{v}}{v}} and
\ifelse{html}{\out{<b>P</b><sub>j</sub>}}{\eqn{\textbf{P}_j}{P_j}}
is the \ifelse{html}{\out{j<sup>th</sup>}}{\eqn{j^{th}}{j^th}} column of the matrix
\ifelse{html}{\out{<b>P</b>}}{\eqn{\textbf{P}}{P}}.

Coefficients are updated for each \eqn{j=1,...,p}{j=1,...,p} from their current value
\ifelse{html}{\out{&beta;<sub>j</sub>}}{\eqn{\beta_j}{beta_j}}
to a new value
\ifelse{html}{\out{&beta;<sub>j</sub>(&alpha;,&lambda;)}}{\eqn{\beta_j(\alpha,\lambda)}{beta_j(alpha,lambda)}},
given \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}} and
\ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}},
by "soft-thresholding" their OLS estimate until convergence as fully described in Friedman (2007).
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  y = as.vector(Y[,"YLD"])  # Response variable
  X = scale(WL)             # Predictors

  # Training and testing sets
  tst = 1:ceiling(0.3*length(y))
  trn = seq_along(y)[-tst]

  # Calculate covariances in training set
  XtX = var(X[trn,])
  Xty = cov(y[trn],X[trn,])
  
  # Run the penalized regression
  fm1 = solveEN(XtX,Xty,alpha=0.5) 
  
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
   fm2 = solveEN(XtX2,Xty,alpha=0.5)  
   max(abs(fm1$beta-fm2$beta))  # Check for discrepances
  }
  }
}
\references{
Friedman J, Hastie T, Höfling H, Tibshirani R (2007). Pathwise coordinate optimization. \emph{The Annals of Applied Statistics}, \bold{1}(2), 302–332.

Hoerl AE, Kennard RW (1970). Ridge Regression: Biased estimation for nonorthogonal problems. \emph{Technometrics}, \bold{12}(1), 55–67.

Tibshirani R (1996). Regression shrinkage and selection via the LASSO. \emph{Journal of the Royal Statistical Society B}, \bold{58}(1), 267–288.

Zou H, Hastie T (2005). Regularization and variable selection via the elastic net. \emph{Journal of the Royal Statistical Society B}, \bold{67}(2), 301–320.
}
\author{
Marco Lopez-Cruz (\email{maraloc@gmail.com}) and Gustavo de los Campos
}
\keyword{solveEN}
