
\name{Methods_LASSO}
\alias{coef.LASSO}
\alias{fitted.LASSO}
\title{LASSO methods}
\usage{
\method{coef}{LASSO}(object, ..., i=NULL)

\method{fitted}{LASSO}(object, ...)
}
\arguments{
\item{object}{An object of the class 'LASSO' returned either by the function 'LARS' or 'solveEN'}

\item{...}{Other arguments: \code{X} (numeric matrix) scores for as many predictors there are in \code{ncol(object$beta)}  (in columns) for a desired number \code{n} of observations (in rows)}

\item{i}{(integer vector) Index columns of matrix 'Gamma' to be considered}
}
\value{
Method \code{coef} returns a matrix that contains, for each value of lambda (in columns), the predicted values corresponding to each row of the matrix \code{X}.

Method \code{fitted} returns fitted values
\ifelse{html}{\out{X&beta;}}{\eqn{\textbf{X}\boldsymbol{\beta}}{X beta}}
}
\description{
Predicted values for a provided matrix of predictors \code{X} 
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  y = as.vector(Y[,"E1"])  # Response variable
  X = scale(X_E1)          # Predictors
  
  # Training and testing sets
  tst = which(Y$trial \%in\% 1:10)
  trn = seq_along(y)[-tst]

  # Calculate covariances in training set
  XtX = var(X[trn,])
  Xty = cov(X[trn,],y[trn])
  
  # Run the penalized regression
  fm = solveEN(XtX,Xty,alpha=0.5)   
  
  # Regression coefficients
  B = coef(fm)
  
  # Predicted values
  yHat1 = fitted(fm, X=X[trn,])  # training data
  yHat2 = fitted(fm, X=X[tst,])  # testing data
  
  # Penalization vs correlation
  plot(-log(fm$lambda[-1]),cor(y[trn],yHat1[,-1]), main="training")
  plot(-log(fm$lambda[-1]),cor(y[tst],yHat2[,-1]), main="testing")
}
\author{
Marco Lopez-Cruz (\email{maraloc@gmail.com}) and Gustavo de los Campos
}
