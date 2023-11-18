
\name{Regression coefficients and fitted values in LASSO-type problems}
\alias{coef.LASSO}
\alias{fitted.LASSO}
\title{LASSO methods}
\usage{
\method{coef}{LASSO}(object, ...)

\method{fitted}{LASSO}(object, ...)
}
\arguments{
\item{object}{An object of the class 'LASSO' returned either by the function 'LARS' or 'solveEN'}

\item{...}{Other arguments: 
 \itemize{
  \item \code{X} (numeric matrix) scores for as many predictors there are in \code{ncol(object$beta)}  (in columns) for a desired number \code{n}  of observations (in rows)
  \item \code{iy} (integer vector) Optional index of columns of the matrix 'Gamma' to be returned in \code{coef} function
    \item \code{ilambda} (integer) Optional to return regression coefficients associated to a specific penalty position
  \item \code{nsup} (numeric) Optional to return regression coefficients associated to a given penalty that yield approximately 'nsup' non-zero coefficients
 }
}
}
\value{
Method \code{coef} returns a matrix that contains the regression coefficients (in rows) associated to each value of lambda (in columns). When the regression was applied to an object \code{Gamma} with more than one column, method \code{coef} returns a list

Method \code{fitted} returns a matrix with fitted values
\ifelse{html}{\out{X&beta;}}{\eqn{\textbf{X}\boldsymbol{\beta}}{X beta}} (in rows)
for each value of lambda (in columns).
}
\description{
Retrieving regression coefficients and predicted values from the 'solveEN' and 'LARS' functions' outputs
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
  dim(coef(fm))
  dim(coef(fm, ilambda=50)) # Coefficients associated to the 50th lambda
  dim(coef(fm, nsup=25))    # Coefficients where around nsup=25 are non-zero
  
  # Predicted values
  yHat1 = fitted(fm, X=X[trn,])  # training data
  yHat2 = fitted(fm, X=X[tst,])  # testing data
  
  # Penalization vs correlation
  plot(-log(fm$lambda[-1]),cor(y[trn],yHat1[,-1]),main="training",type="l")
  plot(-log(fm$lambda[-1]),cor(y[tst],yHat2[,-1]),main="testing",type="l")
}
