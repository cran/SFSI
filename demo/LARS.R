oldwd <- getwd()
setwd(tempdir())

library(SFSI)
data(wheatHTP)
y = as.vector(Y[,"E1"])      # Response variable
X = scale(X_E1)              # Predictors

# Training and testing sets
tst = seq(1,length(y),by=3)
trn = seq_along(y)[-tst]

# Calculate covariances in training set
XtX = var(X[trn,])
Xty = cov(y[trn],X[trn,])

# Run the penalized regression
fm = LARS(XtX,Xty,method="LAR")
fm = LARS(XtX,Xty,method="LAR-LASSO",verbose=TRUE)

# Predicted values
yHat1 = fitted(fm, X=X[trn,])  # training data
yHat2 = fitted(fm, X=X[tst,])  # testing data

# Penalization vs correlation
oldpar <- par(mfrow=c(1,2))
plot(-log(fm$lambda),cor(y[trn],yHat1)[1,], main="training")
plot(-log(fm$lambda),cor(y[tst],yHat2)[1,], main="testing")

par(oldpar)
setwd(oldwd)
