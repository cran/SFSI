
oldwd <- getwd()
setwd(tempdir())

library(SFSI)
data(wheatHTP)
y = as.vector(Y[,"E1"])    # Response variable
X = scale(X_E1)            # Predictors

# Training and testing sets
tst = which(Y$trial %in% 1:10)
trn = seq_along(y)[-tst]

# Calculate covariances in training set
XtX = var(X[trn,])
Xty = cov(X[trn,],y[trn])

# Run the penalized regression
fm = solveEN(XtX,Xty)
fm = solveEN(XtX,Xty,alpha=0.5)

# Predicted values
yHat1 = fitted(fm, X=X[trn,])  # training data
yHat2 = fitted(fm, X=X[tst,])  # testing data

# Penalization vs correlation
oldpar <- par(mfrow=c(1,2))
plot(-log(fm$lambda),cor(y[trn],yHat1)[1,], main="training")
plot(-log(fm$lambda),cor(y[tst],yHat2)[1,], main="testing")

par(oldpar)
setwd(oldwd)
