
oldwd <- getwd()
setwd(tempdir())

library(SFSI)
data(wheatHTP)

X = scale(X[1:400,])/sqrt(ncol(X))  # Subset and scale markers
G = tcrossprod(X)                  # Genomic relationship matrix
y = scale(Y[1:400,"YLD"])           # Subset response variable

# Fit model, whole data
fm1 = fitBLUP(y,K=G)
fm1$varU
fm1$varE
fm1$h2
cor(y,fm1$u)                 # Prediction accuracy

# Same model different parametrization
fm2 = fitBLUP(y,Z=X)
fm2$varU; fm2$varE; fm2$h2
cor(y,X%*%fm2$u)             # Prediction accuracy

# Training and testing sets
tst = sample(seq_along(y),ceiling(0.3*length(y)))
trn = seq_along(y)[-tst]

yNA <- y
yNA[tst] <- NA

# Fit model, split data
fm3 = fitBLUP(yNA,K=G)
plot(y[tst],fm3$u[tst])      # Predicted vs observed values in testing set
cor(y[tst],fm3$u[tst])       # Prediction accuracy in testing set
cor(y[trn],fm3$u[trn])       # Prediction accuracy in training set
fm3$h2                       # Heritability (in training set)

setwd(oldwd)
