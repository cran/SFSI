oldwd <- getwd()
setwd(tempdir())

library(SFSI)
data(wheatHTP)
X = scale(X[1:400,])/sqrt(ncol(X))   # Subset and scale markers
G = tcrossprod(X)                    # Genomic relationship matrix
y = as.vector(scale(Y[1:400,"YLD"])) # Subset response variable

# Predicting a testing set using training set
tst = sample(seq_along(y),ceiling(0.3*length(y)))
trn = (seq_along(y))[-tst]

# Obtain lambda from cross-validation (in traning set)
fm1 = SSI_CV(y,K=G,trn=trn,nFolds=5,nCV=2)
lambda = summary(fm1)$optCOR["mean","lambda"]

# Fit the index with the obtained lambda
fm2 = SSI(y,K=G,h2=NULL,trn=trn,tst=tst,lambda=lambda)
summary(fm2)$accuracy        # Testing set accuracy

# Compare the accuracy with that of the non-sparse index
fm3 = SSI(y,K=G,h2=NULL,trn=trn,tst=tst,lambda=0)
summary(fm3)$accuracy

setwd(oldwd)
