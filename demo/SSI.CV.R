
oldwd <- getwd()
setwd(tempdir())

library(SFSI)
data(wheatHTP)

index = which(Y$trial %in% 1:6)      # Use only a subset of data
Y = Y[index,]
M = scale(M[index,])/sqrt(ncol(M))   # Subset and scale markers
G = tcrossprod(M)                    # Genomic relationship matrix
y = as.vector(scale(Y[,'E1']))       # Subset response variable

# Predicting a testing set using training set
trn_tst = ifelse(Y$trial %in% 2, 0, 1)

# Obtain lambda from cross-validation (in traning set)
fm1 = SSI.CV(y,K=G,trn_tst=trn_tst,nfolds=5,nCV=3)
lambda = summary(fm1)$optCOR["lambda"]

# Fit the index with the obtained lambda
fm2 = SSI(y,K=G,trn_tst=trn_tst,lambda=lambda)
summary(fm2)$accuracy        # Testing set accuracy

# Compare the accuracy with that of the non-sparse index
fm3 = SSI(y,K=G,trn_tst=trn_tst,lambda=0)
summary(fm3)$accuracy

setwd(oldwd)
