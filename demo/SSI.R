
oldwd <- getwd()
setwd(tempdir())

library(SFSI)
data(wheatHTP)

index = which(Y$trial %in% 1:6)      # Use only a subset of data
Y = Y[index,]
M = scale(M[index,])/sqrt(ncol(M))   # Subset and scale markers
G = tcrossprod(M)                    # Genomic relationship matrix
y = as.vector(scale(Y[,'E1']))       # Subset response variable

# Calculate variance components ratio using all data
fm1 = fitBLUP(y,K=G)
varU = fm1$varU
varE = fm1$varE
h2 = varU/(varU + varE)

# Sparse selection index
fm2 = SSI(y,K=G,varU=varU,varE=varE,nlambda=50)

# The same but without passing variance components
fm2 = SSI(y,K=G,nlambda=50)
u2 = fitted(fm2)

plot(fm2)  # Penalization vs accuracy

# Equivalence of the SSI with lambda=0 with G-BLUP
fm3 = SSI(y,K=G,varU=varU,varE=varE,lambda=0,tol=1E-5)

cor(y, fm1$u)        # G-BLUP accuracy
cor(y, fitted(fm3))  # SSI accuracy
cor(fm1$u, fitted(fm3))

# Predicting a testing set using training set
# 0: testing, 1:training
trn_tst = ifelse(Y$trial %in% 2, 0, 1)

# Calculate variance components in training data
yNA = y
yNA[trn_tst==0] = NA
fm1 = fitBLUP(yNA,K=G)
(varU = fm1$varU)
(varE = fm1$varE)
(h2 = varU/(varU + varE))

# Sparse selection index
fm2 = SSI(y,K=G,varU=varU,varE=varE,trn_tst=trn_tst)

# Variance components internaly calculated
fm2 = SSI(y,K=G,trn_tst=trn_tst)
fm2$varU
fm2$varE
fm2$h2

# Effect of the penalization on the accuracy
plot(fm2)

setwd(oldwd)
