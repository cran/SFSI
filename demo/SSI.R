
oldwd <- getwd()
setwd(tempdir())

library(SFSI)
data(wheatHTP)

index = which(Y$CV %in% 1:2)
M = scale(M[index,])/sqrt(ncol(M))   # Subset and scale markers
G = tcrossprod(M)                    # Genomic relationship matrix
y = as.vector(scale(Y[index,"E1"])) # Subset response variable

# Calculate variance components ratio using all data
fm1 = fitBLUP(y,K=G)
theta = fm1$varE/fm1$varU
h2 = fm1$varU/(fm1$varU + fm1$varE)

# Sparse selection index
fm2 = SSI(y,K=G,theta=theta,nLambda=50)

# The same but passing the heritability instead of theta
fm2 = SSI(y,K=G,h2=h2,nLambda=50)
yHat = fitted(fm2)

plot(fm2)  # Penalization vs accuracy

# Equivalence of the SSI with lambda=0 with G-BLUP
fm3 = SSI(y,K=G,theta=theta,lambda=0,tol=1E-5)

cor(y,fm1$u)        # G-BLUP accuracy
cor(y,fitted(fm3))  # SSI accuracy

# Predicting a testing set using training set
tst = seq(1,length(y),by=3)
trn = (seq_along(y))[-tst]

# Calculate variance components in training data
yNA = y
yNA[tst] = NA
fm1 = fitBLUP(yNA,K=G)
(theta = fm1$varE/fm1$varU)
(h2 = fm1$varU/(fm1$varU + fm1$varE))

# Sparse selection index
fm2 = SSI(y,K=G,theta=theta,trn=trn,tst=tst)

# Variance components ratio internaly calculated
fm2 = SSI(y,K=G,theta=NULL,trn=trn,tst=tst)
fm2$theta
fm2$h2

# Effect of the penalization on the accuracy
plot(fm2)

setwd(oldwd)
