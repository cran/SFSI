## ----eval=FALSE---------------------------------------------------------------
#  library(SFSI)
#  if(requireNamespace("BGLR")){
#   data(wheat, package="BGLR")     # Load data from the BGLR package
#  }
#  
#  # Select the environment 1 to work with
#  y <- as.vector(scale(wheat.Y[,1]))
#  
#  # Calculate G matrix
#  G <- tcrossprod(scale(wheat.X))/ncol(wheat.X)
#  
#  # Save data
#  save(y, G, file="geno_pheno.RData")

## ----eval=FALSE---------------------------------------------------------------
#  load("geno_pheno.RData") # Load data
#  
#  # Fit model
#  fm0 <- fitBLUP(y, K=G)
#  fm0$theta <- fm0$varE/fm0$varU          # Residual/genetic variances ratio
#  fm0$h2 <- fm0$varU/(fm0$varU+fm0$varE)  # Heritability
#  c(fm0$varU,fm0$varE,fm0$theta,fm0$h2)   # Print variance components
#  
#  save(fm0, file="varComps.RData")

## ----eval=FALSE---------------------------------------------------------------
#  nPart <- 5                        # Number of partitions
#  load("geno_pheno.RData")          # Load data
#  
#  nTST <- ceiling(0.3*length(y))    # Number of elements in TST set
#  partitions <- matrix(1,nrow=length(y),ncol=nPart)   # Matrix to store partitions
#  seeds <- round(seq(1E3, .Machine$integer.max/10, length=nPart))
#  
#  for(k in 1:nPart)
#  {   set.seed(seeds[k])
#      partitions[sample(1:length(y), nTST),k] <- 2
#  }
#  save(partitions, file="partitions.RData")    # Save partitions

## ----eval=FALSE---------------------------------------------------------------
#  # Load data
#  load("geno_pheno.RData"); load("varComps.RData"); load("partitions.RData")
#  
#  accSSI <- mu <- h2 <- c()       # Objects to store results
#  
#  for(k in 1:ncol(partitions))
#  { cat("  partition = ",k,"\n")
#    trn <- which(partitions[,k] == 1)
#    tst <- which(partitions[,k] == 2)
#    yNA <- y;   yNA[tst] <- NA
#  
#    # G-BLUP model
#    fm1 <- fitBLUP(yNA, K=G)
#    mu[k] <- fm1$b        # Retrieve mu estimate
#    h2[k] <- fm1$h2       # Retrieve h2 estimate
#  
#    # Sparse SI
#    fm2 <- SSI(y,K=G,b=mu[k],h2=h2[k],trn=trn,tst=tst,mc.cores=1,nLambda=100)
#    fm3 <- summary(fm2)   # Useful function to get results
#  
#    accuracy <- c(GBLUP=cor(fm1$u[tst],y[tst]), fm3$accuracy)/sqrt(fm0$h2)
#    lambda <- c(min(fm3$lambda),fm3$lambda)
#    df <- c(max(fm3$df),fm3$df)
#    accSSI <- rbind(accSSI,data.frame(rep=k,SSI=names(accuracy),accuracy,lambda,df))
#  }
#  save(mu,h2,accSSI,file="results_accuracy.RData")

## ----eval=FALSE---------------------------------------------------------------
#  load("results_accuracy.RData")
#  
#  dat <- data.frame(do.call(rbind,lapply(split(accSSI,accSSI$SSI),
#           function(x) apply(x[,-c(1:2)],2,mean))))
#  dat$Model <- unlist(lapply(strsplit(rownames(dat),"\\."),function(x)x[1]))
#  
#  dat2 <- do.call(rbind,lapply(split(dat,dat$Mod),function(x)x[which.max(x$acc),]))
#  
#  if(requireNamespace("ggplot2")){
#   ggplot2::ggplot(dat[dat$df>1,],ggplot2::aes(-log(lambda),accuracy)) +
#     ggplot2::geom_hline(yintercept=dat["GBLUP",]$accuracy, linetype="dashed") +
#     ggplot2::geom_line(ggplot2::aes(color=Model),size=1.1) + ggplot2::theme_bw() +
#     ggplot2::geom_point(data=dat2,ggplot2::aes(color=Model),size=2.5)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  load("geno_pheno.RData");   load("varComps.RData")
#  load("partitions.RData");   load("results_accuracy.RData")
#  
#  lambdaCV <- accSSI_CV <- dfCV <- c()   # Objects to store results
#  
#  for(k in 1:ncol(partitions))
#  {   cat("  partition = ",k,"\n")
#      trn <- which(partitions[,k] == 1)
#      tst <- which(partitions[,k] == 2)
#  
#      # Cross-validating the training set
#      fm1 <- SSI_CV(y,K=G,trn=trn,nLambda=100,mc.cores=1,nFolds=5,nCV=1)
#      lambdaCV[k] <- summary(fm1)$optCOR["mean","lambda"]
#  
#      # Fit a SSI with the estimated lambda
#      fm2 <- SSI(y,K=G,b=mu[k],h2=h2[k],trn=trn,tst=tst,lambda=lambdaCV[k])
#  
#      accSSI_CV[k] <- summary(fm2)$accuracy/sqrt(fm0$h2)
#      dfCV <- cbind(dfCV, fm2$df)
#  }
#  save(accSSI_CV,lambdaCV,dfCV,file="results_accuracyCV.RData")

## ----eval=FALSE---------------------------------------------------------------
#  load("results_accuracy.RData"); load("results_accuracyCV.RData")
#  
#  dat <- data.frame(GBLUP=accSSI[accSSI$SSI=="GBLUP",]$acc,SSI=accSSI_CV)
#  rg <- range(dat)
#  tmp <- c(mean(rg),diff(rg)*0.4)
#  
#  if(requireNamespace("ggplot2")){
#   ggplot2::ggplot(dat,ggplot2::aes(GBLUP,SSI)) +
#    ggplot2::geom_abline(slope=1,linetype="dotted") + ggplot2::theme_bw() +
#    ggplot2::geom_point(shape=21,color="orange") + ggplot2::lims(x=rg,y=rg) +
#    ggplot2::annotate("text",tmp[1],tmp[1]-tmp[2],label=round(mean(dat$GBLUP),3)) +
#    ggplot2::annotate("text",tmp[1]-tmp[2],tmp[1],label=round(mean(dat$SSI),3))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  load("results_accuracyCV.RData")
#  
#  dat <- data.frame(df=as.vector(dfCV))
#  
#  bw <- round(diff(range(dat$df))/40)
#  if(requireNamespace("ggplot2")){
#   ggplot2::ggplot(data=dat,ggplot2::aes(df,stat(count)/length(dfCV))) +
#      ggplot2::theme_bw() +
#      ggplot2::geom_histogram(color="gray20",fill="lightblue",binwidth=bw) +
#      ggplot2::labs(x=bquote("Support set size(" *n[sup]*")"),y="Frequency")
#  }

## ----eval=FALSE---------------------------------------------------------------
#  # Load data
#  load("geno_pheno.RData"); load("partitions.RData"); load("results_accuracyCV.RData")
#  
#  part <- 1      # Choose any partition from 1,…,nPart
#  trn <- partitions[,part] == 1
#  tst <- partitions[,part] == 2
#  
#  # Fit SSI with lambda previously estimated using CV
#  fm <- SSI(y,K=G,trn=trn,tst=tst,lambda=lambdaCV[part])
#  
#  plotNet(fm,K=G,tst=fm$tst[1:16],unified=FALSE,title=NULL,bg.col="white",
#          set.size=c(3,1.5,0.2),point.color="gray40",axis.labels=FALSE)

