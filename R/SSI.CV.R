
SSI.CV <- function(y, X = NULL, b = NULL, Z = NULL, K,
              trn_tst = NULL, varU = NULL, varE = NULL,
              intercept = TRUE, alpha = 1,
              lambda = NULL, nlambda = 100,
              lambda.min = .Machine$double.eps^0.5,
              common.lambda = TRUE, nCV = 1L, nfolds = 5,
              seed = NULL, tol = 1E-4, maxiter = 500,
              method = c("REML","ML"), name = NULL,
              mc.cores = 1L, verbose = TRUE)
{
  method <- match.arg(method)

  # K=G0; trn_tst=trn_tst0; mc.cores=6; method="REML"
  tmp <- c(2,3,4,5,10)
  if(!as.character(nfolds) %in% c(tmp,'n')){
    stop("'nfolds' should be one of ",paste(tmp,collapse=","),",'n'")
  }

  if(length(dim(y)) == 2L){
    y <- as.matrix(y)
  }else{
    y <- matrix(y, ncol=1L)
  }
  dimnames(y) <- NULL
  n <- nrow(y)
  ntraits <- ncol(y)

  if(is.null(trn_tst)){
     trn_tst <- matrix(TRUE, nrow=n, ncol=ntraits)
  }else{
     if(length(dim(trn_tst)) == 2L){
       trn_tst <- as.matrix(trn_tst)
     }else{
       trn_tst <- matrix(trn_tst, ncol=1L)
     }
     stopifnot(length(trn_tst) == (n*ntraits))
     if(storage.mode(trn_tst) %in% c("double","integer")){
       if(any( !unique(as.vector(trn_tst)) %in% c(0,1,NA) )){
         stop("Input 'trn_tst' must contain only 0, 1, or NA")
       }
       trn_tst <- (trn_tst == 1)
     }
     stopifnot(storage.mode(trn_tst) == "logical")
  }
  dimnames(trn_tst) <- NULL

  storage.mode(K) <- "double"

  if(!is.null(Z)){
    if(length(dim(Z)) != 2L){
       stop("Object 'Z' must be a matrix with nrow(Z)=n and ncol(Z)=nrow(K)")
    }
    K <- tcrossprod(Z, tcrossprod(Z,K))   # Z%*%K%*%t(Z)
  }

  trn <- which(as.vector(trn_tst))
  tst <- which(as.vector(!trn_tst))
  MAP <- map_set(n, ntraits, x=trn, y=tst, xlab="trn", ylab="tst")

  name <- ifelse(is.null(name),"SSI.CV",name)
  nTRN <- length(trn)
  isLOOCV <- nfolds=='n'
  nfolds <- ifelse(isLOOCV,nTRN,as.numeric(nfolds))
  mc.cores2 <- ifelse(isLOOCV,1,mc.cores)

  for(k in 1:ncol(trn_tst)){  # Ignore tst entries by setting NA
    tmp <- which(!trn_tst[,k])
    if(length(tmp)>0) trn_tst[tmp,k] <- NA
  }

  compApply <- function(ind)
  {
    trn_tst0 <- trn_tst[]
    tst0 <- trn[folds == ind]
    map <- MAP[tst0,]
    for(i in 1:nrow(map)){
      trn_tst0[map$i[i],map$j[i]] <- FALSE
    }

    fm <- SSI(y, X=X, b=b, K=K, varU=varU, varE=varE,
              intercept=intercept, trn_tst=trn_tst0, alpha=alpha,
              method=method, lambda=lambda, nlambda=nlambda,
              lambda.min=lambda.min, tol=tol, maxiter=maxiter,
              common.lambda=common.lambda, mc.cores=mc.cores2,
              verbose=FALSE)

    if(isLOOCV){
        res <- list(u=as.vector(fitted.SSI(fm)), varU=fm$varU,
                    varE=fm$varE, h2=fm$h2, b=fm$b, tst=tst0,
                    nsup=fm$nsup, lambda=fm$lambda)
    }else{
        ss <- summary.SSI(fm)
        res <- list(varU=fm$varU, varE=fm$varE, h2=fm$h2, b=fm$b,
                    tst=tst0, nsup=ss$nsup, lambda=ss$lambda,
                    nsup_trait=ss$nsup_trait, accuracy=ss$accuracy, MSE=ss$MSE)
    }

    if(verbose){
      message(" CV ",ifelse(isLOOCV,"LOO",k),": Fold ",ind,"/",nfolds," (n=",length(tst0),")")
    }
    res
  }

  if(is.null(seed)){   # Seeds for randomization
    seeds <- round(seq(1E3, .Machine$integer.max/1000, length=nCV))
  }else{
    seeds <- seed
    nCV <- length(seeds)
  }
  stopifnot(nCV == length(seeds))
  nCV <- ifelse(isLOOCV & nCV > 1, 1, nCV)   # No need of repeating CV when LOO

  res <- vector("list",nCV)
  names(res) <- paste0("CV",1:nCV)
  for(k in 1:nCV)
  {
    # Creating folds
    folds <- rep(seq(1:nfolds),ceiling(nTRN/nfolds))[1:nTRN]
    if(!isLOOCV){
      set.seed(seeds[k])
      folds <- sample(folds)
    }

    if(mc.cores > 1L & isLOOCV){
      out <- parallel::mclapply(X=seq(nfolds),FUN=compApply,mc.cores=mc.cores)
    }else{
      out <- lapply(X=seq(nfolds),FUN=compApply)
    }

    tmp <- do.call(rbind,split(data.frame(trn,folds),folds))[,1]
    if(sum(tmp != unlist(lapply(out,function(x)x$tst)))>0){
      stop("Some sub-processes failed. Something went wrong during the analysis.")
    }

    # Calculate accuracy
    nsup_trait <- NULL
    if(isLOOCV){
      nsup_trait <- NULL
      uHat <- do.call(rbind,lapply(out,function(x)x$u))
      accuracy <- suppressWarnings(stats::cor(y[trn],uHat,use="pairwise.complete.obs"))
      MSE <- suppressWarnings(t(apply((y[trn]-uHat)^2,2,sum,na.rm=TRUE)/length(trn)))
      nsup <- t(apply(do.call("rbind",lapply(out,function(x)x$nsup)),2,mean))
      lambda0 <- t(apply(do.call("rbind",lapply(out,function(x)x$lambda)),2,mean))
      b0 <- t(apply(do.call("rbind",lapply(out,function(x)x$b)),2,mean))
      varU0 <- mean(unlist(lapply(out,function(x)x$varU)))
      varE0 <- mean(unlist(lapply(out,function(x)x$varE)))
      h20 <- mean(unlist(lapply(out,function(x)x$h2)))

    }else{
      accuracy <- lapply(out,function(x)x$accuracy)
      MSE <- lapply(out,function(x)x$MSE)
      nsup <- lapply(out,function(x)x$nsup)
      lambda0 <- lapply(out,function(x)x$lambda)
      nsup_trait <- lapply(out,function(x)x$nsup_trait)
      b0 <- do.call("rbind",lapply(out,function(x)x$b))
      varU0 <- lapply(out,function(x)x$varU)
      varE0 <- lapply(out,function(x)x$varE)
      h20 <- lapply(out,function(x)x$h2)
    }

    res0 <- list(folds=data.frame(trn,fold=folds),
                 b=b0, varU=varU0, varE=varE0, h2=h20, accuracy=accuracy,
                 MSE=MSE, nsup=nsup, lambda=lambda0)
    if(ntraits>1) res0$nsup_trait <- nsup_trait

    res[[k]] <- res0
  }

  out <- list(n=n, ntraits=ntraits, nlambda=nlambda, name=name, trn=trn,
              alpha=alpha, nCV=nCV, nfolds=nfolds, seeds=seeds, CV=res)

  class(out) <- "SSI"

  return(out)
}
