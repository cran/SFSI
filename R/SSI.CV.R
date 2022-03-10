 #y=y; X = NULL; Z = D= NULL; K=G; b=NULL; h2 = NULL; trn = trn; alpha = 1;
 #lambda = NULL; nlambda = 100; lambda.min = .Machine$double.eps^0.5; nCV = 1; nfolds = 5;
 #seed = NULL; common.lambda = TRUE; method = c("REML","ML")[1]
 #mc.cores = 1; tol = 1E-4; maxiter = 500; name = NULL; verbose = TRUE

SSI.CV <- function(y, X = NULL, b = NULL, Z = NULL, K, D = NULL,
              theta = NULL, h2 = NULL, trn = seq_along(y), alpha = 1,
              lambda = NULL, nlambda = 100, lambda.min = .Machine$double.eps^0.5,
              nCV = 1, nfolds = 5, seed = NULL, common.lambda = TRUE,
              tol = 1E-4, maxiter = 500, method = c("REML","ML"),
              name = NULL, mc.cores = 1, verbose = TRUE)
{
    method <- match.arg(method)

    tmp <- c(2,3,4,5,10)
    if(!as.character(nfolds) %in% c(tmp,'n')) stop("'nfolds' should be one of ",paste(tmp,collapse=","),",'n'")

    if(is.logical(trn)){
       if(length(y) != length(trn)) stop("Object 'trn' must be of the same length of 'y'")
       trn <- which(trn)
    }
    if(is.character(K)){
        K <- readBinary(K)
    }

    if(!float::storage.mode(K) %in% c("float32","double")) storage.mode(K) <- "double"

    if(!is.null(Z)){
      if(length(dim(Z)) != 2) stop("Object 'Z' must be a matrix with nrow(Z)=n and ncol(Z)=nrow(K)")
      K <- float::tcrossprod(Z,float::tcrossprod(Z,K))   # Z%*%K%*%t(Z)
    }

    name <- ifelse(is.null(name),"SSI.CV",name)
    nTRN <- length(trn)
    isLOOCV <- nfolds=='n'
    nfolds <- ifelse(isLOOCV,nTRN,as.numeric(nfolds))
    mc.cores2 <- ifelse(isLOOCV,1,mc.cores)

    compApply <- function(ind)
    {
      trn0 <- trn[folds != ind]
      tst0 <- trn[folds == ind]

      fm <- SSI(y, X=X, b=b, K=K, D=D, theta=theta, h2=h2, trn=trn0, tst=tst0,
              alpha=alpha, method=method, lambda=lambda, nlambda=nlambda,
              lambda.min=lambda.min, tol=tol, maxiter=maxiter,
              common.lambda=common.lambda, mc.cores=mc.cores2, verbose=FALSE)

      if(isLOOCV){
          rr <- list(u=as.vector(fitted.SSI(fm)), varU=fm$varU, varE=fm$varE,
                     h2=fm$h2, theta=fm$theta, b=fm$b, tst=tst0, df=fm$df, lambda=fm$lambda)
      }else{
          fv <- summary.SSI(fm)
          rr <- list(varU=fm$varU, varE=fm$varE, h2=fm$h2, theta=fm$theta, b=fm$b, tst=tst0,
                     df=fv$df, lambda=fv$lambda, accuracy=fv$accuracy, MSE=fv$MSE)
      }

      if(verbose){
        message("CV ",ifelse(isLOOCV,"LOO",k),": Fold ",ind,"/",nfolds," (n=",length(tst0),")")
      }
      rr
    }

    if(is.null(seed)){   # Seeds for randomization
      seeds <- round(seq(1E3, .Machine$integer.max/1000, length=nCV))
      #if(nCV >1) seeds <- sample(seeds)
    }else{
      seeds <- seed
      nCV <- length(seeds)
    }
    stopifnot(nCV == length(seeds))
    nCV <- ifelse(isLOOCV & nCV > 1,1, nCV)   # No need of repeating CV when LOO

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

      if(mc.cores == 1L | !isLOOCV){
        out = lapply(X=seq(nfolds),FUN=compApply)
      }
      if(mc.cores > 1L & isLOOCV){
        out = parallel::mclapply(X=seq(nfolds),FUN=compApply,mc.cores=mc.cores)
      }

      tmp <- do.call(rbind,split(data.frame(trn,folds),folds))[,1]
      if(sum(tmp != unlist(lapply(out,function(x)x$tst)))>0){
        stop("Some sub-processes failed. Something went wrong during the analysis.")
      }

      # Calculate accuracy
      if(isLOOCV){
        uHat <- do.call(rbind,lapply(out,function(x)x$u))
        accuracy <- suppressWarnings(stats::cor(y[trn],uHat,use="pairwise.complete.obs"))
        MSE <- suppressWarnings(t(apply((y[trn]-uHat)^2,2,sum,na.rm=TRUE)/length(trn)))
        df <- t(apply(do.call("rbind",lapply(out,function(x)x$df)),2,mean))
        lambda0 <- t(apply(do.call("rbind",lapply(out,function(x)x$lambda)),2,mean))
        b0 <- t(apply(do.call("rbind",lapply(out,function(x)x$b)),2,mean))
        varU0 <- mean(unlist(lapply(out,function(x)x$varU)))
        varE0 <- mean(unlist(lapply(out,function(x)x$varE)))
        h20 <- mean(unlist(lapply(out,function(x)x$h2)))
        theta0 <- mean(unlist(lapply(out,function(x)x$theta)))

      }else{
        accuracy <- do.call("rbind",lapply(out,function(x)x$accuracy))
        MSE <- do.call("rbind",lapply(out,function(x)x$MSE))
        df <- do.call("rbind",lapply(out,function(x)x$df))
        lambda0 <- do.call("rbind",lapply(out,function(x)x$lambda))
        b0 <- do.call("rbind",lapply(out,function(x)x$b))
        varU0 <- unlist(lapply(out,function(x)x$varU))
        varE0 <- unlist(lapply(out,function(x)x$varE))
        h20 <- unlist(lapply(out,function(x)x$h2))
        theta0 <- unlist(lapply(out,function(x)x$theta))
      }

      res[[k]] <- list(folds=data.frame(trn,fold=folds),
                   b=b0, varU=varU0, varE=varE0, h2=h20, theta=theta0, accuracy=accuracy,
                   MSE=MSE, df=df, lambda=lambda0)
    }

    out <- list(name=name, y=y, trn=trn, alpha=alpha,
                nCV=nCV, nfolds=nfolds, seeds=seeds, CV=res)
    class(out) <- "SSI"

    return(out)
}
