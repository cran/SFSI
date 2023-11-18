
SSI <- function(y, X = NULL, b = NULL, Z = NULL, K,
                trn_tst = NULL, varU = NULL, varE = NULL,
                intercept = TRUE, alpha = 1, lambda = NULL,
                nlambda = 100, lambda.min = .Machine$double.eps^0.5,
                common.lambda = TRUE, subset = NULL, tol = 1E-4,
                maxiter = 500, method = c("REML","ML"),
                name = NULL, save.at = NULL, mc.cores = 1L,
                precision.format = c("single","double"),
                verbose = TRUE)
{
  method <- match.arg(method)
  precision.format <- match.arg(precision.format)
  # y=Y0[,k]; K=G0; method="REML"; precision.format='double'
  # varU=varU0; varE=varE0; subset=c(2,5); save.at=prefix
  if(length(dim(y)) == 2L){
    y <- as.matrix(y)
  }else{
    y <- matrix(y, ncol=1L)
  }
  n <- nrow(y)
  ntraits <- ncol(y)   # Number of traits

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

  trn <- which(as.vector(trn_tst))

  if(any(!trn_tst)){  # If any testing set
    tst <- which(as.vector(!trn_tst))
  }else{
    if(verbose){
      message(" No testing set was found. The SSI will be fitted to the entire data set")
    }
    tst <- trn[]
  }
  nTST <- length(tst)

  # Design matrix for fixed effects
  BLUE <- ifelse(is.null(X) & !intercept, FALSE, TRUE)
  if(verbose & !BLUE){
    message(" No intercept is estimated. Response is assumed to have mean zero")
  }
  X <- setX(n=n, X=X)
  K <- setK(n=n, Z=Z, K=K)

  labels <- NULL
  if(has_names(K)){
    labels <- rownames(K)
  }

  if(ntraits == 1L)  # Single-trait case
  {
    if(is.null(varU) | is.null(varE))
    {
      # Get variance components and estimate fixed effects
      res <- fitBLUP(y[trn,], X=X[trn, ,drop=FALSE], K=K[trn,trn],
                     intercept=intercept, method=method,
                     BLUP=FALSE, verbose=FALSE)
      if(res$convergence){
        varU <- res$varU
        varE <- res$varE
        b <- res$b

        if(verbose){
          message(" 'varU' and 'varE' were estimated using a BLUP model")
        }
      }else{
        stop("Convergence was not reached in the 'GEMMA' algorithm.\n",
             "       Provide variance components' estimates")
      }
    }else{   # Only estimate fixed effects as GLS
      if(is.null(b) & BLUE){
        b <- fitBLUP(y[trn,], X=X[trn, ,drop=FALSE], K=K[trn,trn],
                     BLUP=FALSE, varU=varU, varE=varE, verbose=FALSE)$b
      }
    }

  }else{   # Multitrait case
    trn0 <- which(apply(trn_tst,1,all))  # Common TRN set
    if(is.null(varU) | is.null(varE)){
      if(length(trn0) == 0){
         stop("No common training to all response variables was found.\n",
              "       Please provide 'varU' and 'varE' estimates")
      }
      res <- getGenCov(y[trn0,], X=X[trn0, ,drop=FALSE], K=K[trn0,trn0],
                       pairwise=TRUE, intercept=intercept, verbose=FALSE)
      varU <- res$varU
      varE <- res$varE
      b <- res$b

      if(verbose){
        message(" 'varU' and 'varE' matrices were pairwise estimated using BLUP",
                " with n=",length(trn0)," common records")
      }

    }else{
      if(is.null(b) & BLUE){
        if(length(trn0) == 0){
           stop("No common training for all response variables was found.\n",
                "       Provide an estimate for 'b'")
        }
        b <- fitBLUP(y[trn0,], X=X[trn0, ,drop=FALSE], K=K[trn0,trn0], BLUP=FALSE,
                     varU=diag(varU), varE=diag(varE), verbose=FALSE)$b
      }
    }

    if(length(dim(varU)) != 2L | length(dim(varE)) != 2L){
      stop("Inputs 'varU' and 'varE' must be matrices of dimension equal to ncol(y)=",ncol(y))
    }
  }

  # Getting K <- (varU*G)[trn,tst] and H <- (varU*G + varE*I)[trn,trn]
  # for multitrait varU*G : Kronecker(varU,K) and varE*I : Kronecker(varE,I))
  H <- tensorEVD::Kronecker_cov(K, Sigma=t(varU), Theta=varE, rows=trn, cols=trn)
  K <- tensorEVD::Kronecker(varU, K, rows=trn, cols=tst)

  h2 <- varU/(varU + varE)
  if(any(diag(t(h2)) < 0.001) & verbose){
     message(" The 'heritability' is too small. Results may be affected")
  }

  # Adjusted training phenotypes
  Xb <- NULL
  if(BLUE){
    if(ifelse(ntraits==1L,length(b),nrow(b)) != ncol(X)){
       stop("Number of fixed effects 'b' must be the same as the number of columns of 'X'")
    }
    if(ntraits==1L){
      Xb <- as.vector(X%*%b)
    }else{
      for(k in 1:ntraits){
        Xb <- c(Xb, as.vector(X%*%b[,k]))
      }
    }
    yTRN <- matrix(y[trn]-Xb[trn], nrow=1)
  }else{
    yTRN <- matrix(y[trn], nrow=1)
  }

  # Standardize matrices
  sdx <-  sqrt(diag(H))
  cov2cor2(H, inplace = TRUE)
  K <- sweep(K, 1L, sdx, FUN = "/")

  lambda <- setLambda(K, alpha=alpha, lambda=lambda, nlambda=nlambda,
                      lambda.min=lambda.min, lambda.max=NULL,
                      common.lambda=common.lambda)

  name <- ifelse(is.null(name),"SSI",name)

  # Split the testing set into subsets. Only the subset provided will be fitted
  if(is.null(subset)){
    fileID <- NULL
    tmp <- ""
  }else{
     if(!is.numeric(subset) & length(subset) != 2L){
       stop("Input 'subset' must be a 2-elements vector")
     }
     sets <- sort(rep(1:subset[2],ceiling(nTST/subset[2]))[1:nTST])
     index <- which(sets == subset[1])
     fileID <- index[]
     tmp <- paste0(" of ",length(tst))
     tst <- tst[index]
     K <- K[,index,drop=FALSE]
     if(ncol(lambda)==nTST){
       lambda <- lambda[,index,drop=FALSE]
     }
  }

  if(verbose){
    message(" Fitting a ",ifelse(ntraits==1L,"SSI",paste("Multi-trait SSI for",ntraits,"traits")),
            " with nTST=",length(tst),tmp," and nTRN=",length(trn))
  }

  # If 'save.at' is not NULL
  if(!is.null(save.at)){
    stopifnot(is.character(save.at))
    save.at <- normalizePath(save.at, mustWork=F)
    prefix <- basename(tempfile(pattern=""))
    if(is.null(subset)){
      outfile <- paste0(save.at,"output.RData")
    }else{
      prefix <- paste0("subset_",subset[1],"_of_",subset[2],"_")
      outfile <- paste0(save.at,prefix,"output.RData")
    }
    unlink(outfile)
  }

  out <- solveEN(Sigma=H, Gamma=K, alpha=alpha, lambda=lambda,
                 nlambda=nlambda, lambda.min=lambda.min,
                 common.lambda=common.lambda, tol=tol,
                 maxiter=maxiter, save.at=save.at, fileID=fileID,
                 verbose=verbose, mc.cores=mc.cores,
                 precision.format=precision.format,
                 scale=FALSE, sdx=sdx)

  if(length(tst) == 1L){
    nsup0 <- matrix(out$nsup, nrow=1)
    lambda0 <- matrix(out$lambda, nrow=1)
  }else{
    nsup0 <- do.call(rbind, out$nsup)
    lambda0 <- do.call(rbind, out$lambda)
  }

  u <- fitted.LASSO(out, yTRN)
  dimnames(u) <- list(tst, paste0("SSI.",1:ncol(u)))
  if(BLUE){
    yHat <- sweep(u, 1L, Xb[tst], FUN="+")
  }else{
    yHat <- u[]
  }

  out <- list(n=n, p=out$p, q=out$q, ntraits=ntraits,
              labels=labels, name=name, nlambda=nlambda,
              y=y, Xb=Xb, u=u, yHat=yHat,
              b=b, varU=varU, varE=varE, h2=h2,
              trn=trn, tst=tst, alpha=alpha,
              nsup = nsup0, lambda = lambda0,
              beta = out$beta,
              file_beta=out$file_beta,
              fileID=out$fileID
              #precision.format=precision.format
            )

  class(out) <- c("SSI")

  # Save outputs if 'save.at' is not NULL
  if(!is.null(save.at)){
    save(out, file=outfile)
  }

  return(out)
}
