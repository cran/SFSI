
SSI <- function(y, X = NULL, b = NULL, Z = NULL, K, D = NULL,
         theta = NULL, h2 = NULL, trn = seq_along(y),
         tst = seq_along(y), subset = NULL, alpha = 1, lambda = NULL,
         nlambda = 100, lambda.min = .Machine$double.eps^0.5,
         common.lambda = TRUE, tol = 1E-4, maxiter = 500,
         method = c("REML","ML"), return.beta = FALSE, save.beta = TRUE,
         save.at = NULL, name = NULL, mc.cores = 1L, verbose = TRUE)
{
  method <- match.arg(method)

  y <- as.vector(y)
  n <- length(y)
  if(is.logical(trn)){
     if(n != length(trn)) stop("Object 'trn' must be of the same length of 'y'\n")
     trn <- which(trn)
  }
  if(is.logical(tst)){
     if(n != length(tst)) stop("Object 'tst' must be of the same length of 'y'\n")
     tst <- which(tst)
  }
  nTST <- length(tst)
  if(any(is.na(y[trn]))){
    stop("All entries in y[trn] must be non-NA")
  }

  if(is.character(K)){
      K <- readBinary(K)
  }

  if(!float::storage.mode(K) %in% c("float32","double")) storage.mode(K) <- "double"

  if(is.null(X))   # Design matrix for fixed effects including the intercept
  {
    X <- model.matrix(~1, data=data.frame(rep(1,n)))
  }else{
    if(length(dim(X)) < 2L){
      X <- stats::model.matrix(~X)
      if(ncol(X) > 2L)  colnames(X)[-1] <- substr(colnames(X)[-1],2,nchar(colnames(X)[-1]))
    }
  }

  if(!is.null(Z)){
    if(length(dim(Z)) != 2L) stop("Object 'Z' must be a matrix with nrow(Z)=n and ncol(Z)=nrow(K)\n")
    K <- float::tcrossprod(Z,float::tcrossprod(Z,K))   # Z%*%K%*%t(Z)
  }

  if(length(dim(K)) != 2L | (length(K) != n^2)){
    stop("Product Z %*% K %*% t(Z) must be a squared matrix with number of rows",
     "\n(and columns) equal to n: number of elements in 'y'")
  }
  id <- NULL
  if(has_names(K)){
    id <- rownames(K)
  }

  if(!is.null(D)){
    if((sum(dim(D))/2)^2 != n^2){
      stop("Object 'D' must be a nxn squared matrix with n: number of elements in 'y'")
    }
    dimnames(D) <- NULL
    D <- D[trn,trn]
  }

  RHS <- K[trn,tst, drop=FALSE]
  K <- K[trn,trn]

  if(is.null(theta) & is.null(h2))
  {
    # Fit LMM to get variance components and estimate fixed effects as GLS
    fm <- fitBLUP(y[trn], X=X[trn, ,drop=FALSE], K=K, method=method)
    if(fm$convergence){
      varU <- fm$varU
      varE <- fm$varE
      h2 <- varU/(varU + varE)
      theta <- varE/varU
      b <- fm$b
    }else stop("Convergence was not reached in the 'GEMMA' algorithm. \n\t",
           "Please provide a heritability estimate in 'h2' parameter")
  }else{   # Only estimate fixed effects as GLS
    varU <- varE <- NA
    if(!is.null(theta) & !is.null(h2)){
      message("Both 'theta' and 'h2' are provided. Only 'theta' will be considered")
      h2 <- NULL
    }
    if(is.null(theta)){
      theta <- (1-h2)/h2
    }else{
      h2 <- 1/(1+theta)
    }
    if(is.null(b)){
      b <- fitBLUP(y[trn], X=X[trn, ,drop=FALSE], K=K, BLUP=FALSE, theta=theta)$b
    }else{
      if(length(b) != ncol(X)) stop("The length of 'b' must be the same as the number of columns of 'X'\n")
    }

  }
  if(h2 < 0.001) warning("The 'heritability' is too small. Results may be affected",immediate.=TRUE)
  Xb <- drop(X%*%b)

  # Adjusted training phenotypes
  yTRN <- matrix(y[trn]-Xb[trn], nrow=1)

  # Adding theta*I
  if(is.null(D)){
    add2diag(K, a=theta, void=TRUE)   #diag(K) <- diag(K) + theta
  }else{
    K <- K + theta*D
  }

  if(is.null(lambda)){
    if(common.lambda){
      sdx <-  sqrt(float::diag(K))
      Cmax <- max(abs(float::sweep(RHS, 1L, sdx, FUN = "/"))/alpha)
      Cmax <- ifelse(alpha > .Machine$double.eps, Cmax, 5)
      lambda <- exp(seq(log(Cmax),log(lambda.min),length=nlambda))
    }
  }else{
    if(length(dim(lambda))==2){
      if(ncol(lambda) >1 & ncol(lambda)<nTST) stop("Number of columns of 'lambda' must be equal to 'length(tst)'")
    }else{
      lambda <- matrix(lambda, ncol=1)
    }
    if(any(apply(lambda, 2, function(x) any(diff(x) >0)))){
        stop("Object 'lambda' must be a matrix (or vector) of decreasing numbers")
    }
  }

  name <- ifelse(is.null(name),"SSI",name)

  # Split the testing set into subsets. Only the subset provided will be fitted
  if(!is.null(subset)){
     if(!is.numeric(subset) & length(subset) != 2){
       stop("Object 'subset' must contain at least a 2-elements vector")
     }
     sets <- sort(rep(1:subset[2],ceiling(nTST/subset[2]))[1:nTST])
     index <- which(sets == subset[1])
     tmp <- paste0(" of ",length(tst))
     tst <- tst[index]
     RHS <- RHS[,index,drop=FALSE]
     subset_index <- do.call(rbind,lapply(1:subset[2], function(k)range(which(sets == k))))
  }else{
    subset_index <- NULL
    tmp <- ""
  }

  if(verbose){
    message(" Fitting SSI for nTST=",length(tst),tmp," and nTRN=",length(trn)," individuals",sep="")
  }

  save.beta <- ifelse(!is.null(save.at), TRUE, save.beta)

  out <- solveEN(K, RHS, X=yTRN, alpha=alpha, lambda=lambda,
                nlambda=nlambda, lambda.min=lambda.min,
                common.lambda=common.lambda,
                tol=tol, maxiter=maxiter,
                mc.cores=mc.cores, return.beta=return.beta,
                save.beta=save.beta,
                verbose=verbose)

  # If 'save.at' is not NULL
  if(!is.null(save.at)){
    if(!file.exists(dirname(save.at)))  dir.create(dirname(save.at),recursive = TRUE)

    if(!is.null(subset)){
       prefix_file_beta <- paste0(save.at,"subset_",subset[1],"_of_",subset[2],"_beta_i_")
       outfile <- paste0(save.at,"subset_",subset[1],"_of_",subset[2],"_output.RData")
    }else{
       prefix_file_beta <- paste0(save.at,"beta_i_")
       outfile <- paste0(save.at,"output.RData")
    }
    unlink(outfile)
    unlink(paste0(prefix_file_beta,"*.RData"))
    for(k in out$name_beta){
      file.copy(gsub("i_\\*.RData",paste0("i_",k,".RData"),out$file_beta),
                paste0(prefix_file_beta,k,".RData"))
    }
    unlink(out$file_beta)
    out$file_beta <- paste0(prefix_file_beta,"*.RData")
  }

  if(out$q == 1L){
    u0 <- out$yHat
    df0 <- matrix(out$df, nrow=1)
    lambda0 <- matrix(out$lambda, nrow=1)
  }else{
    u0 <- do.call(rbind, out$yHat)
    df0 <- do.call(rbind, out$df)
    lambda0 <- do.call(rbind, out$lambda)
  }
  dimnames(u0) <- list(tst, paste0("SSI.",1:ncol(u0)))

  out <- list(name=name, id=id, y=y, Xb=Xb, u=u0,
              b=b, varU=varU, varE=varE,
              theta=theta, h2=h2, trn=trn, tst=tst, alpha=alpha,
              df = df0, lambda = lambda0,
              beta = out$beta,
              file_beta=out$file_beta,
              name_beta=out$name_beta
            )

  if(!is.null(subset)){
    out$subset <- list(subset=subset, index=subset_index)
  }

  class(out) <- "SSI"

  # Save outputs if 'save.at' is not NULL
  if(!is.null(save.at)){
    save(out, file=outfile)
  }
  return(out)
}
