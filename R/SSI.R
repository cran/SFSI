
# X = Z = saveAt = name = subset= NULL; lambda=b=NULL
# alpha = 1; nLambda = 100; commonLambda = TRUE; minLambda = .Machine$double.eps^0.5
# mc.cores = 1; tol = 1E-4; maxIter = 500; verbose = TRUE; method = c("REML","ML")[1]

SSI <- function(y, X = NULL, b = NULL, Z = NULL, K, D = NULL,
         h2 = NULL, trn = seq_along(y), tst = seq_along(y),
         subset = NULL, alpha = 1, lambda = NULL, nLambda = 100,
         minLambda = .Machine$double.eps^0.5, commonLambda = TRUE,
         tol = 1E-4, maxIter = 500, method = c("REML","ML"),
         saveAt = NULL, name = NULL, mc.cores = 1, verbose = TRUE)
{
  method <- match.arg(method)

  n <- length(y)
  if(is.logical(trn)){
     if(n != length(trn)) stop("Object 'trn' must be of the same length of 'y'\n")
     trn <- which(trn)
  }
  if(is.logical(tst)){
     if(n != length(tst)) stop("Object 'tst' must be of the same length of 'y'\n")
     tst <- which(tst)
  }
  nTRN <- length(trn);  nTST <- length(tst)

  if(is.character(K)){
      K <- readBinary(K)
  }

  if(!float::storage.mode(K) %in% c("float32","double")) storage.mode(K) <- "double"
  isFloat <- float::storage.mode(K)=="float32"

  if(is.null(X))   # Design matrix for fixed effects including the intercept
  {
    X <- model.matrix(~1,data=data.frame(rep(1,n)))
  }else{
    if(length(dim(X))<2){
      X <- stats::model.matrix(~X)
      if(ncol(X)>2)  colnames(X)[-1] <- substr(colnames(X)[-1],2,nchar(colnames(X)[-1]))
    }
  }

  if(!is.null(Z)){
    if(length(dim(Z))!=2) stop("Object 'Z' must be a matrix with nrow(Z)=n and ncol(Z)=nrow(K)\n")
    K <- float::tcrossprod(Z,float::tcrossprod(Z,K))   # Z%*%K%*%t(Z)
  }

  if(length(dim(K))!=2 | (length(K) != n^2)){
    stop("Product Z %*% K %*% t(Z) must be a squared matrix with number of rows",
     "\n(and columns) equal to n: number of elements in 'y'")
  }
  dimnames(K) <- NULL

  if(!is.null(D)){
    if((sum(dim(D))/2)^2 != n^2){
      stop("Object 'D' must be a nxn squared matrix with n: number of elements in 'y'")
    }
    dimnames(D) <- NULL
    D <- D[trn,trn]
  }

  RHS <- K[trn,tst,drop=FALSE]
  K <- K[trn,trn]

  if(is.null(h2))
  {
    # Fit LMM to get variance components and estimate fixed effects as GLS
    fm <- fitBLUP(y[trn],X=X[trn, ,drop=FALSE],K=K,method=method)
    if(fm$convergence){
      varU <- fm$varU
      varE <- fm$varE
      h2 <- varU/(varU + varE)
      b <- fm$b
    }else stop("Convergence was not reached in the 'GEMMA' algorithm. \n\t",
           "Please provide a heritability estimate in 'h2' parameter")
  }else{   # Only estimate fixed effects as GLS
    varU <- varE <- NULL
    if(is.null(b)){
      b <- fitBLUP(y[trn],X=X[trn,,drop=FALSE],K=K,BLUP=FALSE,h2=h2)$b
    }else{
      if(length(b) != ncol(X)) stop("The length of 'b' must be the same as the number of columns of 'X'\n")
    }
  }
  if(h2 < 0.001) warning("The 'heritability' is too small. Results may be affected",immediate.=TRUE)
  Xb <- drop(X%*%b)

  # Standardizing
  if(is.null(D)){
    add2diag(K,a=(1-h2)/h2,void=TRUE)   #diag(K) <- diag(K) + (1-h2)/h2
  }else{
    K <- K + ((1-h2)/h2)*D
  }
  sdx <- sqrt(float::diag(K))
  cov2cor2(K,void=TRUE)   # Equal to K=cov2cor(K) but faster
  RHS <- float::sweep(RHS,1L,sdx,FUN="/")  # Scale each row of RHS

  if(is.null(lambda)){
    if(commonLambda){
        Cmax <- ifelse(alpha > .Machine$double.eps,max(abs(RHS)/alpha),5)
        lambda <- exp(seq(log(Cmax),log(minLambda),length=nLambda))
    }
  }else{
    if(length(dim(lambda))==2){
        if(nrow(lambda) != nTST) stop("Object 'lambda' must be a vector or a matrix with nTST rows")
    }else{
      if(length(dim(lambda))==2 | mode(lambda)!="numeric" | any(diff(lambda) >0))
        stop("Object 'lambda' must be a vector of decreasing numbers")
    }
  }

  name <- ifelse(is.null(name),"SSI",name)

  # Split the testing set into subsets. Only the subset provided will be fitted
  if(!is.null(subset)){
     if(!is.numeric(subset) & length(subset) != 2)
      stop("Object 'subset' must contain at least a 2-elements vector")
     sets <- sort(rep(1:subset[2],ceiling(nTST/subset[2]))[1:nTST])
     index <- which(sets == subset[1])
     tmp <- paste0(" of ",length(tst))
     tst <- tst[index]
     RHS <- RHS[,index,drop=FALSE]
     subset_size <- as.vector(table(sets)[as.character(1:subset[2])])
  }else{
    subset_size <- NULL
    tmp <- ""
  }

  if(verbose){
    cat(" Fitting SSI model for nTST=",length(tst),tmp," and nTRN=",nTRN," individuals\n",sep="")
  }

  compApply <- function(ind)
  {
    rhs <- drop(RHS[,ind])
    if(length(dim(lambda))==2){
      lambda0 <- lambda[ind,]
    }else lambda0 <- lambda

    fm <- solveEN(K,rhs,scale=FALSE,lambda=lambda0,nLambda=nLambda,
               minLambda=minLambda,alpha=alpha,tol=tol,maxIter=maxIter)

    # Return betas to the original scale by dividing by sdx
    saveBinary(sweep(fm$beta,1L,as.numeric(sdx),FUN="/"),file=paste0(tmpdir,"/",file_beta,ind,".bin"),
         type=file_type,verbose=FALSE)

    if(verbose){
      cat(1,file=con,append=TRUE)
      utils::setTxtProgressBar(pb, nchar(scan(con,what="character",quiet=TRUE))/length(tst))
    }

    return(list(lambda=fm$lambda,tst=tst[ind],df=fm$df))
  }

  tmpdir <- tempdir()
  file_beta <- paste0(tempfile(tmpdir=""),"_beta_ind_")
  file_type <- ifelse(isFloat,"float","double")
  if(verbose){
     pb = utils::txtProgressBar(style=3)
     con <- tempfile(tmpdir=tmpdir)
  }
  if(mc.cores == 1L) {
    out = lapply(X=seq_along(tst),FUN=compApply)
  }else{
    out = parallel::mclapply(X=seq_along(tst),FUN=compApply,mc.cores=mc.cores)
  }
  if(verbose) {
    close(pb); unlink(con)
  }

  if(sum(tst != unlist(lapply(out,function(x)x$tst)))>0){
    stop("Some sub-processes failed. Something went wrong during the analysis.")
  }

  out <- list(name=name, y=y, Xb=Xb, b=b, varU=varU,
              varE=varE, h2=h2, trn=trn, tst=tst, alpha=alpha,
              df = do.call("rbind",lapply(out,function(x)x$df)),
              lambda = do.call("rbind",lapply(out,function(x)x$lambda)),
              file_beta=paste0(tmpdir,"/",file_beta))
  class(out) <- "SSI"

  # Save outputs if 'saveAt' is not NULL
  if(!is.null(saveAt)){
    if(!is.null(subset)){
       filenames <- paste0(saveAt,c("_B_","_"),subset[1],"_of_",subset[2],c(".bin",".RData"))
    }else  filenames <- paste0(saveAt,c("_B.bin",".RData"))

    if(!file.exists(dirname(filenames[1])))  dir.create(dirname(filenames[1]),recursive = TRUE)

    saveBinary(do.call(cbind,coef.SSI(out)),file=filenames[1],type="float",verbose=FALSE)
    out$file_beta <- filenames[1]
    out$subset <- subset
    out$subset_size <- subset_size
    save(out,file=filenames[2])
  }
  return(out)
}
