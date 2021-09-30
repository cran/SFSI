# Z=U=d=indexK = NULL; BLUP=TRUE; method=c("REML","ML")[1]; h2=NULL
# return.Hinv = FALSE; tol=1E-5; maxIter=1000; interval=c(1E-9,1E9); warn=TRUE

fitBLUP <- function(y, X = NULL, Z = NULL, K = NULL, U = NULL, d = NULL,
                    h2 = NULL, BLUP = TRUE, method = c("REML","ML"),
                    return.Hinv = FALSE, tol = 1E-5, maxIter = 1000,
                    interval = c(1E-9,1E9), warn = TRUE)
{
  method <- match.arg(method)

  if(is.character(K)){
    K <- readBinary(K)
  }
  y <- as.vector(y)
  indexOK <- which(!is.na(y))
  n <- length(indexOK)

  if(is.null(X))
  {
    X <- stats::model.matrix(~1,data=data.frame(rep(1,length(y))))
  }else{
    if(length(dim(X))<2){
      X <- stats::model.matrix(~X)
      if(ncol(X)>2)  colnames(X)[-1] <- substr(colnames(X)[-1],2,nchar(colnames(X)[-1]))
    }
  }
  stopifnot(nrow(X) == length(y))

  isGeigen <- FALSE
  if(is.null(U) & is.null(d))
  {
    if(is.null(Z))
    {
      if(is.null(K)){
        G <- diag(length(y))
      }else  G <- K

    }else{
      if(length(dim(Z)) != 2) stop("Object 'Z' must be a matrix")
      if(is.null(K)){
        G <- float::tcrossprod(Z)  # G = ZKZ'  with K=I
      }else{
        G <- float::tcrossprod(Z,float::tcrossprod(Z,K))  # G = ZKZ'
      }
    }
    stopifnot(nrow(G) == length(y))
    stopifnot(ncol(G) == length(y))

    U <- float::eigen(G[indexOK,indexOK])
    d <- U$values
    U <- U$vectors

  }else{
    isGeigen <- TRUE
    if(is.null(U)) stop("You are providing the eigenvalues, but not the eigenvectors")
    if(is.null(d)) stop("You are providing the eigenvectors, but not the eigenvalues")
    if(n<length(y)) stop("No 'NA' values are allowed when parameters 'U' and 'd' are provided")
  }
  G <- NULL

  tol <- .Machine$double.eps^(3/4)
  if(any(d < tol)){
    if(warn){
      warning("Some eigenvalues are negative or very small:\n  ",sum(d<tol),
      " eigenvalue(s) lie between ",float::dbl(min(d))," and <",tol,".\n",
      "  The corresponding eigenvector(s) will be ignored.",immediate.=TRUE)
    }
    d[d<tol] <- 0
  }

  stopifnot(nrow(U) == n)
  stopifnot(ncol(U) == length(d))

  Uty <- float::crossprod(U,y[indexOK])[,1]
  UtX <- float::crossprod(U,X[indexOK, ,drop=FALSE])

  c0 <- ncol(X)-1

  varP <- var(y[indexOK])*sum(d)/length(d) # mean(d)
  convergence <- lambda0 <- varU <- varE <- bHat <- dbar <- msg <- NA

  if(is.null(h2))
  {
    tt <- searchInt(method,interval=interval,n=n,c0=c0,Uty=Uty,UtX=UtX,d=d,
           maxIter=maxIter,tol=tol,lower=interval[1],upper=interval[2],varP=varP)

    if(is.na(tt$convergence))
    {
      # Divide seeking interval into smaller intervals
      bb <- exp(seq(log(interval[1]),log(interval[2]),length=200))
      tt <- searchInt(method,interval=bb,n=n,c0=c0,Uty=Uty,UtX=UtX,d=d,
             maxIter=maxIter,tol=tol,lower=interval[1],upper=interval[2],varP=varP)

      if(is.na(tt$convergence)){
        # Search in the lower bound
        bb <- exp(seq(log(interval[1]^2),log(interval[1]^0.5),length=200))
        tt <- searchInt(method,interval=bb,n=n,c0=c0,Uty=Uty,UtX=UtX,d=d,
               maxIter=maxIter,tol=tol,lower=interval[1],upper=interval[2],varP=varP)

        if(is.na(tt$convergence)){
          # Search in the upper bound
          bb <- exp(seq(log(interval[2]^0.5),log(interval[2]^2),length=200))
          tt <- searchInt(method,interval=bb,n=n,c0=c0,Uty=Uty,UtX=UtX,d=d,
                 maxIter=maxIter,tol=tol,lower=interval[1],upper=interval[2],varP=varP)
        }
      }
    }
    msg <- tt$msg
    lambda0 <- tt$lambda0
    bHat <- tt$bHat
    dbar <- tt$dbar
    varU <- tt$varU; varE <- tt$varE
    convergence <- ifelse(is.na(tt$convergence),FALSE,tt$convergence)

  }else{
    lambda0 <- h2/(1-h2)
    dbar <- 1/(lambda0*d + 1)
    qq1 <- t(Uty*dbar)%*%UtX
    qq2 <- solve(sweep(t(UtX),2L,dbar,FUN="*")%*%UtX)
    ytPy <- drop(sum(dbar*Uty^2)-qq1%*%qq2%*%t(qq1))
    bHat <- drop(qq2%*%t(qq1))
    varE <- ifelse(method=="REML",ytPy/(n-c0-1),ytPy/n)
    varU <- lambda0*varE
  }


  uHat <- Hinv <- NULL
  if(return.Hinv | (BLUP & !is.na(lambda0) & !isGeigen)){
    if(float::storage.mode(U) == "float32"){
      Hinv <- float::tcrossprod(float::sweep(U,2L,float::fl(lambda0*dbar),FUN="*"),U)
    }else{
      Hinv <- float::tcrossprod(float::sweep(U,2L,lambda0*dbar,FUN="*"),U)
    }
  }

  # Compute BLUP: uHat = KZ'V^{-1} (y-X*b)   with V = varU*ZKZ' + varE*I
  if(BLUP & !is.na(lambda0))
  {
    yStar <- y[indexOK] - X[indexOK ,,drop=FALSE]%*%bHat
    if(isGeigen){
      H <- float::tcrossprod(sweep(U,2L,d*lambda0*dbar,FUN="*"),U)
      uHat <- as.vector(H%*%yStar)
    }else{
      if(is.null(Z) & is.null(K)){  # Z=NULL, K=NULL
        uHat <- rep(0,length(y))
        uHat[indexOK] <- as.vector(Hinv%*%yStar)   # V^{-1}*(y-Xb)

      }else{
        if(is.null(Z)){     # Z=NULL, K=K
          uHat <- as.vector(float::crossprod(K[indexOK,,drop=FALSE],Hinv)%*%yStar)  # K[,trn]*V^{-1}*(y-Xb)
        }else{
          if(is.null(K)){   # Z=Z, K=NULL
              uHat <- as.vector(float::crossprod(Z[indexOK,,drop=FALSE],Hinv)%*%yStar)  # Z[,trn]'*V^{-1}*(y-Xb)
          }else{            # Z=Z, K=K
              ZKt <- float::tcrossprod(Z[indexOK, ,drop=FALSE],K)   # ZK' which is the transpose of KZ'
              uHat <- as.vector(float::crossprod(ZKt,Hinv)%*%yStar)
          }
        }
      }
    }

    if(!return.Hinv) Hinv <- NULL
  }

  if(warn){
    if(ifelse(is.na(convergence),is.na(lambda0),!convergence)){
      warning("Convergence was not reached in the 'GEMMA' algorithm.",immediate.=TRUE)
    }
    if(!is.na(msg)) warning(msg,immediate.=TRUE)
  }

  h2 <- varU/(varU + varE)

  out <- list(varE = varE, varU = varU, h2 = h2, b = bHat, u = uHat,
              Hinv = Hinv, convergence = convergence, method = method)
  return(out)
}
