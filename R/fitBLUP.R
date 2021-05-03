
# X= K= U = d = indexK = NULL; BLUP=TRUE; method="ML"; h2=NULL
# return.Hinv = FALSE; tol=1E-5; maxIter=1000; interval=c(1E-9,1E9)

fitBLUP <- function(y, X = NULL, Z = NULL, K = NULL, U = NULL, d = NULL,
                    indexK = NULL, h2 = NULL, BLUP = TRUE, method = c("REML","ML"),
                    return.Hinv = FALSE, tol = 1E-5, maxIter = 1000,
                    interval = c(1E-9,1E9), warn = TRUE)
{
  method <- match.arg(method)

  if(is.character(K)){
    K <- readBinary(K,indexRow=indexK,indexCol=indexK)
  }
  y <- as.vector(y)
  indexOK <- which(!is.na(y))
  n <- length(indexOK)

  if(is.null(X))
  {
    X <- stats::model.matrix(~1,data=data.frame(rep(1,length(y))))
  }else{
    if(is.vector(X)){
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

    out <- float::eigen(G[indexOK,indexOK])
    d <- out$values
    U <- out$vectors

  }else{
    isGeigen <- TRUE
    if(is.null(U)) stop("You are providing the eigenvalues, but not the eigenvectors")
    if(is.null(d)) stop("You are providing the eigenvectors, but not the eigenvalues")
    if(n<length(y)) stop("No 'NA' values are allowed when parameters 'U' and 'd' are provided")
  }

  if(min(d) <0 & abs(min(d)) > .Machine$double.eps^0.5) stop("Matrix 'K' is not positive semi-definite")
  stopifnot(nrow(U) == n)
  stopifnot(ncol(U) == length(d))

  Uty <- float::crossprod(U,y[indexOK])[,1]
  UtX <- float::crossprod(U,X[indexOK, ,drop=FALSE])
  c0 <- ncol(X)-1

  varP <- var(y[indexOK])*sum(d)/length(d) #mean(d)
  convergence <- lambda0 <- varU <- varE <- bHat <- dbar <- NA
  if(is.null(h2))
  {
    tt <- searchInt(method,interval,n=n,c0=c0,Uty=Uty,UtX=UtX,d=d,maxIter=maxIter,
       tol=tol,lower=interval[1],upper=interval[2],varP=varP)
    convergence <- tt$convergence

    if(is.na(convergence))
    {
      # Divide seeking interval into smaller intervals
      bb <- exp(seq(log(interval[1]),log(interval[2]),length=200))
      tt <- searchInt(method,bb,n=n,c0=c0,Uty=Uty,UtX=UtX,d=d,maxIter=maxIter,
        tol=tol,lower=interval[1],upper=interval[2],varP=varP)
      convergence <- tt$convergence

      if(is.na(convergence)){
        # Search in the lower bound
        bb <- exp(seq(log(interval[1]^2),log(interval[1]),length=200))
        tt <- searchInt(method,bb,n=n,c0=c0,Uty=Uty,UtX=UtX,d=d,maxIter=maxIter,
          tol=tol,lower=interval[1],upper=interval[2],varP=varP)
        convergence <- tt$convergence

        if(is.na(convergence)){
          # Search in the upper bound
          bb <- exp(seq(log(interval[2]),log(interval[2]^2),length=200))
          tt <- searchInt(method,bb,n=n,c0=c0,Uty=Uty,UtX=UtX,d=d,maxIter=maxIter,
            tol=tol,lower=interval[1],upper=interval[2],varP=varP)
          convergence <- tt$convergence
        }
      }
    }
    lambda0 <- tt$lambda0
    bHat <- tt$bHat
    dbar <- tt$dbar
    varU <- tt$varU; varE <- tt$varE

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

  # Compute BLUP: uHat = KZ'V^{-1} (y-X*b)   with V = varU*ZKZ' + varE*I
  uHat <- Hinv <- NULL
  if(BLUP & !is.na(lambda0))
  {
    yStar <- y[indexOK] - X[indexOK ,,drop=FALSE]%*%bHat
    if(isGeigen){
      H <- tcrossprod(sweep(U,2L,d*lambda0*dbar,FUN="*"),U)
      uHat <- drop(H%*%yStar)
    }else{
      Hinv <- float::tcrossprod(float::sweep(U,2L,lambda0*dbar,FUN="*"),U) # Vinv
      if(is.null(Z) & is.null(K)){  # Z=NULL, K=NULL
        uHat <- rep(0,length(y))
        uHat[indexOK] <- drop(Hinv%*%yStar)   # V^{-1}*(y-Xb)

      }else{
        if(is.null(Z)){     # Z=NULL, K=K
          uHat <- drop(float::crossprod(K[indexOK,,drop=FALSE],Hinv)%*%yStar)  # K[,trn]*V^{-1}*(y-Xb)
        }else{
          if(is.null(K)){   # Z=Z, K=NULL
              uHat <- drop(float::crossprod(Z[indexOK,,drop=FALSE],Hinv)%*%yStar)  # Z[,trn]'*V^{-1}*(y-Xb)
          }else{            # Z=Z, K=K
              KZt <- float::tcrossprod(K,Z[indexOK,,drop=FALSE])
              uHat <- drop(KZt%*%Hinv%*%yStar)
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
  }

  h2 <- varU/(varU + varE)

  out <- list(varE = varE, varU = varU, h2 = h2, b = bHat, u = uHat,
              Hinv = Hinv, convergence = convergence, method = method)
  return(out)
}
