
# X <- Z <- K <- U <- d <- NULL; scale = TRUE; pairwise=FALSE; verbose = TRUE
getGenCov <- function(y, X = NULL, Z = NULL, K = NULL,
                      U = NULL, d = NULL, scale = TRUE,
                      pairwise=FALSE, verbose = TRUE, ...)
{
  # K=G0; X = NULL; Z = NULL; U = NULL; d = NULL; scale = TRUE
  eps <- .Machine$double.eps

  if(length(dim(y)) == 2L){
    y <- as.matrix(y)
  }else{
    y <- matrix(y, ncol=1L)
  }

  trn <- which(apply(y,1,function(x)all(!is.na(x))))  # Common TRN set
  if(length(trn) == 0){
    stop("No common training set to all response variables was found")
  }else{
    if(verbose){
      message(" Calculating genetic covariances using n=",length(trn)," training observations")
    }
    y <- y[trn,,drop=FALSE]
  }

  n <- nrow(y)
  p <- ncol(y)

  if(scale){
    sdy <- as.vector(apply(y, 2, sd, na.rm=TRUE))
    y <- scale(y, center=FALSE, scale=sdy)
  }else{
    sdy <- rep(1,p)
  }

  # Create an index for pairwise models
  tmp <- expand.grid(j=1:p,i=1:p)
  INDEX <- data.frame(pos=p+seq(p*(p-1)/2),
                      tmp[tmp$i < tmp$j,c("i","j")])
  if(!pairwise){
    INDEX <- INDEX[INDEX$i==1,]
  }

  Y0 <- cbind(y, matrix(NA, nrow=n, ncol=nrow(INDEX)))

  for(k in 1:nrow(INDEX)){
    Y0[,INDEX$pos[k]] <- y[,INDEX$i[k]] + y[,INDEX$j[k]]
  }

  fm <- fitBLUP(Y0, BLUP=FALSE, X=X, Z=Z, K=K, U=U, d=d, verbose=verbose, ...)

  varUi <- fm$varU[1:p]*(sdy^2)   # Scale using their initial SD
  varEi <- fm$varE[1:p]*(sdy^2)

  # Fixed effects
  if(is.null(fm$b)){
    b <- NULL
  }else{
    b <- sweep(fm$b[,1:p,drop=F],2L,sdy,FUN="*")
  }

  if(pairwise){
    varU <- varE <- matrix(NA, ncol=p, nrow=p)
    diag(varU) <- varUi
    diag(varE) <- varEi

    if(!is.null(colnames(y))){
      dimnames(varU) <- dimnames(varE) <- list(colnames(y),colnames(y))
    }

    out <- list(varU=varU, varE=varE, b=b)
  }else{
    # Only the gencov between the first y[,1] and the
    # remaining y[,2:p] response variables
    varU <- varUi
    varE <- varEi
    covU <- covE <- rep(NA, p-1)

    if(!is.null(colnames(y))){
      names(varU) <- names(varE) <- colnames(y)
      names(covU) <- names(covE) <- colnames(y)[-1]
    }

    out <- list(varU=varU, varE=varE, covU=covU, covE=covE, b=b)
  }

  for(k in 1:nrow(INDEX)){
    pos <- INDEX$pos[k]
    i <- INDEX$i[k]
    j <- INDEX$j[k]

    # Genetic and Environmental covariances
    varUij <- 0.5*sdy[i]*sdy[j]*(fm$varU[pos] - fm$varU[i] - fm$varU[j])
    varEij <- 0.5*sdy[i]*sdy[j]*(fm$varE[pos] - fm$varE[i] - fm$varE[j])

    if(pairwise){
      out$varU[i,j] <- varUij
      out$varE[i,j] <- varEij

      out$varU[j,i] <- out$varU[i,j]
      out$varE[j,i] <- out$varE[i,j]
    }else{
      out$covU[k] <- varUij
      out$covE[k] <- varEij
    }
  }

  return(out)
}
