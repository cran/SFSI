
# y1=y; y2=x; X=Z=U=d=NULL; K=G; mc.cores=5; scale=TRUE
getGenCov <- function(y1, y2, X = NULL, Z = NULL, K = NULL, U = NULL,
           d = NULL, scale = TRUE, mc.cores = 1, warn = FALSE, ...)
{
  if(!is.matrix(y2))
    stop("Object 'y2' must be a matrix with 'nrow(y2)' equal to the number of elements in 'y1'")

  sdy1 <- sd(y1)
  sdy2 <- apply(y2,2,sd)
  if(scale){
    y1 <- as.vector(y1/sdy1)
    y2 <- scale(y2,FALSE,sdy2)
    #y1 <- as.vector(scale(y1))
    #y2 <- scale(y2)
  }else{
    if(any(abs(sdy1 -sdy2) > sqrt(.Machine$double.eps)))
      warning("Variances of y1 and y2 are not equal",immediate.=TRUE)
  }

  if(length(y1) != nrow(y2))
    stop("The number of elements in 'y1' must be equal to the number of rows in 'y2'")

  if(is.null(U) & is.null(d))
  {
    if(is.null(Z))
    {
      if(is.null(K)){
          K <- diag(length(y1))
      }
      G <- K
    }else{
      if(length(dim(Z)) != 2) stop("Object 'Z' must be a matrix")
      if(is.null(K)){
        G <- float::tcrossprod(Z)  # G = ZKZ'  with K=I
      }else{
        G <- float::tcrossprod(Z,float::tcrossprod(Z,K))  # G = ZKZ'
      }
    }

    stopifnot(nrow(G) == length(y1))
    stopifnot(ncol(G) == length(y1))
    out <- float::eigen(G)
    d <- out$values
    U <- out$vectors
  }else{
    if(is.null(U)) stop("You are providing the eigenvalues, but not the eigenvectors")
    if(is.null(d)) stop("You are providing the eigenvectors, but not the eigenvalues")
  }

  fm1 <- fitBLUP(y1,BLUP=FALSE,X=X,U=U,d=d,warn=warn, ...)   # Model for variable 1

  compApply <- function(j)
  {
     fm2 <- fitBLUP(y2[,j],BLUP=FALSE,X=X,U=U,d=d,warn=warn, ...)       # Model for variable 2
     fm3 <- fitBLUP(y1 + y2[,j],BLUP=FALSE,X=X,U=U,d=d,warn=warn, ...)  # Model for variable 3

     cat(1, file = con, append = TRUE)
     utils::setTxtProgressBar(pb,nchar(scan(con,what="character",quiet=TRUE))/ncol(y2))

     c(varU2=fm2$varU, covU=0.5*(fm3$varU - fm1$varU - fm2$varU),
        	varE2=fm2$varE, covE=0.5*(fm3$varE - fm1$varE - fm2$varE))
  }

  pb = utils::txtProgressBar(style = 3); con <- tempfile()
  if(mc.cores == 1L){
    out = lapply(X=1:ncol(y2), FUN=compApply)
  }else{
    out = parallel::mclapply(X=1:ncol(y2), FUN=compApply, mc.cores=mc.cores)
  }
  close(pb)
  unlink(con)

  out <- data.frame(do.call(rbind,out))
  list(varU1=fm1$varU, varE1=fm1$varE, varU2=out$varU2,
      varE2=out$varE2, covU=out$covU,covE=out$covE)
}
