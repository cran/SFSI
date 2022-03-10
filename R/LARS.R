
LARS <- function(Sigma, Gamma, method=c("LAR","LAR-LASSO"), dfmax = NULL,
    eps = .Machine$double.eps, scale = TRUE, verbose = FALSE)
{
  method <- match.arg(method)

  Gamma <- as.vector(Gamma)
  dimnames(Sigma) <- NULL
  p <- length(Gamma)

  if((sum(dim(Sigma))/2)^2 != p^2)
      stop("Object 'Sigma' must be a p*p squared matrix where p=length(Gamma)")

  if(!float::storage.mode(Sigma) %in% c("float32","double")) storage.mode(Sigma) <- "double"
  if(!float::storage.mode(Gamma) %in% c("float32","double")) storage.mode(Gamma) <- "double"
  isFloat <- FALSE
  if(float::storage.mode(Sigma)=="float32" | float::storage.mode(Gamma)=="float32"){
        isFloat <- TRUE
  }

  im <- inactive <- seq(p)
  Sp <- nchar(p)
  textPrint <- c(" Step","\tSec/Step","\tVariable")

  if(scale){
    sdx <- sqrt(float::diag(Sigma))
    Sigma <- cov2cor2(Sigma)  # Equal to Sigma=cov2cor(Sigma) but faster
    Gamma <- Gamma/sdx
  }else{
    sdx <- rep(1,p)
  }

  ignores <- NULL
  if(is.null(dfmax))  dfmax <- p
  beta <- matrix(0, p, dfmax*8)  # p x nLambda
  lambda <- double(dfmax*8)
  active <- NULL
  drops <- FALSE
  Sign <- NULL
  R <- NULL
  k <- 0
  time <- proc.time()[3]
  while((length(active) < dfmax) & (length(active) < (p-length(ignores))))
  {
    covar <- Gamma[inactive]
    Cmax <- max(abs(covar))
    if(Cmax < eps*100){
      if(verbose) cat(" Max |corr| = 0; exiting...\n")
      break
    }
    k <- k+1
    lambda[k] <- ifelse(isFloat,float::dbl(Cmax),Cmax)
    if(!any(drops))
    {
      new <- abs(covar) >= Cmax-eps
      covar <- covar[!new]
      new <- inactive[new]
      for(inew in new)
      {
        R <- upDateR(Sigma[inew,inew],R,as.vector(Sigma[inew,active]),eps=eps)
        if(attr(R,"rank")==length(active))
        {
          nR <- seq(length(active))
          R <- R[nR,nR,drop=FALSE]
          attr(R,"rank") <- length(active)
          ignores <- c(ignores,inew)
          if(verbose){
            cat(" LARS Step ",k,":\t Variable ", inew,"\tcollinear; dropped for good\n",sep="")
          }
        }else{
          active <- c(active,inew)
          Sign <- c(Sign,sign(Gamma[inew]))
          if(verbose){
            cat("--------------------------------------------------------------\n")
            tmp <- proc.time()[3]
            cat(paste(textPrint,"=",c(sprintf("%*d",Sp,k),sprintf('%.*f',4,tmp-time),
                sprintf("%*d",Sp,inew)))," added\n",sep="")
            time <- tmp
          }
        }
      }
    }
    Gi1 <- float::backsolve(R,backsolvet(R,Sign))
    A <- 1/sqrt(sum(Gi1*Sign))
    w <- A*Gi1

    if((length(active) >= (p-length(ignores)))){
        gamhat <- Cmax/A
    }else{
      # a <- drop(w %*% Sigma[active, -c(active,ignores),drop=FALSE])
      a <- float::crossprod(Sigma[active, -c(active,ignores),drop=FALSE],w)[,1]
      gam <- c((Cmax-covar)/(A-a),(Cmax+covar)/(A+a))
      gamhat <- min(gam[gam > eps],Cmax/A)
    }
    if(method == "LAR-LASSO")
    {
      dropid <- NULL
      b1 <- beta[active,k]
      z1 <- -b1/w
      zmin <- min(z1[z1 > eps],gamhat)
      if(zmin < gamhat){
          gamhat <- zmin
          drops <- z1 == zmin
      }else drops <- FALSE
    }
    beta[,k+1] <- beta[,k]
    beta[active,k+1] <- beta[active,k+1] + gamhat*w
    # Gamma <- Gamma - gamhat*Sigma[,active,drop=FALSE]%*%w
    Gamma <- Gamma - gamhat*float::crossprod(Sigma[active, ,drop=FALSE],w)[,1]
    if(method == "LAR-LASSO" && any(drops))
    {
      dropid <- seq(drops)[drops]
      for(id in rev(dropid))
      {
        if(verbose){
          cat("--------------------------------------------------------------\n")
          tmp <- proc.time()[3]
          cat(paste(textPrint,"=",c(sprintf("%*d",Sp,k+1),sprintf('%.*f',4,tmp-time),
            sprintf("%*d",Sp,active[id])))," dropped\n",sep="")
          time <- tmp
        }
        R <- downDateR(R,id)
      }
      dropid <- active[drops]
      beta[dropid,k+1] <- 0
      active <- active[!drops]
      Sign <- Sign[!drops]
    }
    inactive <- im[-c(active, ignores)]
  }
  beta <- beta[,seq(k+1) ,drop = FALSE]
  lambda  <-  c(lambda[seq(k)],0)
  if(scale) beta <- sweep(beta,1L,as.numeric(sdx),FUN="/")
  df <- do.call(c,lapply(1:ncol(beta),function(j) sum(abs(beta[,j])>0)))

  out <- list(method=method,beta=beta,lambda=lambda,df=df)
  class(out) <- "LASSO"
  return(out)
}
