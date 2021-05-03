
lars2 <- function(P, v, method=c("LAR","LAR-LASSO"), maxDF=NULL,
    eps=.Machine$double.eps,scale=TRUE,verbose=FALSE)
{
  method <- match.arg(method)

  v <- as.vector(v)
  dimnames(P) <- NULL
  p <- length(v)

  if((sum(dim(P))/2)^2 != p^2)
      stop("Object 'P' must be a p*p squared matrix where p=length(v)")

  if(!float::storage.mode(P) %in% c("float32","double")) storage.mode(P) <- "double"
  if(!float::storage.mode(v) %in% c("float32","double")) storage.mode(v) <- "double"
  isFloat <- FALSE
  if(float::storage.mode(P)=="float32" | float::storage.mode(v)=="float32"){
        isFloat <- TRUE
  }

  im <- inactive <- seq(p)
  Sp <- nchar(p)
  textPrint <- c(" Step","\tSec/Step","\tVariable")

  if(scale){
    sdx <- sqrt(float::diag(P))
    P <- cov2cor2(P)  # Equal to P=cov2cor(P) but faster
    v <- v/sdx
  }else{
    sdx <- rep(1,p)
  }

  ignores <- NULL
  if(is.null(maxDF))  maxDF <- p
  beta <- matrix(0,maxDF*8,p)
  lambda <- double(maxDF*8)
  active <- NULL
  drops <- FALSE
  Sign <- NULL
  R <- NULL
  k <- 0
  time <- proc.time()[3]
  while((length(active) < maxDF) & (length(active) < (p-length(ignores))))
  {
    covar <- v[inactive]
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
        R <- upDateR(P[inew,inew],R,as.vector(P[inew,active]),eps=eps)
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
          Sign <- c(Sign,sign(v[inew]))
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
      # a <- drop(w %*% P[active, -c(active,ignores),drop=FALSE])
      a <- float::crossprod(P[active, -c(active,ignores),drop=FALSE],w)[,1]
      gam <- c((Cmax-covar)/(A-a),(Cmax+covar)/(A+a))
      gamhat <- min(gam[gam > eps],Cmax/A)
    }
    if(method == "LAR-LASSO")
    {
      dropid <- NULL
      b1 <- beta[k,active]
      z1 <- -b1/w
      zmin <- min(z1[z1 > eps],gamhat)
      if(zmin < gamhat){
          gamhat <- zmin
          drops <- z1 == zmin
      }else drops <- FALSE
    }
    beta[k+1,] <- beta[k,]
    beta[k+1,active] <- beta[k+1,active] + gamhat*w
    # v <- v - gamhat*P[,active,drop=FALSE]%*%w
    v <- v - gamhat*float::crossprod(P[active, ,drop=FALSE],w)[,1]
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
      beta[k+1,dropid] <- 0
      active <- active[!drops]
      Sign <- Sign[!drops]
    }
    inactive <- im[-c(active, ignores)]
  }
  beta <- beta[seq(k+1), ,drop = FALSE]
  lambda  <-  c(lambda[seq(k)],0)
  if(scale) beta <- sweep(beta,2L,as.numeric(sdx),FUN="/")
  df <- do.call(c,lapply(1:nrow(beta),function(i) sum(abs(beta[i,])>0)))

  out <- list(method=method,beta=beta,lambda=lambda,df=df)
  class(out) <- "LASSO"
  return(out)
}
