
LARS <- function(Sigma, Gamma, X = NULL, method=c("LAR","LAR-LASSO"),
                 dfmax = NULL, eps = .Machine$double.eps,
                 scale = TRUE, mc.cores = 1L, return.beta = TRUE,
                 save.beta = FALSE, verbose = FALSE)
{
  method <- match.arg(method)

  if(length(dim(Gamma)) != 2){
      Gamma <- matrix(Gamma, ncol=1)
  }
  dimnames(Sigma) <- NULL
  p <- nrow(Gamma)
  q <- ncol(Gamma)

  if((sum(dim(Sigma))/2)^2 != p^2){
      stop("Object 'Sigma' must be a p*p squared matrix where p=nrow(Gamma)")
  }

  if(!float::storage.mode(Sigma) %in% c("float32","double")) storage.mode(Sigma) <- "double"
  if(!float::storage.mode(Gamma) %in% c("float32","double")) storage.mode(Gamma) <- "double"
  isfloat <- FALSE
  if(float::storage.mode(Sigma)=="float32" | float::storage.mode(Gamma)=="float32"){
        isfloat <- TRUE
  }

  if(!is.null(X)){
    if(length(dim(X)) != 2L){
      X <- matrix(X, nrow=1L)
    }
    stopifnot(ncol(X) == nrow(Gamma))
  }

  Sp <- nchar(p)
  textPrint <- c(" Step","\tSec/Step","\tVariable")

  if(scale){
    sdx <- float::dbl(sqrt(float::diag(Sigma)))
    Sigma <- cov2cor2(Sigma)  # Equal to Sigma=cov2cor(Sigma) but faster
    Gamma <- float::sweep(Gamma, 1L, sdx, FUN = "/")
  }else{
    sdx <- rep(1,p)
  }

  if(return.beta & save.beta){
    message("'return.beta' is set to FALSE when 'save.beta=TRUE'")
    return.beta <- FALSE
  }

  compApply <- function(ind)
  {
    rhs <- as.vector(Gamma[,ind])

    im <- inactive <- seq(p)
    ignores <- NULL
    if(is.null(dfmax))  dfmax <- p
    beta <- Matrix::Matrix(0, nrow=p, ncol=dfmax*8)
    lambda <- double(dfmax*8)
    active <- NULL
    drops <- FALSE
    Sign <- NULL
    R <- NULL
    k <- 0
    time <- proc.time()[3]
    while((length(active) < dfmax) & (length(active) < (p-length(ignores))))
    {
      covar <- rhs[inactive]
      Cmax <- max(abs(covar))
      if(Cmax < eps*100){
        if(verbose & q==1L) cat(" Max |corr| = 0; exiting...\n")
        break
      }
      k <- k+1
      lambda[k] <- ifelse(isfloat, float::dbl(Cmax), Cmax)
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
            if(verbose & q==1L){
              cat(" LARS Step ",k,":\t Variable ", inew,"\tcollinear; dropped for good\n",sep="")
            }
          }else{
            active <- c(active,inew)
            Sign <- c(Sign,sign(rhs[inew]))
            if(verbose & q==1L){
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
      # rhs <- rhs - gamhat*Sigma[,active,drop=FALSE]%*%w
      rhs <- rhs - gamhat*float::crossprod(Sigma[active, ,drop=FALSE],w)[,1]
      if(method == "LAR-LASSO" && any(drops))
      {
        dropid <- seq(drops)[drops]
        for(id in rev(dropid))
        {
          if(verbose & q==1L){
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
    if(scale){
       beta <- sweep(beta, 1L, sdx, FUN="/")
    }
    df <- do.call(c,lapply(1:ncol(beta),function(i) sum(abs(beta[,i])>0)))

    if(!is.null(X)){
      yHat <- X%*%as.matrix(beta)
    }else{
      yHat <- NULL
    }

    if(!return.beta){
      if(save.beta){
        save(beta, file=paste0(prefix_file_beta,ind,".RData"))
      }
      beta <- NULL
    }

    if(verbose & q>1L){
      cat(1,file=con,append=TRUE)
      utils::setTxtProgressBar(pb, nchar(scan(con,what="character",quiet=TRUE))/q)
    }

    list(i=ind, yHat=yHat, beta=beta, lambda=lambda, df=df)
  }

  tmpdir0 <- tempdir()
  prefix_file_beta <- paste0(tempfile(tmpdir=tmpdir0),"_beta_i_")
  unlink(paste0(prefix_file_beta,"*.RData"))

  if(verbose & q>1L){
     pb = utils::txtProgressBar(style=3)
     con <- tempfile(tmpdir=tmpdir0)
  }
  if(mc.cores == 1L){
    out = lapply(X=seq(q),FUN=compApply)
  }else{
    out = parallel::mclapply(X=seq(q),FUN=compApply,mc.cores=mc.cores)
  }
  if(verbose & q>1L) {
    close(pb); unlink(con)
  }

  # Checkpoint
  if(any(seq(q) != unlist(lapply(out,function(x) x$i)) )){
      stop("Some sub-processes failed. Something went wrong during the analysis.")
  }

  out <- list(p=p, q=q, method=method,
              yHat = lapply(out, function(x)x$yHat),
              lambda = lapply(out, function(x)x$lambda),
              df = lapply(out, function(x)x$df),
              beta = lapply(out, function(x)x$beta),
              file_beta = paste0(prefix_file_beta,"*.RData"),
              name_beta = seq(q)
            )
  if(q == 1L){
    out$yHat <- out$yHat[[1]]
    out$df <- out$df[[1]]
    out$lambda <- out$lambda[[1]]
    out$beta <- out$beta[[1]]
  }

  if(is.null(X)){
    out$yHat <- NULL
  }

  if(save.beta | !return.beta){
    out$beta <- NULL
  }
  if(!save.beta){
    out$file_beta <- out$name_beta <- NULL
  }

  class(out) <- "LASSO"
  return(out)
}
