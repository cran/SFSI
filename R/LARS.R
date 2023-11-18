
LARS <- function(Sigma, Gamma, method=c("LAR","LASSO"),
                 nsup.max = NULL, eps = .Machine$double.eps*100,
                 scale = TRUE, sdx = NULL, mc.cores = 1L, save.at = NULL,
                 precision.format = c("double","single"),
                 fileID = NULL, verbose = FALSE)
{
  precision.format <- match.arg(precision.format)
  method <- match.arg(method)

  if(length(dim(Gamma)) == 2L){
    Gamma <- as.matrix(Gamma)
  }else{
    Gamma <- matrix(Gamma, ncol=1L)
  }
  p <- nrow(Gamma)
  q <- ncol(Gamma)

  if((sum(dim(Sigma))/2)^2 != p^2){
    stop("Input 'Sigma' must be a p*p squared matrix where p=nrow(Gamma)")
  }

  scaleb <- TRUE
  if(scale){
    sdx <-  sqrt(diag(Sigma))
    cov2cor2(Sigma, inplace=TRUE)     # Equal to Sigma=cov2cor(Sigma) but faster
    Gamma <- sweep(Gamma, 1L, sdx, FUN="/")
  }else{
    if(is.null(sdx)){
      scaleb <- FALSE
    }else{
      if(length(sdx) != p){
        stop("Input 'sdx' must be a numeric vector of length = ",p)
      }
    }
  }

  nsup.max <- ifelse(is.null(nsup.max), p, nsup.max)
  flagsave <- as.logical(!is.null(save.at))
  isLASSO <- as.logical(method=="LASSO")
  verbose2 <- as.logical(verbose & q==1L)
  mc.cores <- ifelse(q==1L & mc.cores>1L, 1L, mc.cores)
  doubleprecision <- as.logical(precision.format=="double")

  compApply <- function(ind)
  {
    rhs <- as.vector(Gamma[,ind])

    if(flagsave){
      filename <- paste0(file_beta,fileID[ind],".bin")
    }else{
      filename <- NULL
    }

    #dyn.load("c_lasso.so")
    res <- .Call("R_lars", Sigma, rhs, eps, nsup.max,
                 scaleb, sdx, isLASSO, filename,
                 doubleprecision, verbose2)
    #dyn.unload("c_lasso.so")

    if(verbose & q>1L){
      cat(1,file=con,append=TRUE)
      utils::setTxtProgressBar(pb, nchar(scan(con,what="character",quiet=TRUE))/q)
    }

    list(ind=ind, beta=res[[1]], lambda=res[[2]], nsup=res[[3]])
  }

  tmpdir0 <- tempdir()

  file_beta <- NULL
  if(flagsave){
    stopifnot(is.character(save.at))
    save.at <- normalizePath(save.at, mustWork=F)
    file_beta <- paste0(save.at,"beta_i_")

    if(!file.exists(dirname(file_beta))){
      dir.create(dirname(file_beta),recursive=TRUE)
    }

    if(is.null(fileID)){
      fileID <- seq(q)
    }else{
      stopifnot(length(fileID)==q)
    }
  }

  # Run the analysis for 1:ncol(Gamma)
  if(verbose & q>1L){
     pb <- utils::txtProgressBar(style=3)
     con <- tempfile(tmpdir=tmpdir0)
  }
  if(mc.cores == 1L){
    out <- lapply(X=seq(q),FUN=compApply)
  }else{
    out <- parallel::mclapply(X=seq(q),FUN=compApply,mc.cores=mc.cores)
  }
  if(verbose & q>1L) {
    close(pb); unlink(con)
  }

  # Checkpoint
  if(any(seq(q) != unlist(lapply(out,function(x) x$ind)) )){
    stop("Some sub-processes failed. Something went wrong during the analysis.")
  }

  out <- list(p=p, q=q, method=method,
              nlambda=unlist(lapply(out, function(x)length(x$lambda))),
              lambda = lapply(out, function(x)x$lambda),
              nsup = lapply(out, function(x)x$nsup),
              beta = lapply(out, function(x)x$beta)
            )

  if(q == 1L){
    out$nsup <- out$nsup[[1]]
    out$lambda <- out$lambda[[1]]
    out$beta <- as.matrix(out$beta[[1]])
  }

  if(flagsave){
    out$file_beta <- gsub("i_[0-9]+.bin$","i_\\*.bin",
                            normalizePath(paste0(file_beta,fileID[1],".bin")))
    out$fileID <- fileID
    out$beta <- NULL
  }

  class(out) <- "LASSO"
  return(out)
}
