
solveEN <- function(Sigma, Gamma, alpha = 1, lambda = NULL,
                    nlambda = 100, lambda.min = .Machine$double.eps^0.5,
                    lambda.max = NULL, common.lambda = TRUE, beta0 = NULL,
                    nsup.max = NULL, scale = TRUE, sdx = NULL, tol = 1E-5,
                    maxiter = 1000, mc.cores = 1L, save.at = NULL,
                    precision.format = c("double","single"),
                    fileID = NULL, verbose = FALSE)
{
    precision.format <- match.arg(precision.format)
    alpha <- as.numeric(alpha)
    scale <- as.logical(scale)
    tol <- as.numeric(tol)
    maxiter <- as.integer(maxiter)

    if(length(dim(Gamma)) == 2L){
      Gamma <- as.matrix(Gamma)
    }else{
      Gamma <- matrix(Gamma, ncol=1L)
    }
    p <- nrow(Gamma)
    q <- ncol(Gamma)

    if((sum(dim(Sigma))/2)^2 != p^2){
      stop("Input 'Sigma' must be a squared matrix of dimension equal to nrow(Gamma)")
    }

    if(alpha<0 | alpha>1){ stop("Parameter 'alpha' must be a number between 0 and 1")}
    stopifnot(maxiter>0)
    if(tol < .Machine$double.eps){
      stop("Input 'tol' must be > 0")
    }
    nsup.max <- ifelse(is.null(nsup.max), p, as.integer(nsup.max))

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

    # Get lambda grid. Diagonal values in Sigma are assumed to be zero
    lambda <- setLambda(Gamma, alpha=alpha, lambda=lambda, nlambda=nlambda,
                        lambda.min=lambda.min, lambda.max=lambda.max,
                        common.lambda=common.lambda,verbose=FALSE)
    nlambda <- nrow(lambda)

    if(ifelse(is.null(beta0),FALSE,length(beta0)!=p)){
      stop("Input 'beta0' must be a numeric vector of length = ",p)
    }

    flagsave <- as.logical(!is.null(save.at))
    verbose2 <- as.logical(verbose & q==1L)
    mc.cores <- ifelse(q==1L & mc.cores>1L, 1L, mc.cores)
    doubleprecision <- as.logical(precision.format=="double")

    compApply <- function(ind)
    {
      if(ncol(lambda)>1L){
        lambda0 <- lambda[,ind]
      }else{
        lambda0 <- lambda[,1]
      }

      if(flagsave){
        filename <- paste0(file_beta,fileID[ind],".bin")
      }else{
        filename <- NULL
      }

      #dyn.load("c_solveEN.so")
      res <- .Call('R_updatebeta', Sigma, Gamma[,ind],
                  lambda0, alpha, beta0, tol, maxiter, nsup.max,
                  scaleb, sdx, filename,
                  doubleprecision, verbose2)
      #dyn.unload("c_solveEN.so")

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
        dir.create(dirname(file_beta), recursive=TRUE)
      }

      if(is.null(fileID)){
        fileID <- seq(q)
      }else{
        stopifnot(length(fileID) == q)
      }
    }

    # Run the analysis for 1:ncol(Gamma)
    if(verbose & q>1L){
      pb <- utils::txtProgressBar(style=3)
      con <- tempfile(tmpdir=tmpdir0)
    }
    if(mc.cores == 1L){
      out <- lapply(X=seq(q), FUN=compApply)
    }else{
      out <- parallel::mclapply(X=seq(q), FUN=compApply, mc.cores=mc.cores)
    }
    if(verbose & q>1L) {
      close(pb); unlink(con)
    }

    # Checkpoint
    if(any(seq(q) != unlist(lapply(out,function(x) x$ind)) )){
      stop("Some sub-processes failed. Something went wrong during the analysis.")
    }

    out <- list(p=p, q=q, nlambda=nlambda,
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
