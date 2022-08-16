
solveEN <- function(Sigma, Gamma, X = NULL, alpha = 1,
              lambda = NULL, nlambda = 100, common.lambda = TRUE,
              lambda.min = .Machine$double.eps^0.5, dfmax = NULL,
              scale = TRUE, tol = 1E-5, maxiter = 1000, mc.cores = 1L,
              return.beta = TRUE, save.beta = FALSE, verbose = FALSE)
{
    alpha <- as.numeric(alpha)
    scale <- as.logical(scale)
    tol <- as.numeric(tol)
    maxiter <- as.integer(maxiter)

    if(length(dim(Gamma)) != 2L){
      Gamma <- matrix(Gamma, ncol=1L)
    }
    dimnames(Sigma) <- NULL
    p <- nrow(Gamma)
    q <- ncol(Gamma)

    if((sum(dim(Sigma))/2)^2 != p^2){
        stop("Object 'Sigma' must be a p*p squared matrix where p=nrow(Gamma)")
    }

    if(alpha<0 | alpha>1){ stop("Parameter 'alpha' must be a number between 0 and 1")}

    if(!float::storage.mode(Sigma) %in% c("float32","double")) storage.mode(Sigma) <- "double"
    if(!float::storage.mode(Gamma) %in% c("float32","double")) storage.mode(Gamma) <- "double"
    isfloat <- FALSE
    if(float::storage.mode(Sigma)=="float32" | float::storage.mode(Gamma)=="float32"){
          if(float::storage.mode(Sigma)!="float32") Sigma <- float::fl(Sigma)
          if(float::storage.mode(Gamma)!="float32") Gamma <- float::fl(Gamma)
          isfloat <- TRUE
    }

    if(!is.null(X)){
      if(length(dim(X)) != 2L){
        X <- matrix(X, nrow=1L)
      }
      stopifnot(ncol(X) == nrow(Gamma))
    }

    dfmax <- ifelse(is.null(dfmax), p, as.integer(dfmax))

    if(scale){
      sdx <-  sqrt(float::diag(Sigma))
      cov2cor2(Sigma, void=TRUE)     # Equal to Sigma=cov2cor(Sigma) but faster
      Gamma <- float::sweep(Gamma, 1L, sdx, FUN = "/")
    }else{
      sdx <- rep(1, p)
    }
    sdx <-  float::dbl(sdx)

    if(is.null(lambda)){
      if(common.lambda){
        Cmax <- ifelse(alpha > .Machine$double.eps,max(abs(Gamma)/alpha),5)
        lambda <- matrix(exp(seq(log(Cmax),log(lambda.min),length=nlambda)), ncol=1)
      }else{
        lambda <- do.call(cbind,lapply(1:q,function(k){
          Cmax <- ifelse(alpha > .Machine$double.eps,max(abs(Gamma[,k])/alpha),5)
          exp(seq(log(Cmax),log(lambda.min),length=nlambda))
        }))
      }
    }else{
      if(length(dim(lambda))==2L){
        if(ncol(lambda) > 1 & ncol(lambda) < q){
          stop("Number of columns of 'lambda' must be equal to the number of columns in 'Gamma'")
        }
      }else{
        lambda <- matrix(lambda, ncol=1)
      }
      if(any(apply(lambda, 2L, function(x) any(diff(x) > 0)))){
          stop("Object 'lambda' must be a matrix (or vector) of decreasing numbers")
      }
    }
    storage.mode(lambda) <- "double"
    nlambda <- nrow(lambda)
    verbose2 <- (verbose & q==1L)

    if(return.beta & save.beta){
      message("'return.beta' is set to FALSE when 'save.beta=TRUE'")
      return.beta <- FALSE
    }

    compApply <- function(ind)
    {
      if(ncol(lambda)>1L){
        lambda0 <- lambda[,ind]
      }else{
        lambda0 <- lambda[,1]
      }

      #dyn.load("c_utils.so")
      if(isfloat){
        tmp <- .Call('updatebeta', p, Sigma@Data, Gamma@Data[,ind],
                 nlambda, lambda0, alpha, tol, maxiter,
                 dfmax, isfloat, scale, sdx, verbose2)

      }else{
        tmp <- .Call('updatebeta', p, Sigma, Gamma[,ind],
                 nlambda, lambda0, alpha, tol, maxiter,
                 dfmax, isfloat, scale, sdx, verbose2)
      }
      #dyn.unload("c_utils.so")

      df0 <- tmp[[2]]

      if(dfmax < p){
        index <- which(df0 <= dfmax)
        lambda0 <- lambda0[index]
        df0 <- df0[index]
        beta <- Matrix::Matrix(tmp[[1]][, index, drop=FALSE], sparse=TRUE)
      }else{
        beta <- Matrix::Matrix(tmp[[1]], sparse=TRUE)
      }

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

      list(ind=ind, yHat=yHat, beta=beta, lambda=lambda0, df=df0)
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
    if(any(seq(q) != unlist(lapply(out,function(x) x$ind)) )){
        stop("Some sub-processes failed. Something went wrong during the analysis.")
    }

    out <- list(p=p, q=q,
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
