
# Sigma=XtX; Gamma=Xty; alpha = 1; lambda = NULL; nlambda = 100; scale = TRUE
# tol = 1E-5; maxiter = 1000; verbose = FALSE

solveEN <- function(Sigma, Gamma, alpha = 1, lambda = NULL, nlambda = 100,
    lambda.min = .Machine$double.eps^0.5, dfmax = NULL, scale = TRUE,
    tol = 1E-5, maxiter = 1000, verbose = FALSE)
{
    Gamma <- as.vector(Gamma)
    dimnames(Sigma) <- NULL
    p <- length(Gamma)

    if((sum(dim(Sigma))/2)^2 != p^2)
        stop("Object 'Sigma' must be a p*p squared matrix where p=length(Gamma)")

    if(alpha<0 | alpha>1) stop("Parameter 'alpha' must be a number between 0 and 1")

    if(!float::storage.mode(Sigma) %in% c("float32","double")) storage.mode(Sigma) <- "double"
    if(!float::storage.mode(Gamma) %in% c("float32","double")) storage.mode(Gamma) <- "double"
    isFloat <- FALSE
    if(float::storage.mode(Sigma)=="float32" | float::storage.mode(Gamma)=="float32"){
          if(float::storage.mode(Sigma)!="float32") Sigma <- float::fl(Sigma)
          if(float::storage.mode(Gamma)!="float32") Gamma <- float::fl(Gamma)
          isFloat <- TRUE
    }

    if(scale){
      sdx <-  sqrt(float::diag(Sigma))
      Sigma <- cov2cor2(Sigma)  # Equal to Sigma=cov2cor(Sigma) but faster
      Gamma <- Gamma/sdx
    }else{
      sdx <- rep(1,p)
    }

    if(is.null(dfmax)) dfmax <- p

    if(is.null(lambda)){
      Cmax <- ifelse(alpha > .Machine$double.eps, max(abs(Gamma)/alpha), 5)
      lambda <- exp(seq(log(Cmax),log(lambda.min),length=nlambda))
    }else{
      if(length(dim(lambda))==2 | mode(lambda)!="numeric" | any(diff(lambda)>.Machine$double.eps^0.5))
        stop("Object 'lambda' must be a vector of decreasing numbers")
      dfmax <- p
    }
    nlambda <- length(lambda)

    #dyn.load("c_utils.so")
    if(isFloat)
    {
      beta <- .Call('updatebeta',as.integer(p),Sigma@Data,Gamma@Data,
               as.integer(nlambda),as.numeric(lambda),
               as.numeric(alpha), as.numeric(tol),as.integer(maxiter),
               as.integer(dfmax),verbose,isFloat)

    }else{
      beta <- .Call('updatebeta',as.integer(p),Sigma,as.vector(Gamma),
               as.integer(nlambda),as.numeric(lambda),
               as.numeric(alpha),as.numeric(tol),as.integer(maxiter),
               as.integer(dfmax),verbose,isFloat)
    }
    #dyn.unload("c_utils.so")

    df <- beta[[2]]
    beta <- beta[[1]]

    index <- which(df <= dfmax)
    lambda <- lambda[index]
    df <- df[index]
    beta <- beta[, index, drop=FALSE]

    if(scale) beta <- sweep(beta,1L,as.numeric(sdx),FUN="/")

    out <- list(beta=beta,lambda=lambda,df=df)
    class(out) <- "LASSO"
    return(out)
}
