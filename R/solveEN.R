
# P=XtX2; v=Xty; alpha = 1; lambda = NULL; nLambda = 100; scale = TRUE
# tol = 1E-5; maxIter = 1000; verbose = FALSE

solveEN <- function(P, v, alpha = 1, lambda = NULL, nLambda = 100,
    minLambda = .Machine$double.eps^0.5, scale = TRUE,
    tol = 1E-5, maxIter = 1000, verbose = FALSE)
{
    v <- as.vector(v)
    dimnames(P) <- NULL
    p <- length(v)

    if((sum(dim(P))/2)^2 != p^2)
        stop("Object 'P' must be a p*p squared matrix where p=length(v)")

    if(alpha<0 | alpha>1) stop("Parameter 'alpha' must be a number between 0 and 1")

    if(!float::storage.mode(P) %in% c("float32","double")) storage.mode(P) <- "double"
    if(!float::storage.mode(v) %in% c("float32","double")) storage.mode(v) <- "double"
    isFloat <- FALSE
    if(float::storage.mode(P)=="float32" | float::storage.mode(v)=="float32"){
          if(float::storage.mode(P)!="float32") P <- float::fl(P)
          if(float::storage.mode(v)!="float32") v <- float::fl(v)
          isFloat <- TRUE
    }

    if(scale){
      sdx <-  sqrt(float::diag(P))
      P <- cov2cor2(P)  # Equal to P=cov2cor(P) but faster
      v <- v/sdx
    }else{
      sdx <- rep(1,p)
    }

    if(is.null(lambda)){
      Cmax <- ifelse(alpha > .Machine$double.eps, max(abs(v)/alpha), 5)
      lambda <- exp(seq(log(Cmax),log(minLambda),length=nLambda))
    }else{
      if(!is.vector(lambda,mode="numeric") | any(diff(lambda) >0))
        stop("Object 'lambda' must be a vector of decreasing numbers")
    }
    nLambda <- length(lambda)

    #dyn.load("c_utils.so")
    if(isFloat)
    {
      beta <- .Call('updatebeta',as.integer(p),P@Data,v@Data,
               as.integer(nLambda),as.numeric(lambda),
               as.numeric(alpha), as.numeric(tol),
               as.integer(maxIter),verbose,isFloat)[[1]]

    }else{
      beta <- .Call('updatebeta',as.integer(p),P,as.vector(v),
               as.integer(nLambda),as.numeric(lambda),as.numeric(alpha),
               as.numeric(tol),as.integer(maxIter),verbose,isFloat)[[1]]
    }
    #dyn.unload("c_utils.so")

    if(scale) beta <- sweep(beta,2L,as.numeric(sdx),FUN="/")
    df <- do.call(c,lapply(1:nrow(beta),function(i) sum(abs(beta[i,])>0)))

    out <- list(beta=beta,lambda=lambda,df=df)
    class(out) <- "LASSO"
    return(out)
}
