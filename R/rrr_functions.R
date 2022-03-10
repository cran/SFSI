
#====================================================================
# Add the value 'a' to the diagonal of 'V' matix
#====================================================================
add2diag <- function(V, a, void = FALSE)
{
  if((sum(dim(V))/2)^2 != length(V)) stop("Object 'V' must be a squared matrix")
  if(!float::storage.mode(V) %in% c("float32","double")) storage.mode(V) <- "double"

  p <- ncol(V)
  isFloat <- float::storage.mode(V)=="float32"

  #dyn.load("c_utils.so")
  if(void)
  {
    if(isFloat){
      out <- .Call('addvalue2diag',as.integer(p),V@Data,as.numeric(a),isFloat)
    }else{
      out <- .Call('addvalue2diag',as.integer(p),V,as.numeric(a),isFloat)
    }
  }else{
    if(isFloat){
      out <- V@Data[]
    }else out <- V[]

    tmp <- .Call('addvalue2diag',as.integer(p),out,as.numeric(a),isFloat)
    if(isFloat) out <- float::float32(out)
  }
  #dyn.unload("c_utils.so")
  out
}

#====================================================================
# Used by the 'plotPath' function
#====================================================================
getIndexCorrelated <- function(X, maxCor = 0.8)
{
  COV <- stats::cov(X)
  p <- ncol(COV)
  index <- .Call("getCorrelated",as.integer(p),COV,as.numeric(maxCor))
  out <- NULL
  if(index[[2]] > 0L) out <- index[[1]][1:index[[2]]]
  out
}

#====================================================================
# Transpose of the 'backsolve' function
#====================================================================
backsolvet <- function(r, x, k=ncol(r))
{
  float::backsolve(r,x,k,transpose=TRUE)
}

#====================================================================
# Update the lower triangular CHOLESKY decomposition when adding a new column
#====================================================================
# xtx=P[inew,inew]; Xtx=as.vector(P[inew,active])
upDateR <- function(xtx, R = NULL, Xtx, eps = .Machine$double.eps)
{
  norm.xnew <- sqrt(xtx)
  if(is.null(R)) {
    # R <- matrix(norm.xnew, 1, 1)
    R <- float::t(norm.xnew)
    attr(R, "rank") <- 1
    R
  }else{
    r <- backsolvet(R, Xtx)
    rpp <- norm.xnew^2 - sum(r^2)
    rank <- attr(R, "rank")	### check if R is machine singular
    if(rpp <= eps){
      rpp <- eps
    }else{
      rpp <- sqrt(rpp)
      rank <- rank + 1
    }
    if(float::storage.mode(R)=="float32") rpp <- float::fl(rpp)
    R <- cbind(rbind(R, float::fl(0)), c(r, rpp))
    attr(R, "rank") <- rank
    R
  }
}


#====================================================================
# Update the lower triangular CHOLESKY decomposition when deleting one column
#====================================================================
downDateR <- function(R, k = p)
{
	p <- dim(R)[1]
	if(p == 1){
		return(NULL)
  }else{
	   R <- deleteCol(R, rep(1, p), k)[[1]][ - p,  , drop = FALSE]
	   attr(R, "rank") <- p - 1
	   return(R)
  }
}

#====================================================================
# Used by the 'downDateR' function
#====================================================================
deleteCol <- function(R, z, k = p)
{
	p <- dim(R)[1]
	R <- R[, -k, drop = FALSE]
	z <- as.matrix(z)
	pz <- dim(z)[2]
  if(!float::storage.mode(R) %in% c("double","float32")) storage.mode(R) <- "double"
  isFloat <- float::storage.mode(R) == "float32"

  if(isFloat){
    z <- float::fl(z)
  }else storage.mode(z) <- "double"

  #dyn.load("c_utils.so")
	if(isFloat){
  	tmp = .Call("delete_col",R@Data,as.integer(p),as.integer(k),z@Data,as.integer(pz),isFloat)
		return(lapply(tmp,function(x)float::float32(x)))
	}else{
		return(.Call("delete_col",R,as.integer(p),as.integer(k),z,as.integer(pz),isFloat))
	}
  #dyn.unload("c_utils.so")
}

#====================================================================
# Recursive quantities for the GEMMA algorithm
#====================================================================
atPib <- function(i, Uta, Utb, UtX = UtX, dbar = dbar){
  if(i==1) {
    aPw <- sum(Uta*UtX[,i]*dbar)
    bPw <- sum(Utb*UtX[,i]*dbar)
    wPw <- sum(UtX[,i]^2*dbar)
    sum(Uta*Utb*dbar)-aPw*bPw/wPw
  }else{
      atPib(i-1,Uta,Utb,UtX,dbar) - atPib(i-1,Uta,UtX[,i],UtX,dbar)*atPib(i-1,Utb,UtX[,i],UtX,dbar)/atPib(i-1,UtX[,i],UtX[,i],UtX,dbar)
  }
}

atPia <- function(i,Uta, UtX = UtX, dbar = dbar){
  if(i==1) {
    aPw <- sum(Uta*UtX[,i]*dbar)
    wPw <- sum(UtX[,i]^2*dbar)
    sum(Uta^2*dbar)-(aPw^2)/wPw
  }else{
      atPia(i-1,Uta,UtX,dbar) - atPib(i-1,Uta,UtX[,i],UtX,dbar)^2/atPia(i-1,UtX[,i],UtX,dbar)
  }
}

atPiPib <- function(i, Uta, Utb, UtX, dbar){
  if(i==1) {
    aPw <- sum(Uta*UtX[,i]*dbar)
    aPPw <- sum(Uta*UtX[,i]*dbar^2)
    bPw <- sum(Utb*UtX[,i]*dbar)
    bPPw <- sum(Utb*UtX[,i]*dbar^2)
    wPw <- sum(UtX[,i]^2*dbar)
    wPPw <- sum(UtX[,i]^2*dbar^2)
    sum(Uta*Utb*dbar^2)+aPw*bPw*wPPw/(wPw^2) - aPw*bPPw/wPw - bPw*aPPw/wPw
  }else{
      atPiPib(i-1,Uta,Utb,UtX,dbar) +
      atPib(i-1,Uta,UtX[,i],UtX,dbar)*atPib(i-1,Utb,UtX[,i],UtX,dbar)*atPiPib(i-1,UtX[,i],UtX[,i],UtX,dbar)/atPib(i-1,UtX[,i],UtX[,i],UtX,dbar)^2 -
      atPib(i-1,Uta,UtX[,i],UtX,dbar)*atPiPib(i-1,Utb,UtX[,i],UtX,dbar)/atPib(i-1,UtX[,i],UtX[,i],UtX,dbar) -
      atPib(i-1,Utb,UtX[,i],UtX,dbar)*atPiPib(i-1,Uta,UtX[,i],UtX,dbar)/atPib(i-1,UtX[,i],UtX[,i],UtX,dbar)
  }
}

atPiPia <- function(i, Uta, UtX, dbar){
  if(i==1) {
    aPw <- sum(Uta*UtX[,i]*dbar)
    aPPw <- sum(Uta*UtX[,i]*dbar^2)
    wPw <- sum(UtX[,i]^2*dbar)
    wPPw <- sum(UtX[,i]^2*dbar^2)
    sum(Uta^2*dbar^2)+aPw^2*wPPw/(wPw^2) - 2*aPw*aPPw/wPw
  }else{
      atPiPia(i-1,Uta,UtX,dbar) +
      atPib(i-1,Uta,UtX[,i],UtX,dbar)^2*atPiPia(i-1,UtX[,i],UtX,dbar)/atPia(i-1,UtX[,i],UtX,dbar)^2 -
      2*atPib(i-1,Uta,UtX[,i],UtX,dbar)*atPiPib(i-1,Uta,UtX[,i],UtX,dbar)/atPia(i-1,UtX[,i],UtX,dbar)
  }
}

tr_Pi <- function(i, UtX, dbar){
  if(i==1) {
    wPw <- sum(UtX[,i]^2*dbar)
    wPPw <- sum(UtX[,i]^2*dbar^2)
    sum(dbar) - wPPw/wPw
  }else{
    tr_Pi(i-1,UtX,dbar) - atPiPia(i-1,UtX[,i],UtX,dbar)/atPia(i-1,UtX[,i],UtX,dbar)
  }
}

#====================================================================
# Derivative of the Log Likelihood (ML)
#====================================================================
dlogLik <- function(ratio, n, c0, Uty, UtX, d)
{
  dbar <- 1/(ratio*d+1)
  #Tr_Hinv <- sum(dbar)
  Tr_Hinv_G <- (n-sum(dbar))/ratio #(n-Tr_Hinv)/ratio

  ytPy <- atPia(c0+1,Uty,UtX=UtX,dbar=dbar)
  ytPPy <- atPiPia(c0+1,Uty,UtX=UtX,dbar=dbar)

  ytPGPy <- (ytPy-ytPPy)/ratio

  dd <- -0.5*Tr_Hinv_G + 0.5*n * ytPGPy/ytPy

  return(dd)
}

#====================================================================
# Derivative of the Log-restricted Likelihood (REML)
#====================================================================
# tt=dlogResLik(ratio,n=n,c0=c0,Uty=Uty,UtX=UtX,d=d)
dlogResLik <- function(ratio, n, c0, Uty, UtX, d)
{
  dbar <- 1/(ratio*d+1)
  Tr_Px <- tr_Pi(c0+1,UtX=UtX,dbar=dbar)
  Tr_Px_G <- (n-c0-1-Tr_Px)/ratio

  ytPy <-  atPia(c0+1,Uty,UtX=UtX,dbar=dbar)
  ytPPy <- atPiPia(c0+1,Uty,UtX=UtX,dbar=dbar)

  ytPGPy <- (ytPy-ytPPy)/ratio

  dd <- -0.5*Tr_Px_G + 0.5*(n-c0-1)* ytPGPy/ytPy

  return(dd)
}

#====================================================================
# Search for the root ratio=varU/varE in a given interval
#====================================================================
searchInt <- function(method, interval, n, c0, Uty, UtX, d, maxiter, tol, lower, upper, varP)
{
  flag <- TRUE; i <- 1
  convergence <- ratio <- dbar <- varU <- varE <- bHat <- msg <- NA
  while(flag)
  {
    i <- i + 1
    if(method=="REML"){
      tmp <- try(uniroot(f=dlogResLik,interval=c(interval[i-1],interval[i]),n=n,c0=c0,Uty=Uty,
                       UtX=UtX,d=d,tol=tol,maxiter=maxiter,trace=2),
               silent = TRUE)
    }else{
      tmp <- try(uniroot(f=dlogLik,interval=c(interval[i-1],interval[i]),n=n,c0=c0,Uty=Uty,
                        UtX=UtX,d=d,tol=tol,maxiter=maxiter,trace=2),
                silent = TRUE)
    }
    if(inherits(tmp, "list"))
    {
      ratio0 <- tmp$root
      if(ratio0 <= lower){
        ratio0 <- lower
        msg <- paste0("Root varU/varE is the lower bound ",lower)
      }else{
        if(ratio0 >= upper){
          ratio0 <- upper
          msg <- paste0("Root varU/varE is the upper bound ",upper)
        }
      }

      dbar <- 1/(ratio0*d + 1)
      qq1 <- t(Uty*dbar)%*%UtX
      qq2 <- solve(sweep(t(UtX),2L,dbar,FUN="*")%*%UtX)
      ytPy <- drop(sum(dbar*Uty^2)-qq1%*%qq2%*%t(qq1))
      bHat <- drop(qq2%*%t(qq1))

      varE <- ifelse(method=="REML",ytPy/(n-c0-1),ytPy/n)
      varU <- ratio0*varE

      if(varU <= (2)*varP){  # A quality control-like
        convergence <- tmp$iter <= maxiter
        ratio <-  ratio0
      }
    }
    #aa <- rep(NA,3)
    #if(class(tmp) == "list") aa=c(tmp$root,tmp$f.root,tmp$estim.prec)
    #cat("Interval ",i-1,"[",interval[i-1],",",interval[i],"]: root=",aa[1]," f.root=",aa[2]," prec=",aa[3],"\n")
    if(i == length(interval) | !is.na(convergence)) flag <- FALSE
  }
  list(ratio=ratio,varU=varU,varE=varE,convergence=convergence,
       dbar=dbar,bHat=bHat,msg=msg)
}

#====================================================================
# Labels and breaks for the DF axis
#====================================================================
# x = dat$lambda; y = dat$df
get_breaks <- function(x, y, nbreaks = 6, ymin = 1)
{
  p <- max(round(y))
  #yy <- c(ymin,p-1,p)
  yy <- c(ymin,p)
  neglogx <- -log(x)

  fm <- stats::smooth.spline(neglogx,y)

  tmp <- cbind(neglogx,stats::fitted(fm))
  tmp <- tmp[order(tmp[,1]),]
  xxmin <- tmp[min(which(tmp[,2] >= ymin)),1]

  tmp <- seq(min(neglogx), ifelse(xxmin<=0,0,2*xxmin), length=1000)
  tt <- stats::predict(fm,tmp)
  xxmin <- tt$x[min(which(tt$y >= ymin))]

  #xx <- stats::predict(stats::smooth.spline(y, neglogx),yy)$y
  #tmp <- data.frame(x=exp(-xx),neglogx=xx,y=yy)
  #breaks.x <- seq(tmp[1,2], -log(min(x)), length=nbreaks)
  breaks.x <- seq(xxmin, -log(min(x)), length=nbreaks)
  breaks.y <- stats::predict(fm, breaks.x)$y

  return(list(breaks.x=breaks.x,breaks.y=breaks.y))
}


has_names <- function(X){
  if(length(dim(X)) == 2L){
    out <- length(unlist(dimnames(X))) == sum(dim(X))
  }else out <- FALSE

  out
}

#====================================================================
# Obtain layout to for the net.plot function
#====================================================================
get_net <- function(X, K = NULL, xxx = NULL, yyy = NULL, eps = .Machine$double.eps)
{
  labelsAxis <- labels0 <- NULL
  isSymmetric <- isSymmetric(X, tol=1E-6)
  uniqueNames <- unique(c(rownames(X), colnames(X)))

  if(is.null(xxx) | is.null(yyy))
  {
    if(isSymmetric){
      xxx <- yyy <- 1:nrow(X)
    }else{
      if(has_names(X)){
        if(!is.null(K) & has_names(K)){
          if(all(uniqueNames %in% colnames(K))){
            yyy <- match(rownames(X),rownames(K))
            xxx <- match(colnames(X),rownames(K))
          }else{
            yyy <- match(rownames(X),uniqueNames)
            xxx <- match(colnames(X),uniqueNames)
            cat("Some row/column names of 'object' were not found in row names of 'K'.",
                "\nInput 'K' will be ignored\n")
            K <- NULL
          }
        }else{
          yyy <- match(rownames(X),uniqueNames)
          xxx <- match(colnames(X),uniqueNames)
        }
      }else{
       yyy <- 1:nrow(X)
       xxx <- nrow(X)+(1:ncol(X))
      }
    }
  }

  if(!is.null(K)){
    if(all(dim(X) == dim(K))){
      if(has_names(X) & !has_names(K)){
          dimnames(K) <- dimnames(X)
      }
      if(!has_names(X) & has_names(K)){
          dimnames(X) <- dimnames(K)
      }

    }else{
      if((has_names(X) + has_names(K)) <= 1){
        cat("Input 'object' couldn't be linked to 'K' through row/column names.",
            "\nInput 'K' will be ignored\n")
        K <- NULL
      }
    }
    if(has_names(X) & has_names(K)){
      if(any(!uniqueNames %in% rownames(K))){
        cat("Some row/column names of 'object' were not found in row names of 'K'.",
          "\nInput 'K' will be ignored\n")
          K <- NULL
      }
    }
  }

  # Labels
  if(isSymmetric){
    if(has_names(X)){
      if(is.null(K)){
        labels0 <- rownames(X)
      }else labels0 <- rownames(K)
    }else labels0 <- 1:ncol(X)
  }else{
     if(has_names(X)){
       if(is.null(K)){
         labels0 <- rep("",length(uniqueNames))
         labels0[yyy] <- rownames(X)
         labels0[xxx] <- colnames(X)
      }else labels0 <- rownames(K)
     }else{
       labels0 <- rep("",nrow(X)+ncol(X))
       labels0[yyy] <- paste0("R",seq_along(yyy))
       labels0[xxx] <- paste0("C",seq_along(xxx))
     }
  }

  # Get edges
  edges <- vector("list",length(yyy))
  names(edges) <- yyy
  for(i in 1:nrow(X))
  {
    edges[[i]] <- NA
    if(isSymmetric){
      edges[[i]] <- which(abs(X[i,i:ncol(X)]) > eps) + i -1
    }else{
      edges[[i]]  <- xxx[which(abs(X[i, ]) > eps)]
    }
  }

  # Sets: 1=Row, 2=Connected column, 3=Non-connected column, 4=1&2
  set <- rep(3,length(labels0))

  tmp <- unique(unlist(edges))
  set[tmp[!tmp%in%yyy]] <- 2
  set[yyy] <- 1
  tmp <- intersect(tmp,yyy)
  if(!isSymmetric & length(tmp) > 0) set[tmp] <- 4

  if(is.null(K))
  {
    eee <- c()
    if(isSymmetric){
      gr <- igraph::make_empty_graph(n = ncol(X))
      for(i in 1:nrow(X)){
        index <- edges[[i]][-1]
        if(length(index[!is.na(index)]) > 0) eee <- c(eee, as.vector(rbind(i,index)))
      }
    }else{
      if(has_names(X)){
          n0 <- length(uniqueNames)
      }else n0 <- ncol(X)+nrow(X)
      gr <- igraph::make_empty_graph(n = n0)
      for(i in 1:nrow(X)){
        index <- edges[[i]]
        if(length(index[!is.na(index)]) > 0)   eee <- c(eee, as.vector(rbind(yyy[i],index)))
      }
    }

    gr <- igraph::add_edges(gr, eee)
    xy <- igraph::layout_with_fr(gr, dim=2, niter=2000)

  }else{
    tmp <- float::svd(K,nu=2,nv=0)
    d <- tmp$d
    xy <- tmp$u[,1:2]

    expvarPC <- paste0(" (",sprintf('%.1f',100*d/sum(d)),"%)")
    if(length(expvarPC) < 2) expvarPC <- NULL
    labelsAxis <- paste0("PC ",1:2,expvarPC[1:2])
  }

  isEigen <- !is.null(K)
  colnames(xy) <- c("x","y")

  return(list(xy=xy, labels=as.character(labels0), labelsAxis=labelsAxis,
              isSymmetric=isSymmetric, xxx=xxx, yyy=yyy, set=set,
              isEigen=isEigen, edges=edges))
}

#====================================================================
#====================================================================
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("
  |=======================================================================|
  |    ._______. ._______. ._______. ._______.                            |
  |    | ._____| | ._____| | ._____| |__. .__|                            |
  |    | |_____. | |___.   | |_____.    | |                               |
  |    |_____. | | .___|   |_____. |    | |      Authors:                 |
  |    ._____| | | |       ._____| | .__| |__.     Marco Lopez-Cruz       |
  |    |_______| |_|       |_______| |_______|     Gustavo de los Campos  |
  |                                                                       |
  |    Sparse Family and Selection Index. Version 1.1.0 (Mar 10, 2022)    |
  |    Type 'citation('SFSI')' to know how to cite SFSI                   |
  |    Type 'help(package='SFSI',help_type='html')' to see help           |
  |    Type 'browseVignettes('SFSI')' to see documentation                |
  |    Type 'demo(package='SFSI')' to see demos                           |
  |                                                                       |
  |=======================================================================|
  ")
}
