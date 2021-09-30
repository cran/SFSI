#====================================================================
# Covariance matrix to distance matrix (user-level)
#====================================================================
cov2dist <- function(V,void=FALSE)
{
    if((sum(dim(V))/2)^2 != length(V)) stop("Object 'V' must be a squared matrix")
    if(!float::storage.mode(V) %in% c("double","float32")) storage.mode(V) <- "double"

    p <- ncol(V)
    isFloat <- float::storage.mode(V)=="float32"

    #dyn.load("c_utils.so")
    if(void)
    {
     if(isFloat){
      out <- .Call('cov2distance',as.integer(p),V@Data,isFloat)
     }else{
      out <- .Call('cov2distance',as.integer(p),V,isFloat)
     }
    }else{
      if(isFloat){
       out <- V@Data[]
     }else out <- V[]

     tmp <- .Call('cov2distance',as.integer(p),out,isFloat)
     if(isFloat) out <- float::float32(out)
   }
   #dyn.unload("c_utils.so")
   out
}

#====================================================================
# Covariance matrix to correlation matrix (user-level)
#====================================================================
cov2cor2 <- function(V,a=1,void=FALSE)
{
    if((sum(dim(V))/2)^2 != length(V)) stop("Object 'V' must be a squared matrix")
    if(!float::storage.mode(V) %in% c("double","float32")) storage.mode(V) <- "double"

    p <- ncol(V)
    isFloat <- float::storage.mode(V)=="float32"

    #dyn.load("c_utils.so")
    if(void)
    {
      if(isFloat){
       nOK <- .Call('cov2correlation',as.integer(p),V@Data,isFloat,as.numeric(a))[[1]]
      }else{
       nOK <- .Call('cov2correlation',as.integer(p),V,isFloat,as.numeric(a))[[1]]
      }
      out <- NULL
    }else{
      if(isFloat){
       out <- V@Data[]
      }else out <- V[]

     nOK <- .Call('cov2correlation',as.integer(p),out,isFloat,as.numeric(a))[[1]]
     if(isFloat) out <- float::float32(out)
   }
   #dyn.unload("c_utils.so")
   if(nOK != p) warning("Some diagonal values of 'V' are 0 or NA. Results are dobubtful",immediate.=TRUE)
   out
}

#====================================================================
# Add the value 'a' to the diagonal of 'V' matix
#====================================================================
add2diag <- function(V,a,void=FALSE)
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
getIndexCorrelated <- function(X,maxCor=0.8)
{
  COV <- stats::cov(X)
  p <- ncol(COV)
  index <- .Call("getCorrelated",as.integer(p),COV,as.numeric(maxCor))
  out <- NULL
  if(index[[2]]>0) out <- index[[1]][1:index[[2]]]
  out
}

#====================================================================
# Collect all outputs when divided acording to 'subset' parameter (user-level)
#====================================================================
collect <- function(prefix="")
{
  filenames <- Sys.glob(paste0(prefix,"_*_of_*.RData"))
  out <- NULL
  if(length(filenames)>0){
      nFiles <- as.numeric(unlist(lapply(strsplit(filenames,"_"),function(x) gsub(".RData","",x[length(x)]))))
      if(length(unique(nFiles))>1)
        stop(" Different subset output files were found for the given prefix='",prefix,
            "'. Remove old files. No output was collected")

      filenames <- paste0(prefix,"_",1:nFiles[1],"_of_",nFiles[1],".RData")
      if(!all(file.exists(filenames))) stop("Some files are missing for the given prefix='",prefix,"'\n")

      for(i in seq_along(filenames))
      {
        load(filenames[i])
        if(i==1){
          fm <- out
        }else{
          fm$file_beta <- c(fm$file_beta,out$file_beta)
          fm$tst <- c(fm$tst,out$tst)
          fm$df <- rbind(fm$df,out$df)
          fm$lambda <- rbind(fm$lambda,out$lambda)
        }
        cat(" Loaded file: '",filenames[i],"'\n",sep="")
      }
      fm$subset[1] <- NA

  }else stop(" No output files were found for the given prefix='",prefix,"'")

  fm
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
# Save a file as binary (user-level)
#====================================================================
saveBinary <- function(X,file = paste0(tempdir(),"/file.bin"),
              type = c("float","double"), verbose = TRUE)
{
  type <- match.arg(type)

  if(length(dim(X)) != 2L) stop("Object 'X' must be a matrix")
  if(!float::storage.mode(X) %in% c("double","float32")) storage.mode(X) <- "double"

  unlink(file)

  if(float::storage.mode(X)=="float32" & type!='float'){
    type <- 'float'
    warning("Object can be only saved as type='float' when class(X)='float'\n",
            "  Variable type was changed to type='float'",immediate.=TRUE)
  }
  isFloat <- float::storage.mode(X)=="float32"
  size <- ifelse(type=="float",4,8)

  if(isFloat){
    out <- .Call('writeBinFileFloat',file,nrow(X),ncol(X),
           as.integer(size),X@Data,isFloat)
   }else{
    out <- .Call('writeBinFileFloat',file,nrow(X),ncol(X),
           as.integer(size),X,isFloat)
  }

  if(verbose){
    tmp <- c(Gb=1E9,Mb=1E6,Kb=1E3,b=1E0)
    sz <- file.size(file)/tmp[min(which(file.size(file)/tmp>1))]
    cat("Saved file '",file,"'\n")
    cat("     nrow=",nrow(X),", ncol=",ncol(X),", type=",type,", size=",size,"bytes, file.size=",round(sz,2),names(sz),"\n")
  }
}

#====================================================================
# Read a file as binary (user-level)
#====================================================================
readBinary <- function(file = paste0(tempdir(),"/file.bin"),
                  indexRow = NULL, indexCol = NULL, verbose = TRUE)
{
  if(!file.exists(file)){
    stop("File '",file,"' does not exist")
  }

  nsetRow <- as.integer(length(indexRow))
  nsetCol <- as.integer(length(indexCol))

  # Read lines
  X <- .Call("readBinFileFloat",file,nsetRow,nsetCol,
             as.integer(indexRow),as.integer(indexCol))

  n <- X[[1]]; p <- X[[2]]; size <- X[[3]]
  isFloat <- X[[4]]
  nError <- X[[5]]

  if(nError==0){
    if(isFloat | size==4){
      X <- float::float32(X[[6]])
      type <- "float"
    }else{
      X <- X[[6]]
      type <- "double"
    }
    if(verbose){
      tmp <- c(Gb=1E9,Mb=1E6,Kb=1E3,b=1E0)
      sz <- object.size(X)/tmp[min(which(object.size(X)/tmp>1))]
      cat("Loaded file '",file,"'\n")
      cat("     nrow=",n,", ncol=",p,", type=",type,", size=",size,"bytes, object.size=",round(sz,2),names(sz),"\n")
    }

  }else{
    X <- NULL
  }
  return(X)
}

#====================================================================
# Recursive quantities for the GEMMA algorithm
#====================================================================
atPib <- function(i,Uta,Utb,UtX=UtX,dbar=dbar){
  if(i==1) {
    aPw <- sum(Uta*UtX[,i]*dbar)
    bPw <- sum(Utb*UtX[,i]*dbar)
    wPw <- sum(UtX[,i]^2*dbar)
    sum(Uta*Utb*dbar)-aPw*bPw/wPw
  }else{
      atPib(i-1,Uta,Utb,UtX,dbar) - atPib(i-1,Uta,UtX[,i],UtX,dbar)*atPib(i-1,Utb,UtX[,i],UtX,dbar)/atPib(i-1,UtX[,i],UtX[,i],UtX,dbar)
  }
}

atPia <- function(i,Uta,UtX=UtX,dbar=dbar){
  if(i==1) {
    aPw <- sum(Uta*UtX[,i]*dbar)
    wPw <- sum(UtX[,i]^2*dbar)
    sum(Uta^2*dbar)-(aPw^2)/wPw
  }else{
      atPia(i-1,Uta,UtX,dbar) - atPib(i-1,Uta,UtX[,i],UtX,dbar)^2/atPia(i-1,UtX[,i],UtX,dbar)
  }
}

atPiPib <- function(i,Uta,Utb,UtX,dbar){
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

atPiPia <- function(i,Uta,UtX,dbar){
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

tr_Pi <- function(i,UtX,dbar){
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
dlogLik <- function(lambda,n,c0,Uty,UtX,d)
{
  dbar <- 1/(lambda*d+1)
  #Tr_Hinv <- sum(dbar)
  Tr_Hinv_G <- (n-sum(dbar))/lambda #(n-Tr_Hinv)/lambda

  ytPy <- atPia(c0+1,Uty,UtX=UtX,dbar=dbar)
  ytPPy <- atPiPia(c0+1,Uty,UtX=UtX,dbar=dbar)

  ytPGPy <- (ytPy-ytPPy)/lambda

  dd <- -0.5*Tr_Hinv_G + 0.5*n * ytPGPy/ytPy

  return(dd)
}

#====================================================================
# Derivative of the Log-restricted Likelihood (REML)
#====================================================================
# tt=dlogResLik(lambda,n=n,c0=c0,Uty=Uty,UtX=UtX,d=d)
dlogResLik <- function(lambda,n,c0,Uty,UtX,d)
{
  dbar <- 1/(lambda*d+1)
  Tr_Px <- tr_Pi(c0+1,UtX=UtX,dbar=dbar)
  Tr_Px_G <- (n-c0-1-Tr_Px)/lambda

  ytPy <-  atPia(c0+1,Uty,UtX=UtX,dbar=dbar)
  ytPPy <- atPiPia(c0+1,Uty,UtX=UtX,dbar=dbar)

  ytPGPy <- (ytPy-ytPPy)/lambda

  dd <- -0.5*Tr_Px_G + 0.5*(n-c0-1)* ytPGPy/ytPy

  return(dd)
}

#====================================================================
# Search for the root lambda0=varU/varE in a given interval
#====================================================================
searchInt <- function(method,interval,n,c0,Uty,UtX,d,maxIter,tol,lower,upper,varP)
{
  flag <- TRUE; i <- 1
  convergence <- lambda0 <- dbar <- varU <- varE <- bHat <- msg <- NA
  while(flag)
  {
    i <- i + 1
    if(method=="REML"){
      tmp <- try(uniroot(f=dlogResLik,interval=c(interval[i-1],interval[i]),n=n,c0=c0,Uty=Uty,
                       UtX=UtX,d=d,tol=tol,maxiter=maxIter,trace=2),
               silent = TRUE)
    }else{
      tmp <- try(uniroot(f=dlogLik,interval=c(interval[i-1],interval[i]),n=n,c0=c0,Uty=Uty,
                        UtX=UtX,d=d,tol=tol,maxiter=maxIter,trace=2),
                silent = TRUE)
    }
    if(class(tmp) == "list")
    {
      lambda00 <- tmp$root
      if(lambda00 <= lower){
        lambda00 <- lower
        msg <- paste0("Root varU/varE is the lower bound ",lower)
      }else{
        if(lambda00 >= upper){
          lambda00 <- upper
          msg <- paste0("Root varU/varE is the upper bound ",upper)
        }
      }

      dbar <- 1/(lambda00*d + 1)
      qq1 <- t(Uty*dbar)%*%UtX
      qq2 <- solve(sweep(t(UtX),2L,dbar,FUN="*")%*%UtX)
      ytPy <- drop(sum(dbar*Uty^2)-qq1%*%qq2%*%t(qq1))
      bHat <- drop(qq2%*%t(qq1))

      varE <- ifelse(method=="REML",ytPy/(n-c0-1),ytPy/n)
      varU <- lambda00*varE

      if(varU <= (2)*varP){  # A quality control-like
        convergence <- tmp$iter <= maxIter
        lambda0 <-  lambda00
      }
    }
    #aa <- rep(NA,3)
    #if(class(tmp) == "list") aa=c(tmp$root,tmp$f.root,tmp$estim.prec)
    #cat("Interval ",i-1,"[",interval[i-1],",",interval[i],"]: root=",aa[1]," f.root=",aa[2]," prec=",aa[3],"\n")
    if(i == length(interval) | !is.na(convergence)) flag <- FALSE
  }
  list(lambda0=lambda0,varU=varU,varE=varE,convergence=convergence,
       dbar=dbar,bHat=bHat,msg=msg)
}

#====================================================================
# Labels and breaks for the DF axis
#====================================================================
getSecondAxis <- function(lambda,df,maxLength=6)
{
  loglambda <- -log(lambda)
  labels0 <- sort(unique(round(df)))
  if(min(labels0)<1) labels0[which.min(labels0)] <- 1

  if(stats::IQR(df)>0)
  {
    breaks0 <- stats::predict(stats::smooth.spline(df, loglambda),labels0)$y
  }else breaks0 <- NULL

  index <- 1
  while(any((breaks0-breaks0[max(index)])>1)){
    dd <- breaks0-breaks0[max(index)]
    index <- c(index,which(dd > 1)[1])
  }
  breaks0 <- breaks0[index]
  labels0 <- labels0[index]

  if(length(breaks0)>maxLength){
    index <- unique(round(seq(1,length(breaks0),length=maxLength)))
    breaks0 <- breaks0[index]
    labels0 <- labels0[index]
  }

  return(list(breaks=breaks0,labels=labels0))
}

#====================================================================
# Plot the top 2 PCs of the K matrix showing tst and trn points (user-level)
#====================================================================
# Z = NULL; subsetG = tst = U = d = group = group.shape = set.color = set.size = df = NULL
# axis.labels = TRUE; curve = FALSE; bg.color = "gray20"; unified = TRUE; ntst = 36;
# line.color = "gray90"; line.tick = 0.3; legend.pos="right";
# point.color = "gray20"; sets = c("Testing","Supporting","Non-active")

plotNet <- function(fm, B, Z = NULL, K, subsetG = NULL, tst = NULL,
           U = NULL, d = NULL, group = NULL, group.shape = NULL,
           set.color = NULL, set.size = NULL, df = NULL, title, axis.labels = TRUE,
           curve = FALSE, bg.color = "gray20", unified = TRUE, ntst = 36,
           line.color = "gray90", line.tick = 0.3, legend.pos="right",
           point.color = "gray20", sets = c("Testing","Supporting","Non-active"))
{
  set <- PC1 <- PC2 <- PC1_TRN <- PC1_TST <- PC2_TRN <- PC2_TST <- NULL
  legend.pos <- match.arg(legend.pos,
    choices=c("right","bottomright","bottomleft","topleft","topright","none"))

  if(!inherits(fm, "SSI")) stop("Object 'fm' is not of the class 'SSI'")

  if(is.null(U) & is.null(d))
  {
    if(is.character(K)){
      K <- readBinary(K)
    }
    if(is.null(K))
      stop("Matrix 'K' must be a positive semi definite matrix\n")
    if(!is.null(Z)) {
      if(length(dim(Z))!=2) stop("Object 'Z' must be a matrix with ncol(Z)=nrow(K)\n")
      K <- float::tcrossprod(Z,float::tcrossprod(Z,K))   # Z%*%K%*%t(Z)
    }
    tmp <- float::svd(K,nu=2,nv=0)
    d <- tmp$d
    U <- tmp$u
    expvarPC <- 100*d/sum(d)  # 100*d/sum(float::diag(K))
  }else{
    if(is.null(U)){
      stop("You are providing the eigevalues, but not the eigenvectors")
    }else{
      if(is.null(d)){
        message("You are providing the eigenvectors, but not the eigenvalues\n",
                "No variance explained can be calculated")
        expvarPC <- NULL
      }else{
        if(nrow(U) == length(d)){
          expvarPC <- 100*d/sum(d)
        }else expvarPC <- NULL
      }
    }
  }
  tmp <- paste0(" (",sprintf('%.1f',expvarPC),"%)")
  if(length(tmp)<2) tmp <- NULL
  labelsPC <- paste0("PC ",1:2,tmp[1:2])

  if(!is.null(tst)){
      if(any(!tst %in% fm$tst))
          stop("Some elements in 'tst' vector are not contained in set 'fm$tst'")
  }else tst <- fm$tst
  if(!unified & length(tst) >= ntst){
   cat("Large number of testing individuals. Only the first",ntst,"are shown\n")
   tst <- tst[1:ntst]
  }

  justx <- ifelse(length(grep("left",legend.pos))>0,0,1)
  justy <- ifelse(length(grep("bottom",legend.pos))>0,0,1)
  if(!legend.pos %in% c("none","right")) legend.pos <- c(abs(justx-0.01),abs(justy-0.01))

  theme0 <- ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    legend.box.spacing = ggplot2::unit(0.4, "lines"),
    legend.background = ggplot2::element_rect(fill = "gray95"),
    panel.background = ggplot2::element_rect(fill = bg.color),
    legend.justification = c(justx,justy),
    legend.position=legend.pos,
    legend.key.height = ggplot2::unit(0.9,"line"),
    legend.key.width = ggplot2::unit(0.9, "lines"),
    legend.title = ggplot2::element_blank(),
    legend.margin = ggplot2::margin(t=0,b=0.25,l=0.25,r=0.25,unit='line'),
    strip.text = ggplot2::element_blank(), panel.spacing = ggplot2::unit(0.1,"lines")
  )

  if(missing(B)){
    if(is.null(df)) df <- summary.SSI(fm)$optCOR$df
    if(0 > df | df > range(fm$df)[2])
      stop("Parameter 'df' must be greater than zero and no greater than nTRN")
    B <- as.matrix(coef.SSI(fm,df=df))
  }else{
    stopifnot(length(dim(B))==2L)
    df <- mean(do.call(c, lapply(1:nrow(B), function(i) sum(abs(B[i,]) > 0))))
  }

  flagGp <- !is.null(group)
  if(is.null(group)) group <- data.frame(group=rep(1,nrow(U)))
  gpName <- colnames(group)
  if(is.null(subsetG)) subsetG <- 1:nrow(U)

  if(!(class(sets) == "character" & length(sets) == 3))
   stop("Parameter 'sets' must be a triplet of 'character' type")

  dat <- data.frame(id=1:nrow(U),set=sets[3],group=group,float::dbl(U[,1:2]))
  dat$set <- as.character(dat$set)

  # Testing and training (active) set
  dat$set[subsetG[tst]] <- sets[1]
  index <- do.call(c, lapply(1:ncol(B), function(j) any(abs(B[fm$tst %in% tst,,drop=FALSE][,j]) > 0)))
  dat$set[subsetG[fm$trn[index]]] <- sets[2]
  dat$set[subsetG[fm$trn[!index]]] <- sets[3]

  colnames(dat) <- c("id","set","group","PC1","PC2")

  dat$group <- factor(as.character(dat$group))
  dat$set <- factor(dat$set,levels=c(sets))

  # Shape and color for the levels of group
  if(!flagGp) dat$group <- dat$set
  levelsGp <- levels(dat$group)
  if(length(levelsGp) > 5)
   stop("Number of levels of 'group' must be at most 5")

  if(is.null(group.shape)){
    if(flagGp){
      group.shape <- c(21,22,23,24,25)
    }else group.shape <- c(21,21,21)
  }
  group.shape <- group.shape[1:length(levelsGp)]

  if(is.null(set.color)){
    set.color <- c("#E69F00","#56B4E9","#999999")
  }
  set.color <- set.color[1:length(sets)]

  if(is.null(set.size)){
    set.size <- c(2.5,1.5,1)
  }
  set.size <- set.size[1:length(sets)]

  if(any(is.na(group.shape)))
    stop("The number of elements in 'group.shape' must be of length ",length(levelsGp))

  if(any(is.na(set.size)) | any(is.na(set.color)))
    stop("The number of elements in 'set.size' and 'set.color' must be of length ",length(sets))

  if(missing(title)){
     title0 <- bquote(.(fm$name)*". Support set size="*.(round(df)))
     theme0 <- theme0 + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }else{
    title0 <- title
    if(is.null(title)){
      theme0 <- theme0 + ggplot2::theme(plot.title = ggplot2::element_blank())
    }else{
      theme0 <- theme0 + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
  }

  if(!axis.labels){
    theme0 <- theme0 + ggplot2::theme(axis.text=ggplot2::element_blank(),
                        axis.ticks=ggplot2::element_blank())
  }

  names(group.shape) <- levelsGp
  names(set.color) <- names(set.size) <- sets

  # If unified plot
  if(unified)
  {
    pt <- ggplot2::ggplot(dat,ggplot2::aes(x=PC1,y=PC2)) +
           ggplot2::geom_point(data=dat[dat$set==sets[3],],ggplot2::aes(shape=group,fill=set),
                    color=point.color,size=set.size[3])

    for(i in 1:length(tst))
    {
      indexTRN <- which(abs(B[which(fm$tst == tst[i]), ]) > 0)
      if(length(indexTRN)>0)
      {
        dat1 <- dat[subsetG[fm$trn],c("PC1","PC2")][indexTRN,]
        dat2 <- dat[subsetG[tst],c("PC1","PC2")][i,]
        colnames(dat1) <- paste0(colnames(dat1),"_TRN")
        colnames(dat2) <- paste0(colnames(dat2),"_TST")
        dat1 <- data.frame(dat2[rep(1,nrow(dat1)),],dat1)
        if(curve){
          pt <- pt + ggplot2::geom_curve(ggplot2::aes(x=PC1_TST,y=PC2_TST,xend=PC1_TRN,yend=PC2_TRN),
                        data=dat1,alpha=0.4,size=line.tick,color=line.color,curvature=0.4)
        }else{
          pt <- pt + ggplot2::geom_segment(ggplot2::aes(x=PC1_TST,y=PC2_TST,xend=PC1_TRN,yend=PC2_TRN),
                        data=dat1,alpha=0.4,size=line.tick,color=line.color)
        }
      }
    }

    pt <- pt  +
      ggplot2::geom_point(data=dat[dat$set==sets[1],],ggplot2::aes(shape=group,fill=set),color=point.color,size=set.size[1]) +
      ggplot2::geom_point(data=dat[dat$set==sets[2],],ggplot2::aes(shape=group,fill=set),color=point.color,size=set.size[2]) +
      ggplot2::theme_bw() + theme0
  }else{
      dat2 <- c()
      for(i in 1:length(tst))
      {
        indexTRN <- which(abs(B[which(fm$tst == tst[i]), ]) > 0)
        if(length(indexTRN) > 0)
        {
          tmp <- dat[subsetG[fm$trn], ][-indexTRN,]
          tmp$set <- sets[3]
          tmp <- rbind(dat[subsetG[fm$trn], ][indexTRN,], tmp, dat[subsetG[tst], ][i,])
          dat2 <- rbind(dat2,data.frame(tmp, ind = i))
        }
      }

      pt <- ggplot2::ggplot(dat2,ggplot2::aes(x=PC1,y=PC2)) + ggplot2::facet_wrap(~ind) +
             ggplot2::geom_point(data=dat2[dat2$set==sets[3],],ggplot2::aes(fill=set,shape=group),color=point.color,size=set.size[3]) +
             ggplot2::geom_point(data=dat2[dat2$set==sets[2],],ggplot2::aes(fill=set,shape=group),color=point.color,size=set.size[2]) +
             ggplot2::geom_point(data=dat2[dat2$set==sets[1],],ggplot2::aes(fill=set,shape=group),color=point.color,size=set.size[1]) +
             ggplot2::theme_bw() + theme0

  }

  pt <- pt + ggplot2::labs(title=title0, x=labelsPC[1],y=labelsPC[2]) +
    ggplot2::scale_shape_manual(values = group.shape,
              guide=ggplot2::guide_legend(override.aes=list(size=2,fill="white"))) +
    ggplot2::scale_fill_manual(values = set.color,
              guide=ggplot2::guide_legend(override.aes=list(shape=21,size=2)))

 if(!flagGp) pt <- pt + ggplot2::guides(shape="none")

  pt
}


#====================================================================
# Plot the coefficients path in a penalized regression (user-level)
#====================================================================
# Z=NULL; K=NULL; tst=NULL; title=NULL; maxCor=0.85
plotPath <- function(fm, Z=NULL, K=NULL, tst=NULL, title=NULL, maxCor=0.85)
{
  k <- NULL
  flagKinship <- FALSE
  if(!inherits(fm,c("LASSO","SSI"))) stop("Object 'fm' is not of the class 'LASSO' or 'SSI'")

  theme0 <- ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.background = ggplot2::element_rect(fill = "gray95"),
    legend.box.spacing = ggplot2::unit(0.4, "lines"),
    legend.key.height= ggplot2::unit(3,"line"),
    legend.key.width = ggplot2::unit(0.8, "lines")
  )

  if(inherits(fm,"SSI"))
  {
    if(!is.null(K)){
      flagKinship <- TRUE
      if(is.character(K)){
        K <- readBinary(K)
      }
      if(!is.null(Z))
      {
        if(length(dim(Z)) != 2) stop("Object 'Z' must be a matrix")
        K <- float::tcrossprod(Z,float::tcrossprod(Z,K))  # G = ZKZ'
      }
      if(length(dim(K))!=2 | (length(K) != length(fm$y)^2))
        stop("Product Z %*% K %*% t(Z) must be a squared matrix with number of rows (and columns) equal to the number of elements in 'y'")
    }

    beta <- coef.SSI(fm)
    if(!is.null(tst)){
      if(any(!tst %in% fm$tst)) stop("Some elements in 'tst' vector are not contained in 'fm$tst'")
      indexTST <- which(fm$tst %in% tst)
    }else indexTST <- seq_along(fm$tst)
    beta <- beta[indexTST]
    lambda <- apply(fm$lambda,2,mean)
    df <- apply(fm$df,2,mean)

  }else{
    beta <- fm$beta
    lambda <- fm$lambda
    df <- fm$df
  }

  nDF <- length(df)
  if(min(lambda) < .Machine$double.eps*1000)  lambda[which.min(lambda)] <- min(lambda[lambda>0])/2
  if(nDF==1) stop("Coefficients path plot can not be generated for 'nLambda=1'")

  if(inherits(fm,"SSI"))
  {
    dat <- c()
    trim <- length(fm$trn)*length(beta) > 20000
    for(i in seq_along(beta))
    {
      b0 <- as.matrix(beta[[i]])
      if(trim){
        indexOK <- getIndexCorrelated(t(b0),maxCor)
      }else indexOK <- seq_along(fm$trn)

      tmp <- matrix(NA,nrow=1,ncol=length(indexOK))
      if(!is.null(K)){
        tmp <- K[fm$tst[indexTST[i]],fm$trn[indexOK],drop=FALSE]
        if(float::storage.mode(tmp)=='float32') tmp <- float::dbl(tmp)
      }
      dimnames(tmp) <- list(fm$tst[indexTST[i]],fm$trn[indexOK])
      tmp <- reshape2::melt(tmp)
      colnames(tmp) <- c("tst_i","trn_i","value")
      tmp <- tmp[rep(1:nrow(tmp),nDF),]

      df0 <- rep(df,each=length(indexOK))
      lambda0 <- rep(lambda,each=length(indexOK))
      b0 <- as.vector(b0[indexOK,])
      id <- factor(tmp$tst_i):factor(tmp$trn_i)
      dat <- rbind(dat,data.frame(df=df0,lambda=lambda0,beta=float::dbl(b0),k=tmp$value,id=id))
    }

  }else{
    id <- factor(rep(seq(nrow(beta)),ncol(beta)))
    dat <- data.frame(df=rep(df,each=nrow(beta)),lambda=rep(lambda,each=nrow(beta)),beta=as.vector(beta),id=id)
  }

  # Labels and breaks for the DF axis
  ax2 <- getSecondAxis(lambda,df)
  brks0 <- ax2$breaks
  labs0 <- ax2$labels

  title0 <- bquote("Coefficients path. "*.(fm$name))
  if(!is.null(title)) title0 <- title

  if(flagKinship)
  {
    pt <- ggplot2::ggplot(dat,ggplot2::aes(-log(lambda),beta,color=k,group=id)) +
      viridis::scale_color_viridis() + ggplot2::geom_line() + ggplot2::theme_bw() + theme0 +
      ggplot2::labs(title=title0,y=expression(beta),x=expression("-log("*lambda*")"))
  }else{
    pt <- ggplot2::ggplot(dat,ggplot2::aes(-log(lambda),beta,color=id,group=id))+
      ggplot2::geom_line() + ggplot2::theme_bw() + theme0 +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(title=title0,y=expression(beta),x=expression("-log("*lambda*")"))
  }

  if(length(brks0)>3){
    pt <- pt + ggplot2::scale_x_continuous(sec.axis=ggplot2::sec_axis(~.+0,"Number of predictors",breaks=brks0,labels=labs0))
  }
  pt
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
  |    Sparse Family and Selection Index. Version 1.0.0 (Sep 30, 2021)    |
  |    Type 'citation('SFSI')' to know how to cite SFSI                   |
  |    Type 'help(package='SFSI',help_type='html')' to see help           |
  |    Type 'browseVignettes('SFSI')' to see documentation                |
  |    Type 'demo(package='SFSI')' to see demos                           |
  |                                                                       |
  |=======================================================================|
  ")
}
