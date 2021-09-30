
#====================================================================
#====================================================================
coef.SSI <- function(object,...,df=NULL,tst=NULL)
{
  if(!is.null(df)){
      if( df < 0 | df > length(object$trn) )
        stop("Parameter 'df' must be greater than zero and no greater than the number of elements in the training set")
  }

  # Which index is the closest to DF (across all TST individuals)
  which.df <- which.min(abs(apply(object$df,2,mean)-df))

  if(is.null(tst)) tst <- object$tst
  if(sum(!tst %in%object$tst) > 0) stop("Some 'tst' elements were not found in the provided object")

  posTSTind <- NULL
  if(!is.null(object$subset))  posTSTind <- c(0,cumsum(object$subset_size))

  filename <- basename(object$file_beta)
  isSingleFile <- all(substr(filename,nchar(filename)-5,nchar(filename)) == "_B.bin")
  isSingleFile <- isSingleFile & length(filename)==1  # When saved using 'saveAt'
  posBetaInd <- c(0,cumsum(rep(ncol(object$df),length(object$tst))))

  BETA <- vector("list",length(tst))
  for(i in seq_along(tst))
  {
    j <- which(object$tst == tst[i])

    if(is.null(object$subset) & !isSingleFile)
    {
      indexCol <- NULL
      if(!is.null(df))   indexCol <- which.df
      tmp <- readBinary(paste0(object$file_beta,j,".bin"),indexCol=indexCol,verbose=FALSE)
    }else{
      if(isSingleFile){
          indexCol <- (posBetaInd[j]+1):posBetaInd[j+1]
          if(!is.null(df))  indexCol <- posBetaInd[j] + which.df

          tmp <- readBinary(object$file_beta,indexCol=indexCol,verbose=FALSE)
      }else{
        # For the case with subset
        f <- which((posTSTind[-length(posTSTind)] < j) & (j <= posTSTind[-1]))
        posBetaInd0 <- c(0,cumsum(rep(ncol(object$df),object$subset_size[f])))

        j2 <- j - posTSTind[f]
        indexCol <- (posBetaInd0[j2]+1):posBetaInd0[j2+1]
        if(!is.null(df))  indexCol <- posBetaInd[j2] + which.df

        tmp <- readBinary(object$file_beta[f],indexCol=indexCol,verbose=FALSE)
      }
      #cat("i=",i,"j=",j,"j2=",j2,"\n")
    }
    BETA[[i]] <- tmp # Matrix::Matrix(tmp, sparse=TRUE)
  }
  if(!is.null(df)) BETA <- do.call(cbind,BETA)

  BETA
}

#====================================================================
#====================================================================
fitted.LASSO <- function(object,...)
{
  args0 <- list(...)
  if(length(args0)==0) stop("A matrix of predictors must be provided")

  X <- args0[[1]]
  yHat <- X%*%as.matrix(object$beta)
  colnames(yHat) <- paste0("yHat",1:ncol(yHat))
  yHat
}

#====================================================================
#====================================================================
fitted.SSI <- function(object,...)
{
  args0 <- list(...)
  yTRN <- as.vector(object$y[object$trn]-object$Xb[object$trn])
  indexdrop <- is.na(yTRN)
  if(length(indexdrop)>0) yTRN[indexdrop] <- 0

  uHat <- do.call(rbind,lapply(seq_along(object$tst),function(i){
    float::crossprod(coef.SSI(object,tst=object$tst[i])[[1]],yTRN)[,1]
  }))
  dimnames(uHat) <- list(object$tst,paste0("SSI.",1:ncol(uHat)))

  return(uHat)
}

#====================================================================
#====================================================================
plot.SSI <- function(...,title=NULL,py=c("accuracy","MSE"))
{
    PC1 <- PC2 <- PC1_TST <- PC2_TST <- PC1_TRN <- PC2_TRN <- loglambda <- NULL
    k <- model <- y <- trn_tst <- NULL
    py <- match.arg(py)
    args0 <- list(...)

    object <- args0[unlist(lapply(args0,function(x)class(x)=="SSI"))]
    if(length(object)==0) stop("No object of the class 'SSI' was provided")

    # Treat repeated fm$name
    modelNames <- unlist(lapply(object,function(x)x$name))
    index <- table(modelNames)[table(modelNames)>1]
    if(length(index)>0){
      for(i in seq_along(index)){
        tmp <- which(modelNames == names(index[i]))
        for(j in seq_along(tmp)) object[[tmp[j]]]$name <- paste0(object[[tmp[j]]]$name,".",j)
      }
    }

    trn <- do.call(rbind,lapply(object,function(x)x$trn))
    tst <- do.call(rbind,lapply(object,function(x)x$tst))
    if(any(apply(trn,2,function(x)length(unique(x)))!=1)) stop("'Training' set is not same across all 'SSI' objects")
    if(any(apply(tst,2,function(x)length(unique(x)))!=1)) stop("'Testing' set is not same across all 'SSI' objects")
    trn <- as.vector(trn[1,])
    tst <- as.vector(tst[1,])

    nTST <- length(tst)
    nTRN <- length(trn)

    theme0 <- ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.background = ggplot2::element_rect(fill="gray95"),
      #legend.key = element_rect(fill="gray95"),
      legend.box.spacing = ggplot2::unit(0.4, "lines"),
      legend.justification = c(1,ifelse(py=="MSE",1,0)),
      legend.position=c(0.96,ifelse(py=="MSE",0.96,0.04)),
      legend.title = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(t=0,b=0.25,l=0.25,r=0.25,unit='line')
    )

      dat <- c()
      meanopt <- c()
      for(j in 1:length(object))
      {
          fm0 <- object[[j]]
          tmp <- summary.SSI(fm0)
          names(tmp[[py]]) <- paste0("SSI",1:length(tmp[[py]]))

          tt1 <- data.frame(SSI=names(tmp[[py]]),y=tmp[[py]],df=apply(fm0[['df']],2,mean),lambda=apply(fm0[['lambda']],2,mean))
          if(any(tt1$lambda<.Machine$double.eps)){
                tmp <- tt1$lambda[tt1$lambda>=.Machine$double.eps]
                tt1[tt1$lambda<.Machine$double.eps,'lambda'] <- ifelse(length(tmp)>0,min(tmp)/2,1E-6)
          }
          # tt1 <- data.frame(obj=j,model=fm0$name,method=fm0$method,tt1,loglambda=-log(tt1$lambda),stringsAsFactors=FALSE)
          tt1 <- data.frame(obj=j,model=fm0$name,tt1,loglambda=-log(tt1$lambda),stringsAsFactors=FALSE)
          dat <- rbind(dat,tt1)

          index <- ifelse(py=="MSE",which.min(tt1$y),which.max(tt1$y))
          meanopt <- rbind(meanopt,tt1[index,])
      }

      dat$model <- factor(as.character(dat$model))

      # Models with a single SSI (G-BLUP or SSI with a given value of lambda)
      index <- unlist(lapply(split(dat,dat$obj),function(x)length(unique(x$SSI))))
      dat2 <- dat[dat$obj %in% names(index[index==1]),]

      # Remove data from models with a single SSI point
      dat <- dat[!dat$obj %in% names(index[index==1]),]
      meanopt <- meanopt[!meanopt$obj %in% names(index[index==1]),]
      dat <- dat[!(is.na(dat$y) | is.na(dat$loglambda)),]   # Remove NA values
      labY <- ifelse(py=="accuracy",expression('cor(y,'*hat(y)*')'),py)
      labX <- expression("-log("*lambda*")")

      if(nrow(dat)==0 | nrow(meanopt)==0)  stop("The plot can not be generated with the provided data")

      title0 <- paste0("Testing set ",py)
      if(!is.null(title)) title0 <- title

      # Labels and breaks for the DF axis
      ax2 <- getSecondAxis(dat$lambda,dat$df)
      brks0 <- ax2$breaks
      labs0 <- ax2$labels

      pt <- ggplot2::ggplot(dat,ggplot2::aes(loglambda,y,group=model,color=model)) +
            ggplot2::geom_line(size=0.66) +
            ggplot2::geom_hline(data=dat2,ggplot2::aes(yintercept=y,color=model,group=model),size=0.66) +
            ggplot2::labs(title=title0,x=labX,y=labY) + ggplot2::theme_bw() + theme0 +
            #ggplot2::ylim(min(dat$y[dat$df>1]),max(dat$y)) +
            ggplot2::geom_vline(data=meanopt,ggplot2::aes(xintercept=loglambda),size=0.5,linetype="dotted",color="gray50")

      if(length(brks0)>3){
        pt <- pt + ggplot2::scale_x_continuous(sec.axis=ggplot2::sec_axis(~.+0,"Support set size",breaks=brks0,labels=labs0))
      }
      pt
}

#====================================================================
#====================================================================
plot.SSI_CV <- function(...,py=c("accuracy","MSE"), title=NULL,showFolds=FALSE)
{
    py <- match.arg(py)
    args0 <- list(...)
    obj_CV_fold <- negLogLambda <- CV <- fold <- y <- model <- NULL

    object <- args0[unlist(lapply(args0,function(x)class(x)=="SSI_CV"))]

    # Treat repeated fm$name
    modelNames <- unlist(lapply(object,function(x)unique(unlist(lapply(x,function(z)z$name)))))
    index <- table(modelNames)[table(modelNames)>1]
    if(length(index)>0){
      for(i in seq_along(index)){
        tmp <- which(modelNames == names(index[i]))
        for(j in seq_along(tmp)){
          for(k in 1:length(object[[tmp[j]]]))
          object[[tmp[j]]][[k]]$name <- paste0(object[[tmp[j]]][[k]]$name,".",j)
        }
      }
    }
    varNames <- c("y","df","lambda","negLogLambda")
    dat <- c()
    for(j in 1:length(object))
    {
        fm0 <- object[[j]]
        names(fm0) <- paste0("CV",1:length(fm0))
          rawdat <- do.call(rbind,lapply(fm0,function(x)reshape2::melt(x[[py]])))
          colnames(rawdat) <- c("fold","SSI","y")
          rawdat <- data.frame(CV=unlist(lapply(strsplit(rownames(rawdat),"\\."),function(x)x[1])),rawdat)
          rawdat$df <- do.call(rbind,lapply(fm0,function(x)reshape2::melt(x[['df']])))$value
          rawdat$lambda <- do.call(rbind,lapply(fm0,function(x)reshape2::melt(x[['lambda']])))$value
          rawdat$negLogLambda <- -log(rawdat$lambda)
          #rawdat <- data.frame(obj=j,model=fm0[[1]]$name,method=fm0[[1]]$method,rawdat,stringsAsFactors=FALSE)
          rawdat <- data.frame(obj=j,model=fm0[[1]]$name,rawdat,stringsAsFactors=FALSE)
        dat <- rbind(dat,rawdat)
    }

    # Average across folds
    avgdat <- do.call(rbind,lapply(split(dat,paste0(dat$obj,"_",dat$CV,"_",dat$SSI)),function(x){
        x[1,varNames] <- apply(as.matrix(x[,varNames]),2,mean,na.rm=TRUE)
        x$fold <- "mean"
        x[1,]
    }))

    # Average across partitions across folds
    overalldat <- do.call(rbind,lapply(split(avgdat,paste0(avgdat$obj,"_",avgdat$SSI)),function(x){
      x[1,varNames] <- apply(as.matrix(x[,varNames]),2,mean,na.rm=TRUE)
      x$CV <- "mean"
      x[1,]
    }))

    if(showFolds & length(object)==1){
      dat <- rbind(dat,avgdat,overalldat)
    }else dat <- rbind(avgdat,overalldat)

    dat$obj_CV_fold <- factor(paste0(dat$obj,"_",dat$CV,"_",dat$fold))
    dat$obj <- factor(as.character(dat$obj))

    # Optimum INDEX
    optdat <- do.call(rbind,lapply(split(overalldat,overalldat$obj),function(x){
      x[ifelse(py =="accuracy",which.max(x$y),which.min(x$y)),]
    }))
    optdat$obj_CV_fold <- factor(paste0(optdat$obj,"_",optdat$CV,"_",optdat$fold))

    # Models with a single SSI
    index <- unlist(lapply(split(dat,paste0(dat$obj,"_",dat$CV)),function(x)mean(table(x$fold))))
    dat2 <- dat[paste0(dat$obj,"_",dat$CV) %in% names(index[index == 1]),]

    # Remove data from models with a single SSI point
    dat <- dat[paste0(dat$obj,"_",dat$CV) %in% names(index[index > 1]),]
    optdat <- optdat[paste0(optdat$obj,"_",optdat$CV) %in% names(index[index > 1]),]

    dat <- dat[!is.na(dat$y) & !is.na(dat$negLogLambda),]   # Remove NA values
    labY <- ifelse(py=="accuracy",expression('cor(y,'*hat(y)*')'),py)
    labX <- expression("-log("*lambda*")")

    # Labels and breaks for the DF axis
    if(nrow(dat) == 0) stop("The plot cannot be generated with the provided data")
    ax2 <- getSecondAxis(dat$lambda,dat$df)
    brks0 <- ax2$breaks
    labs0 <- ax2$labels

    title0 <- paste0("Average cross-validated ",py)
    if(!is.null(title)) title0 <- title

    theme0 <- ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.background = ggplot2::element_rect(fill="gray95"),
        legend.box.spacing = ggplot2::unit(0.4, "lines"),
        legend.justification = c(1,ifelse(py=="MSE",1,0)),
        legend.position=c(0.97,ifelse(py=="MSE",0.97,0.03)),
        legend.title=ggplot2::element_blank(),
        legend.margin=ggplot2::margin(t=0,b=0.25,l=0.25,r=0.25,unit='line')
    )

    dat <- dat[dat$df > 0.99,]
    index1 <- which(dat$CV == "mean" & dat$fold == "mean")
    index2 <- which(dat$CV != "mean" & dat$fold == "mean")
    index3 <- which(dat$CV != "mean" & dat$fold != "mean")

    pt <- ggplot2::ggplot(dat,ggplot2::aes(group=obj_CV_fold)) +
             ggplot2::geom_hline(data=dat2[dat2$CV == "mean" & dat2$fold == "mean",],
                        ggplot2::aes(yintercept=y,color=model,group=obj_CV_fold),size=0.7) +
             ggplot2::labs(title=title0,x=labX,y=labY) + ggplot2::theme_bw() + theme0 +
             ggplot2::geom_vline(data=optdat,ggplot2::aes(xintercept=negLogLambda),size=0.5,linetype="dotted",color="gray50")

    if(showFolds){
      if(length(object)==1){
        pt <- pt + ggplot2::geom_line(data=dat[index3,],ggplot2::aes(negLogLambda,y,group=obj_CV_fold),color="gray70",size=0.3)
      }else cat("Results for individuals folds are not shown when plotting more than one model\n")
    }

    if(length(object)==1){   # Results from each CV
      pt <- pt + ggplot2::geom_line(data=dat[index2,],ggplot2::aes(negLogLambda,y,group=obj_CV_fold,color=CV),size=0.4) +
                 ggplot2::geom_line(data=dat[index1,],ggplot2::aes(negLogLambda,y,group=obj_CV_fold),color="gray5",size=0.7) +
                 ggplot2::theme(legend.position = "none")
    }else{
      pt <- pt + ggplot2::geom_line(data=dat[index1,],ggplot2::aes(negLogLambda,y,group=obj_CV_fold,color=model),size=0.7)
    }

    if(length(brks0)>3){
      pt <- pt + ggplot2::scale_x_continuous(sec.axis=ggplot2::sec_axis(~.+0,"Support set size",breaks=brks0,labels=labs0))
    }
    pt
}

#====================================================================
#====================================================================
summary.SSI_CV <- function(object, ...)
{
    args0 <- list(...)
    if(!inherits(object, "SSI_CV")) stop("The provided object is not from the class 'SSI_CV'")

    df <- do.call(rbind,lapply(object,function(x)apply(x$df,2,mean,na.rm=TRUE)))
    lambda <- do.call(rbind,lapply(object,function(x)apply(x$lambda,2,mean,na.rm=TRUE)))
    MSE <- do.call(rbind,lapply(object,function(x)apply(x$MSE,2,mean,na.rm=TRUE)))
    accuracy <- do.call(rbind,lapply(object,function(x)apply(x$accuracy,2,mean,na.rm=TRUE)))
    rownames(df) <- rownames(lambda) <- rownames(accuracy) <- rownames(MSE) <- paste0("CV",1:length(object))

    out <- list(df=df,lambda=lambda,accuracy=accuracy,MSE=MSE)

    # Detect maximum accuracy by partition (curve)
    index <- apply(accuracy,1,which.max)
    tmp <- lapply(out,function(x)unlist(lapply(1:nrow(x),function(z)x[z,index[z]])))
    optCOR <- do.call(cbind,tmp)

    # Detect minimum MSE by partition (curve)
    index <- apply(MSE,1,which.min)
    tmp <- lapply(out,function(x)unlist(lapply(1:nrow(x),function(z)x[z,index[z]])))
    optMSE <- do.call(cbind,tmp)

    ##  Maximum accuracy averaging curves
    index <- which.max(apply(accuracy,2,mean,na.rm=TRUE))
    avg.optCOR <- unlist(lapply(out,function(x)apply(x,2,mean,na.rm=TRUE)[index]))

    ##  Minimum MSE averaging curves
    index <- which.min(apply(MSE,2,mean,na.rm=TRUE))
    avg.optMSE <- unlist(lapply(out,function(x)apply(x,2,mean,na.rm=TRUE)[index]))

    rownames(optCOR) <- rownames(optMSE) <- paste0("CV",1:length(object))
    optCOR <- rbind(optCOR,mean=avg.optCOR)
    optMSE <- rbind(optMSE,mean=avg.optMSE)

    do.call(c, list(as.list(out), optCOR=list(optCOR), optMSE=list(optMSE)))
}

#====================================================================
#====================================================================
summary.SSI <- function(object,...)
{
    args0 <- list(...)
    if(!inherits(object, "SSI")) stop("The provided object is not from the class 'SSI'")

    tst <- object$tst
    y <- as.vector(object$y)

    df <- apply(object$df,2,mean)
    lambda <- apply(object$lambda,2,mean)
    uHat <- fitted.SSI(object)
    #if(object$isFloat) uHat <- float::dbl(uHat)

    accuracy <- suppressWarnings(drop(stats::cor(y[tst],uHat,use="pairwise.complete.obs")))
    MSE <- suppressWarnings(apply((y[tst]-uHat)^2,2,sum,na.rm=TRUE)/length(tst))

    out <- data.frame(accuracy=accuracy,MSE=MSE,df=df,lambda=lambda)

    # Detect maximum accuracy
    index <- which.max(out$accuracy)
    if(length(index)==0)
    {
      optCOR <- out[1, ,drop=FALSE]
      #if(nrow(out)>1) optCOR[1,] <- NA
    }else optCOR <- out[index, ,drop=FALSE]

    # Detect minimum MSE
    index <- which.min(out$MSE)
    if(length(index)==0)
    {
      optMSE <- out[1, ,drop=FALSE]
      #if(nrow(out)>1) optMSE[1,] <- NA
    }else optMSE <- out[index, ,drop=FALSE]

    rownames(optCOR) <- rownames(optMSE) <- NULL

    tmp <- as.list(out)
    tmp <- lapply(tmp,function(x){names(x)=rownames(out);x})
    do.call(c, list(tmp, optCOR=list(optCOR), optMSE=list(optMSE)))
}
