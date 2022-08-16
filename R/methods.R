
#====================================================================
#====================================================================
coef.SSI <- function(object, ..., df = NULL, i = NULL)
{
  if("CV" %in% names(object)) stop("'coef' method cannot be applied after cross-validation")
  if(!is.null(df)){
      if( df < 0 | df > length(object$trn) ){
        stop("Parameter 'df' must be between 0 <= df <= length(trn)")
      }
  }

  singleLambda <- ncol(object$lambda) == 1L
  # Which index is the closest to DF (across all TST individuals)
  which.df <- which.min(abs(apply(object$df,2,mean)-df[1]))

  if(is.null(df)){
    if(singleLambda) which.df <- 1
  }else{
    if(!singleLambda) singleLambda <- TRUE
  }

  if(is.null(i)){
    i <- seq_along(object$tst)
  }else{
    if( max(i) > length(object$tst) ){
      stop("All elements in 'i' must be between 0 < i <= length(object$tst)")
    }
  }

  filename <- NULL
  if(!is.null(object$file_beta)){
     filename <- basename(object$file_beta)
  }
  isSingleFile <- length(filename)==1L

  if(is.null(filename)){
    if(is.null(object$beta)){
      stop("No regression coefficients were found for the input object")
    }

    if(length(object$tst) == 1L){
      BETA <- object$beta
    }else{
      BETA <- object$beta[i]
    }

    if(singleLambda){
      BETA <- lapply(BETA, function(x) x[,which.df, drop=FALSE])
    }

  }else{
    BETA <- vector("list",length(i))
    for(j in 1:length(i))
    {
      if(isSingleFile){
        load(gsub("i_\\*.RData",paste0("i_",i[j],".RData"),object$file_beta))
        if(!is.null(df)){
             beta <- beta[, which.df, drop=FALSE]
        }
      }else{
        # For the case with subset
        f <- which(object$subset$index[,1]<= i[j] & object$subset$index[,2]>= i[j])
        i0 <- ifelse(f==1, i[j], i[j]-object$subset$index[f-1,2])
        load(gsub("i_\\*.RData",paste0("i_",i0,".RData"),object$file_beta[f]))
        if(!is.null(df))  beta <- beta[, which.df, drop=FALSE]

        #cat("i=",i,"j=",j,"j2=",j2,"\n")
      }
      BETA[[j]] <- beta
    }
  }

  if(length(BETA) == 1){
    BETA <- BETA[[1]]
  }else{
    if(singleLambda) BETA <- Matrix::t(do.call(cbind,BETA))
  }

  BETA
}

#====================================================================
#====================================================================
coef.LASSO <- function(object, ..., i = NULL)
{
  if(is.null(i)){
    i <- seq(object$q)
  }else{
    if( any(!i %in% seq(object$q)) ){
      stop("All elements in 'i' must be between 0 < i <= ncol(Gamma)")
    }
  }

  if(is.null(object$file_beta)){
    if(is.null(object$beta)){
      stop("No regression coefficients were found for the input object")
    }
    if(object$q ==1L){
      BETA <- object$beta
    }else{
      BETA <- object$beta[i]
    }
  }else{
    BETA <- vector("list",length(i))
    for(k in 1:length(i)){
      load(gsub("i_\\*.RData",paste0("i_",i[k],".RData"),object$file_beta))
      BETA[[k]] <- beta
    }
  }

  if(length(BETA) == 1){
    BETA <- BETA[[1]]
  }

  BETA
}

#====================================================================
#====================================================================
fitted.LASSO <- function(object, ...)
{
  args0 <- list(...)
  if(length(args0) == 0L){
    stop("A matrix of predictors must be provided")
  }
  if('X' %in% names(args0)){
    X <- args0$X
  }else{
    if(length(args0) > 1L){
      message("Only the second argument is considered as the matrix of predictors")
    }
    X <- args0[[1]]
  }
  if(length(dim(X)) != 2L){
    X <- matrix(X, nrow=1L)
  }

  yHat <- lapply(seq(object$q),function(i){
    tmp <- X%*%as.matrix(coef.LASSO(object, i=i))
    colnames(tmp) <- paste0("yHat",1:ncol(tmp))
    tmp
  })
  if(nrow(X)==1L | object$q==1L){
    yHat <- do.call(rbind, yHat)
  }
  yHat
}

#====================================================================
#====================================================================
fitted.SSI <- function(object, ...)
{
  args0 <- list(...)
  if("CV" %in% names(object)) stop("'fitted' method cannot be applied after cross-validation")

  if(length(args0) > 0){
    if('y' %in% names(args0)){
      y <- args0$y
    }else{
      if(length(args0)>1){
        message("Only the second argument is considered as the obvervations vector")
      }
      y <- args0[[1]]
    }

    if(length(y) != length(object$y)){
      stop("Length of the observations vector must be equal to length(object$y)")
    }

    if(any(is.na(y[object$trn]))){
      stop("All entries in y[trn] must be non-NA")
    }

    yTRN <- matrix(as.vector(y[object$trn] - object$Xb[object$trn]), nrow=1)

    u <- do.call(rbind,lapply(seq_along(object$tst),function(i){
      yTRN %*% as.matrix(coef.SSI(object, i=i))
    }))
    dimnames(u) <- list(object$tst, paste0("SSI.",1:ncol(u)))

  }else{
    u <- object$u
  }

  return(u)
}

#====================================================================
#====================================================================
plot.SSI <- function(..., py = c("accuracy","MSE"), nbreaks.x = 6)
{
    name <- y <- obj <- lambda <- NULL
    py <- match.arg(py)
    args0 <- list(...)

    xlab <- "Support set size";  ylab <- py
    main <- NULL
    lwd <- ifelse("lwd" %in% names(args0), args0$lwd, 0.65)
    if("main" %in% names(args0)) main <- args0$main
    if("xlab" %in% names(args0)) xlab <- args0$xlab
    if("ylab" %in% names(args0)) ylab <- args0$ylab
    if("y" %in% names(args0)){
      y <- args0$y
    }

    object <- args0[unlist(lapply(args0,function(x)inherits(x, "SSI")))]
    if(length(object) == 0L){
       stop("No object of the class 'SSI' was provided")
    }
    if(length(object)>1L & !is.null(y)){
      message("Same observation vector is used for all 'SSI' class objects")
    }

    # Treat repeated fm$name
    objectNames <- unlist(lapply(1:length(object),function(k) object[[k]]$name))
    index <- table(objectNames)[table(objectNames)>1]
    if(length(index) > 0L){
      for(i in seq_along(index)){
        tmp <- which(objectNames == names(index[i]))
        for(j in seq_along(tmp)) objectNames[[tmp[j]]] <- paste0(objectNames[[tmp[j]]],"-",j)
      }
    }

    trn <- unlist(lapply(object,function(x)x$trn))
    tmp <- unlist(lapply(object,function(x)all(trn %in% x$trn) & all(x$trn %in% trn)))
    if(any(!tmp)){
       message("'Training' set is not same across all 'SSI' class objects")
    }

    theme0 <- ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.background = ggplot2::element_rect(fill="gray95"),
      legend.box.spacing = ggplot2::unit(0.4, "lines"),
      legend.justification = c(1,ifelse(py=="MSE",1,0)),
      legend.position=c(0.99,ifelse(py=="MSE",0.99,0.01)),
      legend.title = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(t=-0.2,r=0.2,b=0.2,l=0.2,unit='line')
    )

    dat <- c()
    meanopt <- c()
    for(j in 1:length(object))
    {
        fm0 <- object[[j]]
        if(is.null(y)){
          ss <- summary.SSI(fm0)
        }else{
          ss <- summary.SSI(fm0, y=y)
        }
        names(ss[[py]]) <- paste0("SSI.",1:length(ss[[py]]))

        tt <- data.frame(SSI=names(ss[[py]]),y=ss[[py]],df=ss[['df']],lambda=ss[['lambda']])
        if(any(tt$lambda < .Machine$double.eps)){
          tmp <- tt$lambda[tt$lambda >= .Machine$double.eps]
          tt[tt$lambda < .Machine$double.eps,'lambda'] <- ifelse(length(tmp)>0,min(tmp)/10,1E-6)
        }

        tt <- data.frame(obj=j,name=objectNames[j],tt,stringsAsFactors=FALSE)
        dat <- rbind(dat,tt)

        index <- ifelse(py=="MSE",which.min(tt$y),which.max(tt$y))
        meanopt <- rbind(meanopt,tt[index,])
    }

    dat$name <- factor(as.character(dat$name))

    # Models with a single SSI (G-BLUP or SSI with a given value of lambda)
    index <- unlist(lapply(split(dat,dat$obj),function(x)length(unique(x$SSI))))
    dat2 <- dat[dat$obj %in% names(index[index==1]),]

    # Remove data from models with a single SSI point
    dat <- dat[!dat$obj %in% names(index[index==1]),]
    meanopt <- meanopt[!meanopt$obj %in% names(index[index==1]),]

    if(nrow(dat)==0 | nrow(meanopt)==0){
        stop("The plot can not be generated with the provided data")
    }

    dat <- dat[!is.na(dat$y),]   # Remove NA values

    if("xlim" %in% names(args0)){
       xlim <- args0$xlim
    }else{
       xlim <- c(1,max(dat$df, na.rm=TRUE))
    }
    dat <- dat[dat$df >= xlim[1] & dat$df <= xlim[2],]

    if("ylim" %in% names(args0)){
       ylim <- args0$ylim
    }else{
       ylim <- range(dat$y, na.rm=TRUE)
    }

    # Labels and breaks for the DF axis
    tmp <- get_breaks(dat$lambda, dat$df, nbreaks=nbreaks.x, ymin=xlim[1])
    breaks0 <- tmp$breaks.x
    labels0 <- round(tmp$breaks.y)
    labels2 <- sprintf('%.2f', breaks0)

    dat <- dat[dat$y >= ylim[1] & dat$y <= ylim[2],]

    pt <- ggplot2::ggplot(dat,ggplot2::aes(-log(lambda),y,group=obj,color=name)) +
          ggplot2::geom_line(size=lwd) + ggplot2::lims(y=ylim) +
          ggplot2::geom_hline(data=dat2,ggplot2::aes(yintercept=y,group=obj,color=name),size=lwd) +
          ggplot2::labs(title=main,x=xlab,y=ylab) + ggplot2::theme_bw() + theme0 +
          ggplot2::geom_vline(data=meanopt,ggplot2::aes(xintercept=-log(lambda)),size=0.5,linetype="dotted",color="gray50") +
          ggplot2::scale_x_continuous(breaks=breaks0,labels=labels0,
                      sec.axis=ggplot2::sec_axis(~.+0,expression(paste("-log(",lambda,")")),
                      breaks=breaks0, labels=labels2))

    pt
}

#====================================================================
#====================================================================
summary.SSI <- function(object, ...)
{
    args0 <- list(...)

    if(!inherits(object, "SSI")) stop("The input object is not of the class 'SSI'")

    if(length(object$CV) > 0L)
    {
      df <- apply(do.call(rbind,lapply(object$CV,function(x)x$df)),2,mean,na.rm=TRUE)
      lambda <- apply(do.call(rbind,lapply(object$CV,function(x)x$lambda)),2,mean,na.rm=TRUE)
      MSE <- apply(do.call(rbind,lapply(object$CV,function(x)x$MSE)),2,mean,na.rm=TRUE)
      accuracy <- apply(do.call(rbind,lapply(object$CV,function(x)x$accuracy)),2,mean,na.rm=TRUE)
    }else{
      if(length(args0) > 0){
        if('y' %in% names(args0)){
          y <- args0$y
        }else{
          if(length(args0)>1){
            message("Only the second argument is considered as the observations vector")
          }
          y <- as.vector(args0[[1]])
        }
        stopifnot(length(y) == length(object$y))
        u <- fitted.SSI(object, y)
      }else{
        y <- as.vector(object$y)
        u <- object$u
      }
      df <- apply(object$df, 2, mean)
      lambda <- apply(object$lambda, 2, mean)
      accuracy <- suppressWarnings(drop(stats::cor(y[object$tst],u,use="pairwise.complete.obs")))
      MSE <- suppressWarnings(apply((y[object$tst]-u)^2,2,sum,na.rm=TRUE)/length(object$tst))
    }

    out <- data.frame(accuracy=accuracy, MSE=MSE, df=df, lambda=lambda)

    # Detect maximum accuracy
    index <- which.max(out$accuracy)
    if(length(index) == 0L){
      optCOR <- out[1, ,drop=FALSE]
      #if(nrow(out)>1) optCOR[1,] <- NA
    }else{
       optCOR <- out[index, ,drop=FALSE]
    }

    # Detect minimum MSE
    index <- which.min(out$MSE)
    if(length(index)==0){
      optMSE <- out[1, ,drop=FALSE]
      #if(nrow(out)>1) optMSE[1,] <- NA
    }else{
       optMSE <- out[index, ,drop=FALSE]
    }

    optMSE <- as.matrix(optMSE)[1,]
    optCOR <- as.matrix(optCOR)[1,]

    tmp <- as.list(out)
    tmp <- lapply(tmp,function(x){names(x)=rownames(out);x})
    do.call(c, list(nTST=length(object$tst),
                    nTRN=length(object$trn), tmp,
                    optCOR=list(optCOR),
                    optMSE=list(optMSE)))
}
