
#====================================================================
#====================================================================
coef.LASSO <- function(object, ...)
{
  args0 <- list(...)

  iy <- ilambda <- nsup <-  NULL
  if(length(dim(object$nsup)) == 2L){
    nsup0 <- colMeans(object$nsup)
  }else{
    nlambda <- unique(unlist(lapply(object$nsup,length)))
    if(length(nlambda)==1L){
      if(object$q == 1L){
        nsup0 <- object$nsup
      }else{
        nsup0 <- colMeans(do.call(rbind, object$nsup))
      }
    }else{ # This is the case of a LAR-LASSO
      nsup0 <- object$nsup
    }
  }

  if("iy" %in% names(args0)){
    iy <- args0$iy
  }
  if(is.null(iy)){
    iy <- seq(object$q)
  }
  if( any(!iy %in% seq(object$q)) ){
    stop("All elements in 'iy' must be between 1 and ",object$q)
  }

  if("ilambda" %in% names(args0)){
    if(length(args0$ilambda)>1L){
      message(" Only the first element of 'ilambda' is considered")
    }
    ilambda <- rep(args0$ilambda[1],length(iy))
    if( ilambda[1] < 1L | ilambda[1] > min(object$nlambda) ){
      stop("Parameter 'ilambda' must be between 1 and ",min(object$nlambda))
    }
  }

  if("nsup" %in% names(args0)){
    if(is.null(ilambda)){
      if(length(args0$nsup)>1L){
        message(" Only the first element of 'nsup' is considered")
      }
      nsup <- args0$nsup[1]
      if( nsup < 1L | nsup > object$p ){
        stop("Parameter 'nsup' must be between 1 and ", object$p)
      }
      if(is.list(nsup0)){ # LAR-LASSO case
        ilambda <- lapply(nsup0,function(x)which.min(abs(x-nsup)))
        message(" A different column was selected for each 'iy' for the case of a LAR-LASSO")
      }else{
        ilambda <- rep(which.min(abs(nsup0-nsup)),length(iy))
      }
    }else{
      message(" Parameter 'nsup' is ignored when 'ilambda' is provided")
    }
  }

  #message("'nsup' is NULL: ",is.null(nsup))
  #message("'ilambda' is NULL: ",is.null(ilambda))
  #message("'ilambda' length: ",length(ilambda))
  #message("ilambda=: ",ilambda)

  if(is.null(object$file_beta)){
    if(is.null(object$beta)){
      stop("No regression coefficients were found for the input object")
    }
    if(object$q == 1L){
      BETA <- object$beta
      if(!is.null(ilambda)){
        BETA <- BETA[,ilambda, drop=FALSE]
      }
    }else{
      BETA <- object$beta[iy]
      if(!is.null(ilambda)){
        BETA <- lapply(1:length(BETA), function(k) BETA[[k]][,ilambda[k], drop=FALSE])
      }
    }

  }else{
    BETA <- vector("list",length(iy))
    for(k in 1:length(iy)){
      tmp <- paste0("i_",object$fileID[iy[k]],".bin")
      indexcol <- NULL
      if(!is.null(ilambda)){
        indexcol <- ilambda[k]
      }
      BETA[[k]] <- readBinary(gsub("i_\\*.bin$",tmp,object$file_beta),
                              cols=indexcol, verbose=FALSE)
    }
  }

  if(is.list(BETA)){
    if(length(BETA) == 1L){
      BETA <- BETA[[1]]
    }else{
      if(!is.null(ilambda)){
        BETA <- t(do.call(cbind,BETA))
      }
    }
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
  }else{
    if('X' %in% names(args0)){
      X <- args0$X
    }else{
      if(length(args0) > 1L){
        message(" Only the second argument is considered as the matrix of predictors")
      }
      X <- args0[[1]]
    }
    if(length(dim(X)) == 2L){
      X <- as.matrix(X)
    }else{
      X <- matrix(X, nrow=1L)
    }

    yHat <- lapply(seq(object$q),function(i){
      tmp <- X%*%coef.LASSO(object, iy=i)
      colnames(tmp) <- paste0("yHat",1:ncol(tmp))
      tmp
    })
    if(nrow(X)==1L | object$q==1L){
      yHat <- do.call(rbind, yHat)
    }
  }
  return(yHat)
}

#====================================================================
#====================================================================

coef.SSI <- function(object, ...){
  coef.LASSO(object, ...)
}

#====================================================================
#====================================================================
fitted.SSI <- function(object, ...)
{
  args0 <- list(...)
  if(!inherits(object, "SSI")){
     stop("The input object is not of the class 'SSI'")
  }
  if("CV" %in% names(object)){
     stop("'fitted' method cannot be applied after cross-validation")
  }

  if(length(args0) == 0L){
    yHat <- object$yHat

  }else{
    if('y' %in% names(args0)){
      y0 <- as.vector(args0$y)
    }else{
      if(length(args0)>1L){
        message(" Only the second argument is considered as the response matrix")
      }
      y0 <- as.vector(args0[[1]])
    }

    if(length(y0) != (object$n * object$ntraits)){
      stop("Length of the response matrix must be equal to length(object$y)")
    }
    if(any(is.na(y0[object$trn]))){
      stop("All entries in y[trn] must be non-NA")
    }

    if(is.null(object$b)){
      yTRN <- as.vector(y0[object$trn])
    }else{
      yTRN <- as.vector(y0[object$trn] - object$Xb[object$trn])
    }

    u <- fitted.LASSO(object, yTRN)
    dimnames(u) <- list(object$tst, paste0("SSI.",1:ncol(u)))

    if(is.null(object$b)){
      yHat <- u[]
    }else{
      yHat <- sweep(u, 1L, object$Xb[object$tst], FUN="+")
    }
  }

  return(yHat)
}

#====================================================================
#====================================================================
summary.SSI <- function(object, ...)
{
  args0 <- list(...)

  if(!inherits(object, "SSI")){
    stop("The input object is not of the class 'SSI'")
  }

  map <- map_trn <- map_tst <- nsup_trait <- NULL

  nTST <- length(object$tst)
  nTRN <- length(object$trn)

  flag_accuracy <- TRUE
  if(length(object$CV) > 0L)
  {
    nfolds <- object$nfolds
    nCV <- object$nCV

    # Average across all folds
    nsup <- lapply(object$CV,function(x)Reduce("+",x$nsup)/nfolds)
    lambda <- lapply(object$CV,function(x)Reduce("+",x$lambda)/nfolds)
    accuracy <- lapply(object$CV,function(x)Reduce("+",x$accuracy)/nfolds)
    MSE <- lapply(object$CV,function(x)Reduce("+",x$MSE)/nfolds)

    # Average across all CV repetitions
    nsup <- Reduce("+",nsup)/nCV
    lambda <- Reduce("+",lambda)/nCV
    accuracy <- Reduce("+",accuracy)/nCV
    MSE <- Reduce("+",MSE)/nCV

    if(object$ntraits > 1L){
      names_nsup <- paste0("nsup_",1:object$ntraits)
      nsup_trait <- lapply(object$CV,function(x){
        Reduce("+",lapply(x$nsup_trait,function(z)z[,names_nsup]))/nfolds
      })
      nsup_trait <- Reduce("+",nsup_trait)/nCV
      tmp <- object$CV[[1]]$nsup_trait[[1]][,c("SSI","trait")]
      nsup_trait <- data.frame(tmp, nsup_trait)

      map <- map_set(object$n, object$ntraits, x=object$trn, y=NULL,
                   xlab="trn", ylab="tst")
      map_trn <- map[object$trn,]

      tt <- factor(as.character(map_trn$j), levels=seq(object$ntraits))
      nTRN <- c("Across"=nTRN, table(tt))
    }

  }else{
    if(length(args0) == 0L){
      y <- object$y
      yHat <- object$yHat # Should be equal to fitted.SSI(object)

    }else{
      if('y' %in% names(args0)){
        y <- args0$y
      }else{
        if(length(args0)>1L){
          message(" Only the second argument is considered as the response matrix")
        }
        y <- args0[[1]]
      }
      yHat <- fitted.SSI(object, y)
    }
    if(length(dim(y)) == 2L){
      stopifnot(nrow(y) == object$n)
      stopifnot(ncol(y) == object$ntraits)
    }else{
      if(object$ntraits > 1L){
        stop("Response matrix 'y' should contain ",object$n," rows and ",object$ntraits," columns")
      }
    }
    y0 <- as.vector(y)

    if(sum(is.na(y0[object$tst])) > 0L){
       message(" Some testing entries in the response matrix 'y' are NA.\n",
               " Provide a full 'y' matrix to compute accuracy in testing data")
       flag_accuracy <- FALSE
    }

    yTST <- y0[object$tst]
    tmp <- list(colnames(yHat), "Across")
    nsup <- matrix(colMeans(object$nsup),ncol=1,dimnames=tmp)
    lambda <- matrix(colMeans(object$lambda),ncol=1,dimnames=tmp)
    accuracy <- matrix(suppressWarnings(stats::cor(yHat,yTST)),ncol=1,dimnames=tmp)
    MSE <- matrix(suppressWarnings(apply(sweep(yHat,1L,yTST,FUN="-")^2,2,mean)),
                  ncol=1,dimnames=tmp)

    nsup1 <- lambda1 <- accuracy1 <- MSE1 <- NULL
    if(object$ntraits > 1L){  # Calculate MSE and accuracy within response variable
      map <- map_set(object$n, object$ntraits, x=object$trn, y=object$tst,
                     xlab="trn", ylab="tst")
      map_trn <- map[object$trn,]
      map_tst <- map[object$tst,]

      nsup1 <- lambda1 <- accuracy1 <- MSE1 <- matrix(NA,
                                                    nrow=object$nlambda,
                                                    ncol=object$ntraits)
      colnames(nsup1) <- colnames(lambda1) <- seq(object$ntraits)
      colnames(accuracy1) <- colnames(MSE1) <- seq(object$ntraits)

      for(j in 1:object$ntraits){
        index <- which(map_tst$j == j)
        map0 <- map_tst[index,]
        nsup1[,j] <- colMeans(object$nsup[index,,drop=F])
        lambda1[,j] <- colMeans(object$lambda[index,,drop=F])

        yHat0 <- yHat[map0$index_set,,drop=F]
        accuracy1[,j] <- drop(suppressWarnings(stats::cor(yHat0,y0[map0$index])))
        MSE1[,j] <- suppressWarnings(apply(sweep(yHat0,1L,y0[map0$index],FUN="-")^2,2,mean))
      }

      # Get nsup_trait: nsup for each trait in tst corresponding to each trait in trn
      nsup_trait <- get_summary_nsup(object, map=map)

      nsup <- cbind(nsup, nsup1)
      lambda <- cbind(lambda, lambda1)
      accuracy <- cbind(accuracy, accuracy1)
      MSE <- cbind(MSE, MSE1)

      tt <- factor(as.character(map_trn$j), levels=seq(object$ntraits))
      nTRN <- c("Across"=nTRN, table(tt))
      tt <- factor(as.character(map_tst$j), levels=seq(object$ntraits))
      nTST <- c("Across"=nTST, table(tt))
    }
  }

  index <- which(colnames(accuracy)=="Across")
  out <- data.frame(accuracy=accuracy[,index], MSE=MSE[,index],
                    nsup=nsup[,index], lambda=lambda[,index])

  # Detect maximum accuracy
  index <- which.max(out$accuracy)
  if(length(index) == 0L){
    optCOR <- out[1, ,drop=FALSE]
    if(nrow(out)>1) optCOR[1,] <- NA
  }else{
    optCOR <- out[index, ,drop=FALSE]
  }

  # Detect minimum MSE
  index <- which.min(out$MSE)
  if(length(index)==0){
    optMSE <- out[1, ,drop=FALSE]
    if(nrow(out)>1) optMSE[1,] <- NA
  }else{
    optMSE <- out[index, ,drop=FALSE]
  }

  optMSE <- as.matrix(optMSE)[1,]
  optCOR <- as.matrix(optCOR)[1,]

  out <- list(accuracy=accuracy, MSE=MSE, nsup=nsup, lambda=lambda)
  if(object$ntraits > 1L){
    out$nsup_trait <- nsup_trait
  }

  if(!flag_accuracy){
    out$accuracy <- NULL
    out$MSE <- NULL
  }

  tmp <- list(n=object$n, q=object$q,
              ntraits=object$ntraits,
              nTST=list(nTST), nTRN=list(nTRN),
              out,
              optCOR=list(optCOR),
              optMSE=list(optMSE))

  return(do.call(c, tmp))
}

#====================================================================
#====================================================================
# x.stat="nsup"; y.stat="accuracy"; nbreaks.x=7
plot.SSI <- function(..., x.stat = c("nsup","lambda"),
                     y.stat = c("accuracy","MSE"),
                     nbreaks.x = 7)
{
    x <- y <- name <- obj <- lambda <- NULL

    args0 <- list(...)
    x.stat <- match.arg(x.stat)
    y.stat <- match.arg(y.stat)

    #if(!inherits(x, "SSI")){
    #  stop("Input 'x' is not of the class 'SSI'")
    #}
    object <- list()
    #if(!missing(y)){
    #  if(inherits(y, "SSI")) object[[length(object)+1]] <- y
    #}
    if(length(args0) > 0L){
      for(i in 1:length(args0)){
        if(inherits(args0[[i]], "SSI")) object[[length(object)+1]] <- args0[[i]]
      }
    }
    if(length(object) == 0L){
       stop("No input object of the class 'SSI' was provided")
    }

    xlab <- ifelse(x.stat=="nsup","Support set size",expression(paste("-log(",lambda,")")))
    ylab <- capitalize(y.stat)
    lwd <- ifelse("lwd" %in% names(args0), args0$lwd, 0.65)
    if("xlab" %in% names(args0)) xlab <- args0$xlab
    if("ylab" %in% names(args0)) ylab <- args0$ylab

    # Treat repeated fm$name
    objectNames <- unlist(lapply(1:length(object),function(k) object[[k]]$name))
    index <- table(objectNames)[table(objectNames)>1L]
    if(length(index) > 0L){
      for(i in seq_along(index)){
        tmp <- which(objectNames == names(index[i]))
        for(j in seq_along(tmp)) objectNames[[tmp[j]]] <- paste0(objectNames[[tmp[j]]],"-",j)
      }
    }

    nTRN <- unlist(lapply(object,function(x)length(x$trn)))
    nTST <- unlist(lapply(object,function(x)length(x$tst)))
    isCV <- unlist(lapply(object,function(x)length(x$CV)>0))
    isTRN_TST <- unlist(lapply(object,function(x){
     ifelse(length(x$trn)==length(x$tst), all(x$trn==x$tst), FALSE)
    }))

    if(length(unique(nTRN))>1L){
       stop("Training set size is not same across all 'SSI' class objects")
    }
    if(length(unique(nTST))>1L){
       stop("Testing set size is not same across all 'SSI' class objects")
    }
    if(length(table(isCV))>1L){
       stop("All 'SSI' class objects must be of the same type:\n",
            "  either a trn-tst prediction or a cross-validation")
    }
    isCV <- isCV[1]

    if(isCV){
      main <- bquote("SSI CV ("*n[trn]==.(nTRN[1])*")")
    }else{
      if(isTRN_TST[1]){
        main <- bquote("SSI ("*n[trn]==.(nTST[1])*")")
      }else{
        main <- bquote("SSI ("*n[tst]==.(nTST[1])*")")
      }
    }
    if("main" %in% names(args0)) main <- args0$main

    theme0 <- theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.background = element_rect(fill="gray95"),
      legend.box.spacing = unit(0.4, "lines"),
      legend.key.size = unit(0.85, "lines"),
      legend.text = element_text(size=8),
      legend.justification = c(1,ifelse(tolower(y.stat)=="mse",1,0)),
      legend.position=c(0.99,ifelse(tolower(y.stat)=="mse",0.99,0.01)),
      legend.title = element_blank(),
      legend.margin = margin(t=-0.2,r=0.2,b=0.2,l=0.2,unit='line'),
      strip.text.x = element_text(size=8.5, margin=margin(t=1.5,b=1.5))
    )

    eps <- .Machine$double.eps
    dat <- data.frame(matrix(nrow=0,ncol=9))
    nlambda <- ntraits <- rep(NA, length(object))
    for(k in 1:length(object))
    {
      fm0 <- object[[k]]
      ntraits[k] <- fm0$ntraits
      nlambda[k] <- fm0$nlambda
      ss <- summary.SSI(fm0)

      if(nlambda[k] > 1L & !is.null(ss$accuracy) & !is.null(ss$MSE))
      {
        tt <- reshape2::melt(ss$lambda)
        colnames(tt) <- c("SSI","trait","lambda")
        tt$x <- -log(tt$lambda)
        tmp <- ss[[which(tolower(names(ss)) == tolower(y.stat))]]
        tt$y <- reshape2::melt(tmp)$value
        tt$nsup <- reshape2::melt(ss$nsup)$value
        tt <- tt[as.character(tt$trait) == "Across",]

        if(any(tt$lambda < eps)){
          tmp <- tt$lambda[tt$lambda >= eps]
          tt[tt$lambda < eps,'lambda'] <- ifelse(length(tmp)>0,min(tmp)/10,1E-6)
        }

        if(isCV){
          tt$n0 <- ss$nTRN[as.character(tt$trait)]
        }else{
          tt$n0 <- ss$nTST[as.character(tt$trait)]
        }

        tt <- data.frame(object=k,name=objectNames[k],tt,stringsAsFactors=FALSE)
        dat <- rbind(dat,tt)
      }
    }

    dat$name <- factor(as.character(dat$name))

    if(any(nlambda == 1L)){
      message(" Object(s) ",paste(which(nlambda==1L),collapse=",")," contain a single SSI point.",
              " They are excluded from the plot")
    }

    if(nrow(dat) == 0){
       stop("The plot can not be generated with the provided data")
    }

    dat <- dat[!is.na(dat$y),]   # Remove NA values
    dat$trait <- factor(as.character(dat$trait))

    if("ylim" %in% names(args0)){
       ylim <- args0$ylim
    }else{
       ylim <- c(NA, NA)
    }

    breaksx <- labelsx <- NULL
    if(x.stat=="nsup")
    {
      if("xlim" %in% names(args0)){
         xlim <- args0$xlim
      }else{
         xlim <- c(1,max(dat$nsup, na.rm=TRUE))
      }

      index <- dat$nsup >= xlim[1] & dat$nsup <= xlim[2]
      dat <- dat[index,]

      # Labels and breaks for the nsup axis
      breaks0 <- get_breaks(x=dat$x, y=dat$nsup, nbreaks.x=nbreaks.x, ymin=xlim[1])
      breaksx <- breaks0$breaks.x
      labelsx <- breaks0$breaks.y

    }else{
      if("xlim" %in% names(args0)){
        xlim <- args0$xlim
      }else{
        tmp <- dat[dat$nsup>=1,]
        tmp <- tmp[abs(tmp$nsup-min(tmp$nsup))<1E-8,,drop=F]
        xlim <- c(mean(tmp$x, na.rm=T), max(dat$x))
      }

      dat <- dat[dat$x >= xlim[1] & dat$x <= xlim[2],]
    }

    dat2 <- do.call(rbind,lapply(split(dat, paste(dat$trait,dat$object)),function(x){
      x[ifelse(tolower(y.stat)=="mse",which.min(x$y),which.max(x$y)),]
    }))

    pp <- ggplot(dat, aes(x,y,group=object,color=name)) +
          geom_line(size=lwd) +
          labs(title=main, x=xlab, y=ylab) +
          theme_bw() + theme0 +
          geom_vline(data=dat2, aes(xintercept=x),
                      size=0.5,linetype="dotted",color="gray50") +
          scale_y_continuous(limits=ylim)

    if(x.stat=="nsup"){
       pp <- pp + scale_x_continuous(breaks=breaksx, labels=round(labelsx))
    }else{
       pp <- pp +
        scale_x_continuous(breaks=scales::extended_breaks(n=nbreaks.x), limits=xlim)
    }

    return(pp)
}
