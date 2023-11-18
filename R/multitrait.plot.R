#====================================================================
#====================================================================
# x.stat="nsup"; y.stat="accuracy"; nbreaks.x = 7; line.color="orange"
multitrait.plot <- function(object, x.stat = c("nsup","lambda"),
                            y.stat = c("accuracy","MSE"),
                            line.color = "orange",
                            nbreaks.x = 7, ...)
{
    x <- y <- nsupmin <- nsupmax <- nsup_trait <- label <- x1 <- NULL
    args0 <- list(...)
    x.stat <- match.arg(x.stat)
    y.stat <- match.arg(y.stat)

    if(!inherits(object, "SSI")){
      stop("Input 'object' is not of the class 'SSI'")
    }
    if(length(args0) > 0L){
      tmp <- unlist(lapply(args0, function(x)inherits(x, "SSI")))
      if(sum(tmp)>0){
        message(" More than one 'SSI' class objects were provided. Only the first one is plotted")
      }
    }

    xlab <- ifelse(x.stat=="nsup","Support set size",expression(paste("-log(",lambda,")")))
    ylab <- capitalize(y.stat)
    ylab2 <- "Proportion of support set"
    lwd <- ifelse("lwd" %in% names(args0), args0$lwd, 0.65)
    if("xlab" %in% names(args0)) xlab <- args0$xlab
    if("ylab" %in% names(args0)) ylab <- args0$ylab
    if("ylab2" %in% names(args0)) ylab2 <- args0$ylab2

    ntraits <- object$ntraits
    if(ntraits == 1L){
      stop("Input 'object' is not from a Multi-trait SSI")
    }
    nTRN <- length(object$trn)
    nTST <- length(object$tst)
    isCV <- length(object$CV) > 0
    isTRN_TST <- ifelse(length(object$trn)==length(object$tst),
                        all(object$trn==object$tst), FALSE)

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

    tmp <- ifelse(isCV,'trn',ifelse(isTRN_TST[1],'trn','tst'))
    facet_lab <- paste0("' * ' (' * n[",tmp,"] * ' = ' * ")
    theme0 <- theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.title = element_text(hjust=0.5),
      legend.background = element_rect(fill="gray95"),
      legend.key.size = unit(0.75, "lines"),
      legend.text = element_text(size=8),
      legend.justification = c(0,1), #ifelse(tolower(y.stat)=="mse",1,0)),
      legend.position="right",#c(0.99,ifelse(tolower(y.stat)=="mse",0.99,0.01)),
      legend.title = element_blank(),
      legend.margin = margin(t=-0.2,r=0.2,b=0.2,l=0.2, unit='lines'),
      legend.box.margin = margin(l=-10),
      strip.text.x = element_text(size=8.5, margin=margin(t=1.5,b=1.5)),
      axis.title.y.right = element_text(vjust = 1.6)
    )

    eps <- .Machine$double.eps
    dat <- data.frame(matrix(nrow=0,ncol=8))
    nlambda <- object$nlambda
    ss <- summary.SSI(object)

    # Get count of support set by trait
    names_nsup <- paste0("nsup_",1:ntraits)
    dat1 <- ss$nsup_trait
    dat1$nsup <- apply(dat1[,names_nsup],1,sum)
    dat1[,names_nsup] <- dat1[,names_nsup]/dat1$nsup

    # Data to plot
    if(nlambda > 1L & !is.null(ss$accuracy) & !is.null(ss$MSE))
    {
      tt <- reshape2::melt(ss$lambda)
      colnames(tt) <- c("SSI","trait","lambda")
      tt$x <- -log(tt$lambda)
      tmp <- ss[[which(tolower(names(ss)) == tolower(y.stat))]]
      tt$y <- reshape2::melt(tmp)$value
      tt$nsup <- reshape2::melt(ss$nsup)$value
      tt <- tt[as.character(tt$trait) != "Across",]

      if(any(tt$lambda < eps)){
        tmp <- tt$lambda[tt$lambda >= eps]
        tt[tt$lambda < eps,'lambda'] <- ifelse(length(tmp)>0,min(tmp)/10,1E-6)
      }

      if(isCV){
        tt$n0 <- ss$nTRN[as.character(tt$trait)]
      }else{
        tt$n0 <- ss$nTST[as.character(tt$trait)]
      }
      dat <- data.frame(tt,stringsAsFactors=FALSE)
    }
    stopifnot(paste(dat1$SSI,dat1$trait) == paste(dat$SSI,dat$trait))
    dat <- data.frame(dat, dat1[,names_nsup])

    if(nrow(dat) == 0){
       stop("The plot can not be generated with the provided data")
    }

    dat <- dat[!is.na(dat$y),]   # Remove NA values
    dat$trait2 <- factor(paste0("'Trait ",as.character(dat$trait), facet_lab, dat$n0,"L * ')'"))

    if("ylim" %in% names(args0)){
       ylim <- args0$ylim
    }else{
       ylim <- c(NA, NA)
    }

    if(x.stat=="nsup")
    {
      xd <- 5
      levels0 <- levels(dat$trait2)
      breaksx <- labelsx <- index <- x2 <- c()
      maxx <- 0
      for(k in 1:length(levels0))
      {
        index0 <- which(as.character(dat$trait2) == levels0[k])
        dat0 <- dat[index0,]
        if("xlim" %in% names(args0)){
           xlim <- args0$xlim
        }else{
           xlim <- c(1,max(dat0$nsup, na.rm=TRUE))
        }

        tmp <- dat0$nsup >= xlim[1] & dat0$nsup <= xlim[2]
        index0 <- index0[tmp]
        dat0 <- dat0[tmp,]

        # Labels and breaks for the nsup axis
        breaks0 <- get_breaks(x=dat0$x, y=dat0$nsup, nbreaks.x=nbreaks.x, ymin=xlim[1])

        dd <- maxx + (k-1)*xd
        maxx <- maxx + max(dat0$x)
        x2 <- c(x2, dd + dat0$x)
        breaksx <- c(breaksx, dd + breaks0$breaks.x)
        labelsx <- c(labelsx, breaks0$breaks.y)
        index <- c(index, index0)
      }
      dat <- dat[index,]
      dat$x <- x2
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


    tmp <- split(dat, dat$trait2)
    # Data for optimum SSI
    dat2 <- do.call(rbind,lapply(tmp,function(x){
      x[ifelse(tolower(y.stat)=="mse",which.min(x$y),which.max(x$y)),]
    }))

    ylim2 <- c(0,1)
    # Data for areas
    names0 <- colnames(dat)[!colnames(dat)%in%names_nsup]
    dat3 <- do.call(rbind,lapply(tmp,function(x){
      ylim1 <- range(x$y)
      b <- diff(ylim1)/diff(ylim2)
      a <- ylim1[1] - b*ylim2[1]
      tt <- b*as.matrix(x[,names_nsup])
      x2 <- c()
      for(i in 1:nrow(x)){
        cc0 <- a + c(0,cumsum(tt[i,]))
        cc1 <- cbind(nsupmin=cc0[1:length(names_nsup)],nsupmax=cc0[-1])
        x2 <- rbind(x2, data.frame(x[rep(i,length(names_nsup)),names0],
                                   nsup_trait=names_nsup,cc1))
      }
      rownames(x2) <- NULL
      x2
    }))
    dat3$nsup_trait <- gsub("nsup_","Trait ",dat3$nsup_trait)

    expand.x <- 0.025
    breaksy2 <- seq(0,1,length=ntraits+1)
    dat4 <- do.call(rbind,lapply(tmp,function(x){
      xlim <- range(x$x)
      ylim1 <- range(x$y)
      b <- diff(ylim1)/diff(ylim2)
      a <- ylim1[1] - b*ylim2[1]
      data.frame(x[1,c("trait","trait2"),drop=T],x=xlim[2],
                 x1=xlim[2] + 1*expand.x*diff(xlim),
                 x2=xlim[2] + 1*expand.x*diff(xlim),
                 y=a + b*breaksy2, label=sprintf('%.2f',breaksy2))
    }))
    rownames(dat2) <- rownames(dat3) <- rownames(dat4) <- NULL

    pp <- ggplot(dat, aes(x=x,y=y)) +
          geom_ribbon(data=dat3, aes(ymin=nsupmin,ymax=nsupmax,fill=nsup_trait),
                     color="white", size=0.2, alpha=0.25) +
          geom_text(data=dat4, aes(label=label), size=2.5, hjust=1.0) +
          geom_rect(data=dat4, aes(xmin=x1, xmax=Inf, ymin=y, ymax=y),
                    fill=NA, color="black") +
          #coord_cartesian(clip = 'off') +
          geom_line(size=lwd, color=line.color) +
          labs(title=main, x=xlab, y=ylab) +
          theme_bw() + theme0 +
          geom_vline(data=dat2, aes(xintercept=x),
                      size=0.5,linetype="dotted",color="gray50") +
          facet_wrap(~trait2, scales="free", labeller=label_parsed) +
          scale_y_continuous(limits=ylim, expand=expansion(mult=c(0.03)),
                             sec.axis=sec_axis(~., ylab2, breaks=NULL)) +
          guides(fill = guide_legend(override.aes = list(alpha=0.5)))

    if(x.stat=="nsup"){
       pp <- pp +
             scale_x_continuous(breaks=breaksx, labels=round(labelsx),
                                expand=expansion(mult=expand.x))
    }else{
       pp <- pp +
             scale_x_continuous(breaks=scales::extended_breaks(n=nbreaks.x),
                                limits=xlim,expand=expansion(mult=expand.x))
    }

    return(pp)
}
