#====================================================================
# Plot the top 2 PCs of the K matrix showing tst and trn points
#====================================================================
# Z = NULL; K=R0; group = group.shape = set.color = set.size = tst= df=NULL; eps=.Machine$double.eps
# axis.labels = TRUE; curve = FALSE; bg.color = "white"; unified = TRUE; ntst = 36;
# line.color = "gray90"; line.tick = 0.3; legend.pos="right"; show.names = TRUE
# point.color = "gray20"; sets = c("Testing","Supporting","Non-active")

net.plot <- function(object, Z = NULL, K = NULL, tst = NULL,
           show.names = FALSE, group = NULL, group.shape = NULL,
           set.color = NULL, set.size = NULL, df = NULL, main,
           axis.labels = TRUE, curve = FALSE, bg.color = "white",
           unified = TRUE, ntst = 36, line.color = "gray80",
           line.tick = 0.3, legend.pos = "right", point.color = "gray20",
           sets = c("Testing","Supporting","Non-active"),
           eps = .Machine$double.eps)
{
  set <- set_name <- label <- x <- y <- x_TRN <- x_TST <- y_TRN <- y_TST <- NULL
  xxx <- yyy <- isEigen <- NULL
  legend.pos <- match.arg(legend.pos,
    choices=c("right","bottomright","bottomleft","topleft","topright","none"))

  if(!is.null(K)){
    if(is.character(K))  K <- readBinary(K)
    if(length(dim(K)) != 2L | (length(K) != nrow(K)^2)) {
        stop("Input 'K' must be a squared matrix")
    }
    if(!is.null(Z)) {
      if(length(dim(Z)) != 2L) stop("Object 'Z' must be a matrix with ncol(Z)=nrow(K)\n")
      K <- float::tcrossprod(Z,float::tcrossprod(Z,K))   # Z%*%K%*%t(Z)
    }
  }

  isSSI <- FALSE
  if(inherits(object, "SSI")){
    X <- NULL
    isSSI <- TRUE
  }else{
    if(length(dim(object)) == 2L){
      X <- object
    }else stop("The input object is not of the class 'SSI' or a matrix")
  }

  if(isSSI){
    if(is.null(df)) df <- summary.SSI(object)$optCOR['df']
    if(0 > df | df > range(object$df)[2])
      stop("Parameter 'df' must be greater than zero and no greater than 'trn' size")
    X <- as.matrix(coef.SSI(object, df=df))

    if(is.null(tst)){
      yyy <- object$tst
    }else{
      yyy <- tst
      if(any(!tst %in% object$tst))
        stop("Some elements in 'tst' are not contained in 'object$tst'")
    }

    X <- X[object$tst %in% yyy, ,drop=FALSE]
    xxx <- object$trn

    if(!is.null(object$id)){
      rownames(X) <- object$id[yyy]
      colnames(X) <- object$id[xxx]
    }else{
      rownames(X) <- yyy
      colnames(X) <- xxx
    }

    if(!is.null(K)){
      if(length(object$y) != nrow(K)){
        if(!is.null(object$id) & all(unlist(dimnames(X)) %in% rownames(K))){
          xxx <- yyy <- NULL
        }else{
          cat("Input 'object' couldn't be linked to 'K'. Input 'K' will be ignored\n")
          K <- xxx <- yyy <- NULL
        }
      }else{
        if(is.null(object$id)){
          dimnames(K) <- list(1:nrow(K),1:ncol(K))
        }else dimnames(K) <- list(object$id,object$id)
        if(any(!yyy %in% 1:nrow(K))) stop("Some 'object$tst' indices are larger than 'ncol(K)'")
        if(any(!xxx %in% 1:nrow(K))) stop("Some 'object$trn' indices are larger than 'ncol(K)'")
      }
    }else{
      xxx <- yyy <- NULL
    }
  }

  net <- get_net(X, K=K, xxx=xxx, yyy=yyy, eps=eps)
  isSymmetric <- net$isSymmetric
  labelsAxis <- net$labelsAxis
  xxx <- net$xxx
  yyy <- net$yyy
  isEigen <- net$isEigen

  if(!isSSI & isSymmetric){
      legend.pos <- "none"
      if(!unified){
        cat("Only an 'unified' plot can be produced with the input object data\n")
        unified <- TRUE
      }
  }

  if(!unified & length(yyy) >= ntst){
    cat("Large number of testing individuals. Only the first",ntst,"are shown\n")
    yyy <- yyy[1:ntst]
  }

  justx <- ifelse(length(grep("left",legend.pos))>0,0,1)
  justy <- ifelse(length(grep("bottom",legend.pos))>0,0,1)
  if(!legend.pos %in% c("none","right")) legend.pos <- c(abs(justx-0.01),abs(justy-0.01))

  flagGp <- !is.null(group)
  if(is.null(group)) group <- data.frame(group=rep(1,nrow(net$xy)))
  gpName <- colnames(group)

  if(!(class(sets) == "character" & length(sets) == 3))
   stop("Parameter 'sets' must be a triplet of 'character' type")

  dat <- data.frame(id=1:nrow(net$xy),label=net$labels,set=net$set,
                    set_name=sets[net$set],group=group,float::dbl(net$xy))
  if(any(net$set==4)) dat$set_name[net$set==4] <- sets[1]

  dat$group <- factor(as.character(dat$group))
  dat$set_name <- factor(as.factor(dat$set_name),levels=c(sets))

  if(length(show.names)==1L) show.names <- rep(show.names, 3)

  # Shape and color for the levels of group
  if(!flagGp) dat$group <- dat$set_name
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
  }else if(length(set.color)==1) set.color <- rep(set.color,length(sets))
  set.color <- set.color[1:length(sets)]

  if(is.null(set.size)){
    set.size <- c(3.1, 2.1, 0.8)
    if(any(show.names)) set.size[show.names] <- 3.1
  }else if(length(set.size)==1L) set.size <- rep(set.size,length(sets))
  set.size <- set.size[1:length(sets)]

  if(any(is.na(group.shape)))
    stop("The number of elements in 'group.shape' must be of length ",length(levelsGp))

  if(any(is.na(set.size)) | any(is.na(set.color)))
    stop("The number of elements in 'set.size' and 'set.color' must be of length ",length(sets))

  theme0 <- ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
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

  main0 <- NULL
  if(missing(main)){
     if(isSSI){
       main0 <- bquote(.(object$name)*". Support set size="*.(round(df)))
     }else theme0 <- theme0 + ggplot2::theme(plot.title = ggplot2::element_blank())

  }else main0 <- main

  if(is.null(main0)){
    theme0 <- theme0 + ggplot2::theme(plot.title = ggplot2::element_blank())
  }

  if(!axis.labels){
    theme0 <- theme0 + ggplot2::theme(axis.text=ggplot2::element_blank(),
                        axis.ticks=ggplot2::element_blank())
  }

  if(!isEigen){
    theme0 <- theme0 + ggplot2::theme(axis.title=ggplot2::element_blank(),
                               axis.text=ggplot2::element_blank(),
                               axis.ticks=ggplot2::element_blank())
  }

  names(group.shape) <- levelsGp
  names(set.color) <- names(set.size) <- sets

  # If unified plot
  if(unified){
    pt <- ggplot2::ggplot(dat,ggplot2::aes(x=x,y=y))
    if(show.names[3]){
      pt <- pt + ggplot2::geom_label(data=dat[dat$set==3,],ggplot2::aes(label=label,fill=set_name),
               label.padding=ggplot2::unit(0.15,"lines"),color=point.color,size=set.size[3])

    }else{
           pt <- pt + ggplot2::geom_point(data=dat[dat$set==3,],ggplot2::aes(shape=group,fill=set_name),
                    color=point.color,size=set.size[3])
    }

    for(i in 1:length(yyy))
    {
      xxx0 <- net$edges[[i]]
      if(length(xxx0)>0)
      {
        dat1 <- dat[xxx0, c("x","y")]
        dat2 <- dat[yyy[i], c("x","y")]
        colnames(dat1) <- paste0(colnames(dat1),"_TRN")
        colnames(dat2) <- paste0(colnames(dat2),"_TST")
        dat1 <- data.frame(dat2[rep(1,nrow(dat1)),],dat1)
        if(curve){
          pt <- pt + ggplot2::geom_curve(ggplot2::aes(x=x_TST,y=y_TST,xend=x_TRN,yend=y_TRN),
                        data=dat1,alpha=0.4,size=line.tick,color=line.color,curvature=0.4)
        }else{
          pt <- pt + ggplot2::geom_segment(ggplot2::aes(x=x_TST,y=y_TST,xend=x_TRN,yend=y_TRN),
                        data=dat1,alpha=0.4,size=line.tick,color=line.color)
        }
      }
    }
    # Nodes in rows
    if(show.names[1]){
      pt <- pt + ggplot2::geom_label(data=dat[dat$set%in%c(1,4),],ggplot2::aes(label=label,fill=set_name),
                          label.padding=ggplot2::unit(0.15,"lines"),color=point.color,size=set.size[1])
    }else{
      pt <- pt  +
        ggplot2::geom_point(data=dat[dat$set%in%c(1,4),],ggplot2::aes(shape=group,fill=set_name),
                            color=point.color,size=set.size[1])
    }
    # Nodes in columns
    if(show.names[2]){
      pt <- pt + ggplot2::geom_label(data=dat[dat$set==2,],ggplot2::aes(label=label,fill=set_name),
                          label.padding=ggplot2::unit(0.15,"lines"),color=point.color,size=set.size[2])
    }else{
      pt <- pt  + ggplot2::geom_point(data=dat[dat$set==2,],ggplot2::aes(shape=group,fill=set_name),
                            color=point.color,size=set.size[2])
    }
    # Nodes that are in both rows and columns (based on row/column names)
    if(any(dat$set==4)){
      if(show.names[1] | show.names[2]){
        pt <- pt +
              ggplot2::geom_label(data=dat[dat$set==4,],label=" ",fill=set.color[sets[2]],
                  label.padding=ggplot2::unit(0.135,"lines"),label.r=ggplot2::unit(0.35,"lines"),
                  color=set.color[sets[1]],size=set.size[2]) +
              ggplot2::geom_text(data=dat[dat$set==4,],ggplot2::aes(label=label),
                  color=point.color,size=set.size[2])
      }else{
        pt <- pt  +
           ggplot2::geom_point(data=dat[dat$set==4,],ggplot2::aes(shape=group),fill=set.color[sets[2]],
                             color=set.color[sets[1]],size=set.size[1]*0.55)
      }

    }
    pt <- pt + ggplot2::theme_bw() + theme0

  }else{
      set.size <- 0.7*set.size
      dat2 <- c()
      for(i in 1:length(yyy)){
        xxx0 <- net$edges[[i]]
        if(length(xxx0) > 0){
          tmp <- dat[-xxx0,]
          tmp$set <- 3; tmp$set_name <- sets[3]
          tmp2 <- dat[xxx0, ]
          tmp2$set <- 2; tmp2$set_name <- sets[2]
          tmp <- rbind(tmp, tmp2, dat[yyy[i], ])
          dat2 <- rbind(dat2,data.frame(tmp, ind = i))
        }
      }

      pt <- ggplot2::ggplot(dat2,ggplot2::aes(x=x,y=y)) + ggplot2::facet_wrap(~ind) +
            ggplot2::geom_point(data=dat2[dat2$set_name==sets[3],],ggplot2::aes(fill=set_name,shape=group),color=point.color,size=set.size[3]) +
            ggplot2::geom_point(data=dat2[dat2$set_name==sets[2],],ggplot2::aes(fill=set_name,shape=group),color=point.color,size=set.size[2]) +
            ggplot2::geom_point(data=dat2[dat2$set_name==sets[1],],ggplot2::aes(fill=set_name,shape=group),color=point.color,size=set.size[1]) +
            ggplot2::theme_bw() + theme0

  }

  pt <- pt + ggplot2::labs(title=main0, x=labelsAxis[1],y=labelsAxis[2]) +
    ggplot2::scale_shape_manual(values = group.shape,
              guide=ggplot2::guide_legend(override.aes=list(size=2,fill="white"))) +
    ggplot2::scale_fill_manual(values = set.color,
              guide=ggplot2::guide_legend(override.aes=list(shape=21,size=2)))

  if(!flagGp) pt <- pt + ggplot2::guides(shape="none")

  pt
}
