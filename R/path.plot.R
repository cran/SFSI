
# Plot the coefficients path in a penalized regression (user-level)

# Z=NULL; K=NULL; tst=NULL; main=NULL; cor.max=0.85

path.plot <- function(object, Z = NULL, K = NULL, tst = NULL,
  cor.max = 0.85, nbreaks.x = 6, ...)
{
  k <- NULL
  flagKinship <- FALSE
  args0 <- list(...)

  if(!inherits(object, c("LASSO","SSI"))) stop("The input object is not of the class 'LASSO' or 'SSI'")

  if(inherits(object, "SSI") & "CV" %in% names(object))
    stop("'path.plot' cannot be applied after cross-validation")

  xlab <- "Number of active predictors"
  ylab <- "beta"
  main <- NULL
  lwd <- ifelse("lwd" %in% names(args0), args0$lwd, 0.4)
  if("main" %in% names(args0)) main <- args0$main
  if("xlab" %in% names(args0)) xlab <- args0$xlab
  if("ylab" %in% names(args0)) ylab <- args0$ylab

  theme0 <- ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.background = ggplot2::element_rect(fill = "gray95"),
    legend.box.spacing = ggplot2::unit(0.4, "lines"),
    legend.key.height= ggplot2::unit(3,"line"),
    legend.key.width = ggplot2::unit(0.8, "lines")
  )

  if(inherits(object, "SSI"))
  {
    if(!is.null(K)){
      flagKinship <- TRUE
      if(is.character(K)){
        K <- readBinary(K)
      }
      if(!is.null(Z))
      {
        if(length(dim(Z)) != 2L) stop("Object 'Z' must be a matrix")
        K <- float::tcrossprod(Z,float::tcrossprod(Z,K))  # G = ZKZ'
      }
      if(length(dim(K)) != 2L | (length(K) != length(object$y)^2))
        stop("Product Z %*% K %*% t(Z) must be a squared matrix with number of rows (and columns) equal to the number of elements in 'y'")
    }

    beta <- coef.SSI(object)
    if(!is.null(tst)){
      if(any(!tst %in% object$tst)) stop("Some elements in 'tst' are not contained in 'object$tst'")
      indexTST <- which(object$tst %in% tst)
    }else indexTST <- seq_along(object$tst)
    beta <- beta[indexTST]
    lambda <- apply(object$lambda, 2, mean)
    df <- apply(object$df, 2, mean)

  }else{
    beta <- object$beta
    lambda <- object$lambda
    df <- object$df
  }

  ndf <- length(df)
  if(min(lambda) < .Machine$double.eps*1000)  lambda[which.min(lambda)] <- min(lambda[lambda>0])/2
  if(ndf == 1) stop("Coefficients path plot can not be generated for 'nlambda=1'")

  if(inherits(object, "SSI"))
  {
    dat <- c()
    trim <- length(object$trn)*length(beta) > 20000
    for(i in seq_along(beta))
    {
      b0 <- as.matrix(beta[[i]])
      if(trim){
        indexOK <- getIndexCorrelated(t(b0), cor.max)
      }else indexOK <- seq_along(object$trn)

      tmp <- matrix(NA,nrow=1,ncol=length(indexOK))
      if(!is.null(K)){
        tmp <- K[object$tst[indexTST[i]],object$trn[indexOK],drop=FALSE]
        if(float::storage.mode(tmp)=='float32') tmp <- float::dbl(tmp)
      }
      dimnames(tmp) <- list(object$tst[indexTST[i]], object$trn[indexOK])
      tmp <- reshape2::melt(tmp)
      colnames(tmp) <- c("tst_i","trn_i","value")
      tmp <- tmp[rep(1:nrow(tmp), ndf),]

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
  tmp <- get_breaks(lambda, df, nbreaks=nbreaks.x, ymin=0)
  breaks0 <- tmp$breaks.x
  labels0 <- round(tmp$breaks.y)
  labels2 <- sprintf('%.2f', breaks0)

  if(flagKinship)
  {
    pt <- ggplot2::ggplot(dat,ggplot2::aes(-log(lambda),beta,color=k,group=id)) +
      viridis::scale_color_viridis() + ggplot2::geom_line(size=lwd) + ggplot2::theme_bw() + theme0
  }else{
    pt <- ggplot2::ggplot(dat,ggplot2::aes(-log(lambda),beta,color=id,group=id))+
      ggplot2::geom_line(size=lwd) + ggplot2::theme_bw() + theme0 +
      ggplot2::theme(legend.position = "none")
  }

  pt <- pt + ggplot2::labs(title=main,y=ylab,x=xlab) +
             ggplot2::scale_x_continuous(breaks=breaks0,labels=labels0,
                      sec.axis=ggplot2::sec_axis(~.+0,expression(paste("-log(",lambda,")")),
                      breaks=breaks0, labels=labels2))
  pt
}
