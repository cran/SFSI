
path.plot <- function(object, Z = NULL, K = NULL, i = NULL,
                      prune = FALSE, cor.max = 0.97,
                      lambda.min = .Machine$double.eps^0.5,
                      nbreaks.x = 6, ...)
{
  flagKinship <- FALSE
  eps <- .Machine$double.eps
  args0 <- list(...)

  if(!inherits(object, c("LASSO","SSI"))) stop("The input object is not of the class 'LASSO' or 'SSI'")

  if(inherits(object, "SSI") & "CV" %in% names(object)){
    stop("'path.plot' cannot be applied after cross-validation")
  }

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
    object$lambda <- lapply(1:nrow(object$lambda),function(k)object$lambda[k,])
    object$df <- lapply(1:nrow(object$df),function(k)object$df[k,])

    if(length(object$tst) == 1L){
      beta <- list(beta)
    }

    if(!is.null(i)){
      if( max(i) > length(object$tst) ){
        stop("All elements in 'i' must be between 0 < i <= length(object$tst)")
      }
      indexTST <- i

    }else indexTST <- seq_along(object$tst)

  }else{
    beta <- coef.LASSO(object)
    object$trn <- seq(object$p)
    object$tst <- seq(object$q)

    if(object$q == 1L){
      beta <- list(beta)
      object$lambda <- list(object$lambda)
      object$df <- list(object$df)
    }
    if(!is.null(i)){
      if( max(i) > object$q ){
        stop("All elements in 'i' must be between 0 < i <= q, where q=ncol(Gamma)")
      }
      indexTST <- i
    }else{
      indexTST <- seq_along(object$tst)
    }
  }

  beta <- beta[indexTST]
  object$lambda <- object$lambda[indexTST]
  object$df <- object$df[indexTST]

  nLambda <- unlist(lapply(object$lambda,length))
  if(all(nLambda < 5L)){
    stop("Coefficients path plot can be generated for at least 5 lambda values")
  }

  if(any(unlist(lapply(object$lambda,min)) < lambda.min)){
    min0 <- min(unlist(lapply(object$lambda, function(x)min(x[x>lambda.min+eps]))))
    for(k in 1:length(object$lambda)){
      tmp <- object$lambda[[k]]
      object$lambda[[k]] <- ifelse(tmp < lambda.min, min0/10, tmp)
    }
  }

  if(length(unique(nLambda)) > 1L){
    INDEX1 <- matrix(seq_along(indexTST),ncol=1)
    if(prune){
      message("'pruning' is applied for each response variable") }
  }else{
    if(length(object$trn)*length(beta) > 10000){
      nTST0 <- ceiling(10000/length(object$trn))
      message("The number of paths is very large. Only ",nTST0*length(object$trn),
              " paths corresponding to the first ",nTST0,
              ifelse(inherits(object, "SSI")," testing elements"," response variables"),
              " are considered.")
      message("You can select specific paths through 'i' argument")

      beta <- beta[seq(nTST0)]
      object$lambda <- object$lambda[seq(nTST0)]
      object$df <- object$df[seq(nTST0)]
      indexTST <- indexTST[seq(nTST0)]
    }
    nTSTi <- ceiling(5000/length(object$trn))
    INDEX1 <- matrix(seq(nTSTi*ceiling(length(indexTST)/nTSTi)),ncol=nTSTi, byrow=TRUE)
    if(prune & nrow(INDEX1)>1L){
      message("'pruning' is applied in groups of ",nTSTi*length(object$trn)," paths")
    }
  }

  id <- Kij <- lambda <- NULL
  dat <- c()
  for(k in 1:nrow(INDEX1))
  {
    tst0 <- INDEX1[k,][INDEX1[k,] <= length(indexTST)]
    b0 <- as.matrix(do.call(rbind,beta[tst0]))
    INDEX2 <- cbind(rep(tst0, each=length(object$trn)),
                    rep(seq_along(object$trn),length(tst0)))

    indexOK <- 1:nrow(b0)
    indexdrop <- which(apply(b0, 1, function(x)all(abs(x) <= eps)))
    if(length(indexdrop) > 1){  # leave one that has all bij==0
      b0 <- b0[-indexdrop[-1], ,drop=FALSE]
      indexOK <- indexOK[-indexdrop[-1]]
    }

    if(prune & nrow(b0)>1){
      tmp <- Prune(R=cor(t(b0))^2, threshold=cor.max^2)$prune.in
      indexOK <- indexOK[tmp]
      b0 <- b0[tmp, ,drop=FALSE]
    }

    dat0 <- do.call(rbind,lapply(1:length(indexOK), function(j){
          tmp <- INDEX2[indexOK[j],]
          id <- paste0(object$tst[indexTST[tmp[1]]],":",object$trn[tmp[2]])

          if(flagKinship){
            Kij <- float::dbl(K[object$tst[indexTST[tmp[1]]],object$trn[tmp[2]]])
          }else{
            Kij <- NA
          }
          data.frame(df=object$df[[tmp[1]]],lambda=object$lambda[[tmp[1]]],
                     beta=b0[j,], Kij=Kij, id=id)
    }))

    dat <- rbind(dat, dat0)
  }
  dat$id <- factor(as.character(dat$id))

  # Labels and breaks for the DF axis
  tmp <- get_breaks(unlist(object$lambda), unlist(object$df),
                    nbreaks=nbreaks.x, ymin=0)
  breaks0 <- tmp$breaks.x
  labels0 <- round(tmp$breaks.y)
  labels2 <- sprintf('%.2f', breaks0)

  if(flagKinship){
    pt <- ggplot2::ggplot(dat,ggplot2::aes(-log(lambda),beta,color=Kij,group=id)) +
          ggplot2::labs(color=expression(k[ij])) +
          viridis::scale_color_viridis() + ggplot2::geom_line(size=lwd) +
          ggplot2::theme_bw() + theme0

  }else{
    pt <- ggplot2::ggplot(dat,ggplot2::aes(-log(lambda),beta,color=id,group=id)) +
      ggplot2::geom_line(size=lwd) + ggplot2::theme_bw() + theme0 +
      ggplot2::theme(legend.position = "none")
  }

  pt <- pt + ggplot2::labs(title=main,y=ylab,x=xlab) +
             ggplot2::scale_x_continuous(breaks=breaks0,labels=labels0,
                      sec.axis=ggplot2::sec_axis(~.+0,expression(paste("-log(",lambda,")")),
                      breaks=breaks0, labels=labels2))
  pt

}
