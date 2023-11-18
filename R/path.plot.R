
path.plot <- function(object, Z = NULL, K = NULL, i = NULL,
                      prune = FALSE, cor.max = 0.97,
                      lambda.min = .Machine$double.eps^0.5,
                      nbreaks.x = 6, ...)
{
  # object=fm1; K=G; prune=TRUE; cor.max=0.9
  flagKinship <- FALSE
  eps <- .Machine$double.eps
  args0 <- list(...)

  if(!inherits(object, c("LASSO","SSI"))){
     stop("The input object is not of the class 'LASSO' or 'SSI'")
  }

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

  theme0 <- theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.background = element_rect(fill = "gray95"),
    legend.box.spacing = unit(0.4, "lines"),
    legend.key.height= unit(3,"line"),
    legend.key.width = unit(0.8, "lines")
  )

  if(inherits(object, "SSI"))
  {
    if(!is.null(K)){
      flagKinship <- TRUE
      if(!is.null(Z)){
        if(length(dim(Z)) != 2L) stop("Object 'Z' must be a matrix")
        K <- tcrossprod(Z,tcrossprod(Z,K))  # G = ZKZ'
      }
      if(length(dim(K)) != 2L | (length(K) != object$n^2)){
        stop("Product Z %*% K %*% t(Z) must be a squared matrix with dimension equal to the number of rows in 'Y'")
      }
    }

    beta <- coef.SSI(object)
    # Coarce nsup and lambda to lists
    object$lambda <- lapply(1:nrow(object$lambda),function(k)object$lambda[k,])
    object$nsup <- lapply(1:nrow(object$nsup),function(k)object$nsup[k,])

    if(length(object$tst) == 1L){
      beta <- list(beta)
    }

    if(!is.null(i)){
      if( max(i) > length(object$tst) ){
        stop("All elements in 'i' must be between 0 < i <= length(object$tst)")
      }
      indexTST <- i

    }else{ indexTST <- seq_along(object$tst) }

  }else{
    beta <- coef.LASSO(object)
    object$trn <- seq(object$p)
    object$tst <- seq(object$q)

    if(object$q == 1L){
      beta <- list(beta)
      object$lambda <- list(object$lambda)
      object$nsup <- list(object$nsup)
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
  object$nsup <- object$nsup[indexTST]

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
      message("'pruning' is applied for each response variable")
    }
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
      object$nsup <- object$nsup[seq(nTST0)]
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
    b0 <- do.call(rbind, beta[tst0])
    INDEX2 <- cbind(rep(tst0, each=length(object$trn)),
                    rep(seq_along(object$trn),length(tst0)))

    indexOK <- 1:nrow(b0)
    indexdrop <- which(apply(b0, 1L, function(x)all(abs(x) <= eps)))
    if(length(indexdrop) > 1){  # leave one that has all bij==0
      b0 <- b0[-indexdrop[-1], ,drop=FALSE]
      indexOK <- indexOK[-indexdrop[-1]]
    }

    if(prune & nrow(b0)>1){
      tmp <- Prune(t(b0), alpha=cor.max^2)$prune.in
      indexOK <- indexOK[tmp]
      b0 <- b0[tmp, ,drop=FALSE]
    }

    dat0 <- do.call(rbind,lapply(1:length(indexOK), function(j){
          tmp <- INDEX2[indexOK[j],]
          id <- paste0(object$tst[indexTST[tmp[1]]],":",object$trn[tmp[2]])

          if(flagKinship){
            Kij <- K[object$tst[indexTST[tmp[1]]],object$trn[tmp[2]]]
          }else{
            Kij <- NA
          }
          data.frame(nsup=object$nsup[[tmp[1]]],lambda=object$lambda[[tmp[1]]],
                     beta=b0[j,], Kij=Kij, id=id)
    }))

    dat <- rbind(dat, dat0)
  }
  dat$id <- factor(as.character(dat$id))

  # Labels and breaks for the nsup axis
  tmp <- get_breaks(x=-log(unlist(object$lambda)), y=unlist(object$nsup),
                    nbreaks.x=nbreaks.x, ymin=0)
  breaks0 <- tmp$breaks.x
  labels0 <- round(tmp$breaks.y)
  labels2 <- sprintf('%.2f', breaks0)

  if(flagKinship){
    pt <- ggplot(dat, aes(-log(lambda),beta,color=Kij,group=id)) +
          labs(color=expression(k[ij])) +
          viridis::scale_color_viridis() + geom_line(size=lwd) +
          theme_bw() + theme0

  }else{
    pt <- ggplot(dat, aes(-log(lambda),beta,color=id,group=id)) +
          geom_line(size=lwd) + theme_bw() + theme0 +
          theme(legend.position = "none")
  }

  pt <- pt + labs(title=main,y=ylab,x=xlab) +
        scale_x_continuous(breaks=breaks0,labels=labels0,
                    sec.axis=sec_axis(~.+0,expression(paste("-log(",lambda,")")),
                    breaks=breaks0, labels=labels2))
  pt

}
