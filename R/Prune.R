#====================================================================
# Prunne elements based on pairwise similarity R (correlation, R^2, etc)
#====================================================================

Prune <- function(R, threshold=0.95, D=NULL, d.max=NULL, verbose=FALSE)
{
    nc <- nchar(ncol(R))
    PrFn <- function(x, nc=3)sprintf(paste0('%',nc,'d'),x)

    A <- (R > threshold)
    diag(A) <- FALSE
    pruneIn <- c()
    pruneOut <- c()

    if(!is.null(d.max) & is.null(D)){
      cat("Maximum distance 'd.max' is ignored when matrix 'D=NULL'\n")
    }
    if(is.null(d.max) & !is.null(D)){
      cat("Distance matrix 'D' is ignored when 'd.max=NULL'\n")
      D <- NULL
    }
    if(!is.null(d.max) & !is.null(D)){
      colnames(D) <- rownames(D) <- NULL
      D2 <- (D <= d.max)
      diag(D2) <- FALSE
      A <- A & D2
    }

    if(verbose){
      tmp <- ifelse(is.null(d.max),"",paste0(" within ",d.max," Mb"))
      cat("Pruning ", ncol(R)," subjects",tmp," ...\n",sep="")
    }

    ID <- 1:ncol(R)

    CON <- as.vector(rowSums(A))
    nConn <- CON[]
    remain <- 1:ncol(R)
    mc <- nchar(max(CON))

    cont <- 0
    if(any(CON==0)){
        tmp <- which(CON==0)
        pruneIn <- ID[tmp]
        remain <- remain[-tmp]
        CON <- CON[-tmp]
        ID <- ID[-tmp]
        if(verbose){
          cont <- cont + 1
          cat("--------------------------------------------------------\n")
          cat(" S",PrFn(cont),". nConn=",PrFn(0,mc), ". In: n=",PrFn(length(tmp),nc),
              ". Out: n=",PrFn(0,mc),". Remain: n=",PrFn(length(remain),nc),"\n",sep="")
        }
    }

    CON0 <- rep(0,length(CON))
    flag <- any(CON>0)
    while(flag){
      if(any(CON > 0)){
        tmp <- which.max(CON)
        pruneIn <- c(pruneIn, ID[tmp])
        remove <- which(A[remain[tmp],remain]) # conections with the maximum

        pruneOut <- c(pruneOut, ID[remove])

        if(verbose){
          cont <- cont + 1
          cat("--------------------------------------------------------\n")
          cat(" S",PrFn(cont),". nConn=",PrFn(CON[tmp],mc),". In: i=",PrFn(tmp,nc),
                ". Out: n=",PrFn(length(remove),mc),". Remain: n=",
                PrFn(ncol(R)-length(pruneOut)-length(pruneIn),nc),"\n",sep="")
        }

        CON0 <- rowSums(A[remain[-remove], remain[remove], drop=FALSE])
        remain <- remain[-remove]
        CON <- CON[-remove]
        ID <- ID[-remove]
        CON <- CON-CON0
        flag=(length(remain) > 1)
      }else{
         flag <- FALSE
         tmp <- which(!ID %in% pruneIn)
         if(length(tmp)>0){
          pruneIn <- c(pruneIn, ID[tmp])
          if(verbose){
            cont <- cont + 1
            cat("--------------------------------------------------------\n")
            cat(" S",PrFn(cont),". nConn=",PrFn(max(CON),mc),
                ". In: n=",PrFn(length(tmp),nc),". Out: n=",PrFn(0,mc),". Remain: n=",
                PrFn(ncol(R)-length(pruneOut)-length(pruneIn),nc),"\n",sep="")
          }
        }
      }
    }
    return(list(nConn=nConn, prune.in=sort(pruneIn), prune.out=sort(pruneOut)))
}
