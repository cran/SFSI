
# Read a file as binary

readBinary <- function(file = paste0(tempdir(),"/file.bin"),
                  index.row = NULL, index.col = NULL, verbose = TRUE)
{
  if(!file.exists(file)){
    stop("File '",file,"' does not exist")
  }

  nsetRow <- as.integer(length(index.row))
  nsetCol <- as.integer(length(index.col))

  # Read lines
  X <- .Call("readBinFileFloat",file,nsetRow,nsetCol,
             as.integer(index.row),as.integer(index.col))

  n <- X[[1]]; p <- X[[2]]; size <- X[[3]]
  isFloat <- X[[4]]
  nError <- X[[5]]

  if(nError == 0L){
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
