
# Read a file as binary

readBinary <- function(file = paste0(tempdir(),"/file.bin"),
                  index.row = NULL, index.col = NULL, verbose = TRUE)
{
  if(!file.exists(file)){
    stop("File '",file,"' does not exist")
  }

  nsetRow <- as.integer(length(index.row))
  nsetCol <- as.integer(length(index.col))

  maxsetRow <- as.integer(ifelse(nsetRow>0, max(index.row), 0))
  maxsetCol <- as.integer(ifelse(nsetCol>0, max(index.col), 0))

  # Read lines
  #dyn.load("c_utils.so")
  X <- .Call("readBinFile", file, nsetRow, nsetCol,
             maxsetRow, maxsetCol,
             as.integer(index.row),as.integer(index.col))
  #dyn.unload("c_utils.so")

  n <- X[[1]]
  p <- X[[2]]
  isfloat <- as.logical(X[[3]])
  vartype <- X[[4]]
  size <- X[[5]]
  nError <- X[[6]]

  type <- ifelse(vartype==1L,"integer",ifelse(vartype==2L,"logical","double"))
  if(nError == 0L){
    if(isfloat | (type=="double" & size==4)){
      X <- float::float32(X[[7]])
    }else{
      X <- X[[7]]
    }
    if(type=="logical"){
      storage.mode(X) <- "logical"
    }
    type <- ifelse(isfloat, "float32", type)
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
