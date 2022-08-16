
# Save a file as binary

saveBinary <- function(X, file = paste0(tempdir(),"/file.bin"),
              precision.format = c("double","single"), verbose = TRUE)
{
  precision.format <- match.arg(precision.format)
  type <- float::storage.mode(X)

  if(length(dim(X)) != 2L) stop("Object 'X' must be a matrix")
  if(!type %in% c("logical","double","integer","float32")){
     stop("'storage.mode(X)' must be either 'logical', 'float32', 'double' , 'integer'")
  }

  unlink(file)

  if(type=="float32" & precision.format!='single'){
    precision.format <- 'single'
    cat("Variable type is set to precision.format='single' when storage.mode(X)='float'\n")
  }
  isfloat <- as.logical(type=="float32")
  precision <- ifelse(precision.format=="single", 1L, 2L)

  #dyn.load("c_utils.so")
  if(isfloat){
    out <- .Call('writeBinFile', file, nrow(X), ncol(X), X@Data, isfloat, precision)
   }else{
    out <- .Call('writeBinFile', file, nrow(X), ncol(X), X, isfloat, precision)
  }
  #dyn.unload("c_utils.so")

  if(verbose){
    size <- out[[5]]
    tmp <- c(Gb=1E9,Mb=1E6,Kb=1E3,b=1E0)
    sz <- file.size(file)/tmp[min(which(file.size(file)/tmp>1))]
    cat("Saved file '",file,"'\n")
    cat("     nrow=",nrow(X),", ncol=",ncol(X),", type=",type,", size=",size,"bytes, file.size=",round(sz,2),names(sz),"\n")
  }
}
