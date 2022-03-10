
# Save a file as binary

saveBinary <- function(X, file = paste0(tempdir(),"/file.bin"),
              type = c("float","double"), verbose = TRUE)
{
  type <- match.arg(type)

  if(length(dim(X)) != 2L) stop("Object 'X' must be a matrix")
  if(!float::storage.mode(X) %in% c("double","float32")) storage.mode(X) <- "double"

  unlink(file)

  if(float::storage.mode(X)=="float32" & type!='float'){
    type <- 'float'
    warning("Object can be only saved as type='float' when class(X)='float'\n",
            "  Variable type was changed to type='float'",immediate.=TRUE)
  }
  isFloat <- float::storage.mode(X)=="float32"
  size <- ifelse(type=="float",4,8)

  if(isFloat){
    out <- .Call('writeBinFileFloat',file,nrow(X),ncol(X),
           as.integer(size),X@Data,isFloat)
   }else{
    out <- .Call('writeBinFileFloat',file,nrow(X),ncol(X),
           as.integer(size),X,isFloat)
  }

  if(verbose){
    tmp <- c(Gb=1E9,Mb=1E6,Kb=1E3,b=1E0)
    sz <- file.size(file)/tmp[min(which(file.size(file)/tmp>1))]
    cat("Saved file '",file,"'\n")
    cat("     nrow=",nrow(X),", ncol=",ncol(X),", type=",type,", size=",size,"bytes, file.size=",round(sz,2),names(sz),"\n")
  }
}
