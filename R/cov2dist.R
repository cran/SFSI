
# Covariance matrix to distance matrix

cov2dist <- function(V, void = FALSE)
{
    if((sum(dim(V))/2)^2 != length(V)) stop("Object 'V' must be a squared matrix")
    if(!float::storage.mode(V) %in% c("double","float32")) storage.mode(V) <- "double"

    p <- ncol(V)
    isFloat <- float::storage.mode(V)=="float32"

    #dyn.load("c_utils.so")
    if(void)
    {
     if(isFloat){
      out <- .Call('cov2distance',as.integer(p),V@Data,isFloat)
     }else{
      out <- .Call('cov2distance',as.integer(p),V,isFloat)
     }
    }else{
      if(isFloat){
       out <- V@Data[]
     }else out <- V[]

     tmp <- .Call('cov2distance',as.integer(p),out,isFloat)
     if(isFloat) out <- float::float32(out)
   }
   #dyn.unload("c_utils.so")
   out
}
