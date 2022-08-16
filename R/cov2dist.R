
# Covariance matrix to distance matrix

cov2dist <- function(V, void = FALSE)
{
    if((sum(dim(V))/2)^2 != length(V)) stop("Object 'V' must be a squared matrix")
    if(!float::storage.mode(V) %in% c("double","float32")) storage.mode(V) <- "double"

    p <- ncol(V)
    isfloat <- as.logical(float::storage.mode(V)=="float32")

    #dyn.load("c_utils.so")
    if(void){
      if(isfloat){
       out <- .Call('cov2distance', p, V@Data, isfloat)
      }else{
       out <- .Call('cov2distance', p, V, isfloat)
      }
    }else{
      if(isfloat){
       out <- V@Data[]
     }else{
       out <- V[]
     }

     tmp <- .Call('cov2distance', p, out, isfloat)
     if(isfloat){
        out <- float::float32(out)
     }
   }
   #dyn.unload("c_utils.so")
   out
}
