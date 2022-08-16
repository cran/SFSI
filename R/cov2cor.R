
# Covariance matrix to correlation matrix

cov2cor2 <- function(V, a = 1, void = FALSE)
{
    if((sum(dim(V))/2)^2 != length(V)) stop("Object 'V' must be a squared matrix")
    if(!float::storage.mode(V) %in% c("double","float32")) storage.mode(V) <- "double"

    p <- ncol(V)
    isfloat <- as.logical(float::storage.mode(V)=="float32")

    #dyn.load("c_utils.so")
    if(void){
      if(isfloat){
       nOK <- .Call('cov2correlation', p, V@Data, isfloat, as.numeric(a))[[1]]
      }else{
       nOK <- .Call('cov2correlation', p, V, isfloat, as.numeric(a))[[1]]
      }
      out <- NULL
    }else{
      if(isfloat){
       out <- V@Data[]
      }else{
        out <- V[]
      }

     nOK <- .Call('cov2correlation', p, out, isfloat, as.numeric(a))[[1]]
     if(isfloat){
        out <- float::float32(out)
     }
   }
   #dyn.unload("c_utils.so")
   if(nOK != p){
      warning("Some diagonal values of 'V' are 0 or NA. Results are dobubtful",immediate.=TRUE)
   }
   out
}
