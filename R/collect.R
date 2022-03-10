
# Collect all outputs when divided acording to 'subset' parameter

collect <- function(prefix = "")
{
  filenames <- Sys.glob(paste0(prefix,"_*_of_*.RData"))
  out <- NULL
  if(length(filenames) > 0){
      nFiles <- as.numeric(unlist(lapply(strsplit(filenames,"_"),function(x) gsub(".RData","",x[length(x)]))))
      if(length(unique(nFiles))>1)
        stop(" Different subset output files were found for the given prefix='",prefix,
            "'. Remove old files. No output was collected")

      filenames <- paste0(prefix,"_",1:nFiles[1],"_of_",nFiles[1],".RData")
      if(!all(file.exists(filenames))) stop("Some files are missing for the given prefix='",prefix,"'\n")

      for(i in seq_along(filenames))
      {
        load(filenames[i])
        if(i==1){
          fm <- out
        }else{
          fm$file_beta <- c(fm$file_beta,out$file_beta)
          fm$tst <- c(fm$tst,out$tst)
          fm$df <- rbind(fm$df,out$df)
          fm$lambda <- rbind(fm$lambda,out$lambda)
        }
        cat(" Loaded file: '",filenames[i],"'\n",sep="")
      }
      fm$subset[1] <- NA

  }else stop(" No output files were found for the given prefix='",prefix,"'")

  fm
}
