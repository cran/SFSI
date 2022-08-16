
# Collect all outputs when divided acording to 'subset' parameter

collect <- function(prefix = "")
{
  filenames <- Sys.glob(paste0(prefix,"subset_*_of_*_output.RData"))
  out <- NULL
  if(length(filenames) > 0){
      tmp <- strsplit(basename(filenames),"_")
      nFiles <- as.numeric(unlist(lapply(tmp,function(x) x[length(x)-1])))
      if(length(unique(nFiles))>1){
        stop(" Different subset output files were found for the given prefix='",prefix,
            "'. Remove old files. No output was collected")
      }

      filenames <- paste0(prefix,"subset_",1:nFiles[1],"_of_",nFiles[1],"_output.RData")
      if(!all(file.exists(filenames))) stop("Some files are missing for the given prefix='",prefix,"'\n")

      for(i in seq_along(filenames))
      {
        load(filenames[i])
        if(i==1){
          fm <- out
          fm$name_beta <- vector('list',nFiles[1])
          fm$name_beta[[1]] <- out$name_beta

        }else{
          fm$file_beta <- c(fm$file_beta, out$file_beta)
          fm$name_beta[[i]] <- out$name_beta
          fm$tst <- c(fm$tst, out$tst)
          fm$u <- rbind(fm$u, out$u)
          fm$df <- rbind(fm$df, out$df)
          fm$lambda <- rbind(fm$lambda, out$lambda)
        }
        cat(" Loaded file: '",filenames[i],"'\n",sep="")
      }
      fm$subset$subset[1] <- NA

  }else stop(" No output files were found for the given prefix='",prefix,"'")

  fm
}
