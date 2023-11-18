
# Collect all outputs when divided acording to 'subset' parameter

collect <- function(path = ".")
{

  pattern <- "subset_[0-9]+_of_[0-9]+_output.RData$"
  isFolder <- FALSE
  if(dir.exists(path)){   # If is a folder
    infolder <- normalizePath(path)
    prefix <- NULL
    isFolder <- TRUE
  }else{
    infolder <- dirname(path)
    prefix <- basename(path)
    pattern <- paste0(prefix,pattern)
  }
  pattern <- paste0("^",pattern)
  listfiles <- list.files(infolder)
  filenames <- grep(pattern, listfiles, value=TRUE)
  #filenames <- normalizePath(paste0(infolder,"/",filenames))

  out <- NULL
  if(length(filenames) > 0){
      tmp <- strsplit(basename(filenames),"subset_|_of_|_output.RData")
      indexFile <- as.numeric(unlist(lapply(tmp,function(x) x[2])))
      nFiles <- as.numeric(unlist(lapply(tmp,function(x) x[length(x)])))

      if(length(unique(nFiles)) == 1L){
        if(length(nFiles) != nFiles[1]){
          stop("Found ",length(nFiles)," output file(s), ",nFiles[1]," are needed ",
               "\n  Remove old files. No output was collected")
        }
      }else{
        stop("Different output files were found for the path='",path,"'",
            "\n  Remove old files. No output was collected")
      }

      filenames <- normalizePath(paste0(infolder,"/",filenames[order(indexFile)]))
      if(!all(file.exists(filenames))){
        stop("Some files are missing for the path='",path,"'")
      }

      for(i in seq_along(filenames))
      {
        load(filenames[i])
        if(i==1){
          fm <- out

        }else{
          fm$q <- fm$q + out$q
          fm$file_beta <- c(fm$file_beta, out$file_beta)
          fm$tst <- c(fm$tst, out$tst)
          fm$u <- rbind(fm$u, out$u)
          fm$yHat <- rbind(fm$yHat, out$yHat)
          fm$nsup <- rbind(fm$nsup, out$nsup)
          fm$lambda <- rbind(fm$lambda, out$lambda)
          fm$fileID <- c(fm$fileID, out$fileID)
        }
        message(" Loaded file: '",basename(filenames[i]),"'")
      }
      if(length(fm$file_beta)>0){
        stopifnot(length(unique(fm$file_beta))==1L)
        fm$file_beta <- fm$file_beta[1]
      }

  }else{
     stop("No output files were found for the path='",path,"'")
  }
  fm
}
