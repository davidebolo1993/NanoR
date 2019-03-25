############################################ NanoPrepareM ############################################


#' @title Prepares basecalled data 
#' @description NanoPrepareM prepares MinION and GridION X5 basecalled data for other functions from this package. Name of the function is inherited from NanoR previous version as MinION and GridION X5 had different default-output formats and only MinION outputted basecalled .fast5 files. From MinION release 18.12 and GridION 18.12.1 outputs are the same.
#' @param DataPass Path to passed .fast5 files folder
#' @param DataFail Path to failed .fast5 files folder. Default to NA
#' @param DataSkip Path to skipped .fast5 files folder. Default to NA
#' @param Label Label to identify the experiment. A folder with this name will be created in DataOut directory
#' @param MultiRead Logical. If TRUE, enable multi-read .fast5 files support. Default to FALSE
#' @details NanoPreareM can find .fast5 files recursively, so be careful to specify path to folder contaning a unique data type (e.g. only passed .fast5 files for the experiment). DataFail and DataSkip can be omitted.
#' @return Object of class list
#' @examples
#' #do not run
#' DataPass<-"/path/to/fast5_pass"
#' DataFail<-"/path/to/fast5_fail" #can be omitted
#' # PathSkip. Useful for old MinION data versions 
#' Label<-"Exp"
#' #single-read .fast5 files
#' List<-NanoPrepareM(DataPass,DataFail,Label=Label)
#' #multi-read .fast5 files
#' List<-NanoPrepareM(DataPass,DataFail, Label=Label,MultiRead=TRUE)



NanoPrepareM<-function(DataPass,DataFail=NA,DataSkip=NA, Label, MultiRead=FALSE) { #simple and fast way to store informations
  
  PassFiles<-list.files(DataPass, full.names=TRUE, recursive = TRUE, pattern=".fast5")
  
  if (MultiRead == FALSE) {

    message(length(PassFiles), " .fast5 files specified as passed")

  }

  else {

    message(length(PassFiles), " multi-read .fast5 files specified as passed")

  }

  if (is.na(DataFail)) {
    FailFilesLength<-0
    message("No failed .fast5 files path specified")
  }

  else {

    if (MultiRead==FALSE) {

      FailFilesLength<-length(list.files(DataFail, recursive = TRUE, pattern=".fast5"))
      message(FailFilesLength, " .fast5 files specified as failed")
    }

    else {

      library(rhdf5)

      FailFilesLength<-0
      FailFiles<-list.files(DataFail, recursive = TRUE, pattern=".fast5", full.names=TRUE)
      message(length(FailFiles), " multi-read .fast5 files specified as failed")


      FailFilesOrdered<-FailFiles[order(as.numeric(gsub("[^0-9]+", "", FailFiles)))]

      First<-FailFilesOrdered[1]
      FileOpen<-H5Fopen(First)
      howmany<-nrow(h5ls(FileOpen, recursive=FALSE, datasetinfo=FALSE))
      H5Fclose(FileOpen)

      Last<-FailFilesOrdered[length(FailFilesOrdered)]
      FileOpen<-H5Fopen(Last)
      howmany_<--nrow(h5ls(FileOpen, recursive=FALSE, datasetinfo=FALSE))
      H5Fclose(FileOpen)

      FailFilesLength<-(howmany*(length(FailFiles)-1) + howmany_)

    }


  }
  
  if (is.na(DataSkip)) { #no data skip for multiline? 
    SkipFilesLength<-0
    message("No skipped .fast5 files path specified")
  }

  else {

    if (MultiRead==FALSE) {

      SkipFilesLength<-length(list.files(DataSkip, recursive = TRUE, pattern=".fast5"))
      message(SkipFilesLength, " .fast5 files specified as skipped")

    }

    else { ## don't think is needed

      library(rhdf5)

      SkipFilesLength<-0
      SkipFiles<-list.files(DataSkip, recursive = TRUE, pattern=".fast5", full.names=TRUE)
      message(length(SkipFiles), " multi-read .fast5 files specified as skipped")


      SkipFilesOrdered<-SkipFiles[order(as.numeric(gsub("[^0-9]+", "", SkipFiles)))]

      First<-SkipFilesOrdered[1]
      FileOpen<-H5Fopen(First)
      howmany<-nrow(h5ls(FileOpen, recursive=FALSE, datasetinfo=FALSE))
      H5Fclose(FileOpen)

      Last<-SkipFilesOrdered[length(SkipFilesOrdered)]
      FileOpen<-H5Fopen(Last)
      howmany_<--nrow(h5ls(FileOpen, recursive=FALSE, datasetinfo=FALSE))
      H5Fclose(FileOpen)

      SkipFilesLength<-(howmany*(length(SkipFiles)-1) + howmany_)

    }

  }

  Label<-as.character(Label)
  
  List<-list(PassFiles,FailFilesLength,SkipFilesLength,Label,MultiRead)

  message("Done")  

  return(List)

}