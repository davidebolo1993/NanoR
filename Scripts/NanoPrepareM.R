 
############################################ NanoPrepareM ############################################


#' @title Prepares MinION data for your analyses with NanoR
#' @description NanoPrepareM generates an object of class list that contains informations required by other functions from NanoR when analyzing MinION data
#' @param DataPass Path to MinION passed .fast5 files folder
#' @param DataFail Path to MinION failes .fast5 files folder
#' @param DataSkip Path to MinION skipped .fast5 files folder
#' @param Label Label to identify your MinION experiment
#' @details NanoPreareM can find .fast5 files recursivel. DataFail and DataSkip can be omitted (MinKNOW generates passed, failes and skipped .fast5 files folders but failed and skipped .fast5 files are taken into account only for calculating their number and percentage)
#' @return Object of class list
#' @examples
#' #do not run
#' PathPass<-"/Path/To/PassFast5"
#' Lab<-"Exp"
#' NanoMList<-NanoPrepareM(DataPass=PathPass, Label=Lab)


NanoPrepareM<-function(DataPass,DataFail=NA,DataSkip=NA,Label) { 
  
  PassFiles<-list.files(DataPass, full.names=TRUE, recursive = TRUE, pattern=".fast5")
  
  message("Found ", length(PassFiles), " .fast5 files specified as passed!")

  if (is.na(DataFail) == FALSE) {
    FailFilesLength<-length(list.files(DataFail, recursive = TRUE, pattern=".fast5"))
    message("Found ", FailFilesLength, " .fast5 files specified as failed!")
  }
  if (is.na(DataFail)) {
    FailFilesLength<-0
    message("No failed .fast5 files path specified.")
  }
  
  if (is.na(DataSkip) == FALSE) {
    SkipFilesLength<-length(list.files(DataSkip, recursive = TRUE, pattern=".fast5"))
    message("Found ", SkipFilesLength, " .fast5 files specified as skipped!")
  }
  if (is.na(DataSkip)) {
    SkipFilesLength<-0
    message("No skipped .fast5 files path specified.")
  }
  
  Label<-as.character(Label)
  
  List<-list(PassFiles,FailFilesLength,SkipFilesLength,Label)
  return(List)
}

