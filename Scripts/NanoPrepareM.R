 
############################################ NanoPrepareM ############################################


#' @title Prepares MinION data for your analyses with NanoR
#' @description NanoPrepareM generates an object of class "list" that contains informations required by other functions from NanoR when analyzing MinION data.
#' @param DataPass Path to MinION "Pass" .fast5 files folder (can find .fast5 files recursively)
#' @param DataFail Path to MinION "Fail" .fast5 files folder (can find .fast5 files recursively)
#' @param DataSkip Path to MinION "Skip" .fast5 files folder (can find .fast5 files recursively)
#' @param Label One label that identifies your MinION experiment
#' @details DataFail and DataSkip can be omitted as, altough MinION local basecaller (MinKNOW) generates "Pass", "Fail" and "Skip" .fast5 files folders, "Fail" and "Skip" .fast5 files are taken into account only for calculating their number and percentage.
#' @return Object of class "list" containing informations required by NanoTableM and NanoStatsM functions.
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

