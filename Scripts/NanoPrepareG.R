

############################################ NanoPrepareG ############################################



#' @title Prepares GridION X5 data for your analyses with NanoR
#' @description NanoPrepareG generates an object of class "list" that contains informations required by other functions from NanoR when analyzing GridION X5 data.
#' @param BasecalledFast5 Logical. Have your experiment returned basecalled .fast5 files (BasecalledFast5=TRUE)? By default GridION X5 outputs .fastq, sequencing summary .txt files and non-basecalled .fast5 files (that's why BasecalledFast5 is, by default, set to FALSE).
#' @param Data Path to the GridION X5 folder that includes .fastq and sequencing summary files (if BasecalledFast5 = FALSE) or to basecalled .fast5 files (if BasecalledFast5 = TRUE)
#' @param DataSkip Path to "Skip" .fast5 files folder (can find .fast5 files recursively); non-basecalled reads (considered "Skipped") are stored in /data/reads/[FlowCellId]/[ExperimentId]/fast5/[FolderNumber]
#' @param DataFail Path to "Fail" .fast5 files folder (can find .fast5 files recursively); can be provided, for example, if working with .fast5 basecalled with ONT's Albacore
#' @param Cores The number of cores to be used to accelerate sequencing summary files reading (useful only if BasecalledFast5 is set to FALSE)
#' @param Label A label that will be used, together with the Flow Cell identifier extracted from the inputted data, to identify your experiment: do not use labels that contains the underscore symbles ("_")
#' @details If working with .fastq files and sequencing summary files, DataSkip can be omitted (in that case the number of skipped data is considered to be 0) while the number of failed reads is automatically computed by NanoR. If working with basecalled .fast5 files, DataFail can be omitted but, in this case, the number of failed .fast5 files is considered to be 0.
#' @return Object of class "list" containing informations required by NanoTableG and NanoStatsG functions.
#' @examples
#' #do not run
#' #when working with sequencing summary files and .fastq files
#' NanoGList<-NanoPrepareG(BasecalledFast5=FALSE, Data="/data/basecalled/ExperimentName/FlowCellId", DataSkip="/data/reads/[FlowCellId]/[ExperimentId]/fast5/" Label="Ex", Cores=3)
#' #when working with basecalled .fast5 files
#' Pass<-"Path/to/workspace/pass"
#' Fail<-"Path/to/workspace/fail"
#' NanoGList<-NanoPrepareG(BasecalledFast5=TRUE, Data=Pass, DataFail=Fail, Label="Ex")



NanoPrepareG<-function(BasecalledFast5=FALSE,Data,DataFail=NA,DataSkip=NA,Cores=1, Label) {
  
  if(BasecalledFast5 == FALSE) {

  	library(parallel)
    
    FastqFiles<-list.files(Data,pattern=".fastq",full.names=TRUE, recursive=TRUE)
    FastqFilesPathOrdered<-FastqFiles[order(as.numeric(gsub("[^0-9]+", "", FastqFiles)))]
    message("Found ",length(FastqFiles), " .fastq files in folder!")
    
    
    SummariesFiles<-list.files(Data,full.names=TRUE,pattern="sequencing_summary", recursive=TRUE)
    SummariesFilesOrdered<-SummariesFiles[order(as.numeric(gsub("[^0-9]+", "", SummariesFiles)))]
    
    message("Found ",length(SummariesFiles), " sequencing summary files in folder!")
    message("Reading and organizing sequencing summary files...")
    
    Read_Table_Summary<-function(i,File) {
      Table<-read.table(File[i],header=FALSE,sep="\t",skip=1)
      RealativeTimeToAdd<-(as.numeric(Table[,8])+as.numeric(Table[,10]))## calculate a relative time that will be rescaled
      SummaryTable<-cbind(Table,RealativeTimeToAdd)
      Flowcell_ID_Label<-unlist(strsplit(as.character(Table[,1]),"_"))[4]
      Flowcell_ID<-rep(Flowcell_ID_Label,dim(Table)[1])
      Read_Id<-as.character(Table[,2])
      Channel<-as.numeric(Table[,4])
      Length<-as.numeric(Table[,11])
      Qscore<-as.numeric(Table[,12])
      Relative_Time<-as.numeric(SummaryTable[,14])
      Table<-cbind(Flowcell_ID,Read_Id,Channel,Relative_Time,Length,Qscore)
      return(Table)
    }
    
    
    cl <- makeCluster(as.numeric(Cores)) 
    clusterExport(cl, c("Read_Table_Summary","SummariesFiles"),envir=environment())
    List<-parLapply(cl, c(1:length(SummariesFiles)),Read_Table_Summary,SummariesFiles)
    stopCluster(cl)
    
    SummaryTable<-do.call(rbind,List)
    
    colnames(SummaryTable)<-c("Flowcell ID","Read Id","Channel Number","Relative Time","Length of Read","Quality")
    Fast5Data<-NA
    DataSkip_Length<-0
    DataFail_Length<-0
    label<-as.character(Label)
    
    
    message("Done!")
  }
  
  if (BasecalledFast5 == TRUE) {
    FastqFilesPathOrdered<-NA
    SummaryTable<-NA
    Fast5Data<-list.files(Data,recursive=TRUE,full.names=TRUE,pattern=".fast5")
    message("Found ",length(Fast5Data), " .fast5 files in folder!")
    if(is.na(DataFail)) {
      DataFail_Length<-0
      message("No failed .fast5 files folder specified.")
    }
    else {
      DataFail_Length<-length(list.files(DataFail,recursive=TRUE,full.names=TRUE,pattern=".fast5"))
      message("Found ",DataFail_Length, " failed .fast5 files in folder!")
    }
    if (is.na(DataSkip)) {
      DataSkip_Length<-0
      message("No skipped .fast5 files folder specified.")
    }
    else {
      DataSkip_Length<-length(list.files(DataSkip,recursive=TRUE,full.names=TRUE,pattern=".fast5"))
      message("Found ",length(SkippedData), " skipped .fast5 files in folder!")
    }
    label<-as.character(Label)
    message("Done!")
  }
  
  List<-list(FastqFilesPathOrdered,SummaryTable,Fast5Data,DataFail_Length,DataSkip_Length,label)
  return(List)
  
}

