
############################################ NanoFastqG ############################################

#' @title Filter .fastq files
#' @description Retain only high-quality .fastq from .fastq files (optionally, converts .fastq to .fasta as well).
#' @param DataSummary Path to sequencing summary file/s folder
#' @param DataFastq Path to passed .fastq folder
#' @param DataOut Where the .fastq file will be saved
#' @param Label Label to identify the experiment. A folder with this name will be created in DataOut directory.
#' @param Cores Number of cores to use to accelerate sequencing summary files reading when dealing with multiple ones. Default to 1
##' @param FASTQTOT Logical. If TRUE, combine all the .fastq together and store in the DataOut folder. Default to FALSE. #removed because there are faster way to concatenate. In bash, use "cat".
#' @param FASTA Logical. If TRUE translate .fastq to .fasta as well. Default to FALSE
#' @param Minquality Minimum quality to retain the .fastq sequence. Default to 7
#' @return filtered .fastq file and, optionally, convert to .fastq file
#' @examples
#' #do not run
#' new behaviour
#' DataSummary<-"/path/to/sequencing_summary"
#' DataFastq<-"/path/to/fastq_pass"
#' DataOut <- "/path/to/DataOut"
#' Label<-'Exp'
#' Filter on higher treshold
#' NanoFastqG(DataSummary, DataFastq, DataOut, Label, Minquality=10)
#' At the moment, NanoFastqM do not offer gz compression, as is it very slow to do in the R environment.


NanoFastqG<-function(DataSummary, DataFastq, DataOut,Label,Cores=1,FASTA=FALSE, Minquality=7) {
  
  library(ShortRead)

  Label<-as.character(Label)

  FastqFiles<-list.files(DataFastq,pattern=".fastq",full.names=TRUE, recursive=TRUE)
  SummariesFiles<-list.files(DataSummary,full.names=TRUE,pattern="sequencing_summary", recursive=TRUE)
  
  Minimal_Read_Table_Summary<-function(i,File) {
    Table<-read.table(File[i],header=FALSE,sep="\t",skip=1)
    #Flowcell_ID_Label<-unlist(strsplit(as.character(Table[,1]),"_"))[4]
    #Flowcell_ID<-rep(Flowcell_ID_Label,dim(Table)[1])
    Read_Id<-as.character(Table[,2])
    Qscore<-as.numeric(Table[,12])
    Table<-cbind(Read_Id,Qscore)
    return(Table)
  }

  Minimal_Read_Table_Summary_SC<-function(File) {
    Table<-read.table(File,header=FALSE,sep="\t",skip=1)
    #Flowcell_ID_Label<-unlist(strsplit(as.character(Table[,1]),"_"))[4]
    #Flowcell_ID<-rep(Flowcell_ID_Label,dim(Table)[1])
    Read_Id<-as.character(Table[,2])
    Qscore<-as.numeric(Table[,12])
    Table<-cbind(Read_Id,Qscore)
    return(Table)
  }

  Minimal_Read_Table_New<-function(File) {
    Table<-read.table(File,header=FALSE,sep="\t",skip=1)
    #Flowcell_ID_Label<-unlist(strsplit(as.character(Table[,1]),"_"))[4]
    #Flowcell_ID<-rep(Flowcell_ID_Label,dim(Table)[1])
    Read_Id<-as.character(Table[,3])
    Qscore<-as.numeric(Table[,15])
    Table<-cbind(Read_Id,Qscore)
    return(Table)
  }


  if (length(SummariesFiles) == 1) {

    SummaryTable<-Minimal_Read_Table_New(SummariesFiles) #if only one table, assume is the new format

  }

  else { #old behaviour

    if (Cores > 1) {
    
      cl <- makeCluster(as.numeric(Cores)) 
      clusterExport(cl, c("Minimal_Read_Table_Summary","SummariesFiles"),envir=environment())
      List<-parLapply(cl, c(1:length(SummariesFiles)),Minimal_Read_Table_Summary,SummariesFiles)
      stopCluster(cl)

    }

    else {

      List <- lapply(SummariesFiles,Minimal_Read_Table_Summary_SC)

    }

    SummaryTable<-do.call(rbind,List)

  }
    
  
  #Flowcell_ID_Label<-as.character(SummaryTable[1,1])
  
  #if (Flowcell_ID_Label == "GA10000") {
    #Flowcell_ID_Label<-as.character("FC1")
  #}
  #if (Flowcell_ID_Label == "GA20000") {
    #Flowcell_ID_Label<-as.character("FC2")
  #}
  #if (Flowcell_ID_Label == "GA30000") {
    #Flowcell_ID_Label<-as.character("FC3")
  #}
  #if (Flowcell_ID_Label == "GA40000") {
    #Flowcell_ID_Label<-as.character("FC4")
  #}
  #if (Flowcell_ID_Label == "GA50000") {
    #Flowcell_ID_Label<-as.character("FC5")
  #}

  Directory<-file.path(DataOut,Label)
  dir.create(Directory, showWarnings=FALSE)

  TablePass<-which(as.numeric(SummaryTable[,2]) >= Minquality) #higher the treshold, lower the ER in final bamfiles!
  IdPass<-as.character(SummaryTable[TablePass,][,1])
  
  message("Filtering .fastq files...")
  
  for(i in 1:length(FastqFiles)) {
    
    fqFile<-FastqFile(FastqFiles[i])
    
    BuildTotalFastq<-tryCatch({
      readFastq(fqFile)},
      error = function(cond) {
        return(NULL)},
      warning = function(cond) {
        message(cond)
        return(NULL)}
    )

    if (is.null(BuildTotalFastq)) {
        
      warning("Ill-formatted .fastq file: ", FastqFiles[i], " .Skipped")
    
    }
    else {

      Id<-id(BuildTotalFastq)
      CharId<-as.character(Id)
      NewReadsSplitted<-unlist(strsplit(CharId," "))
      IdNames<-NewReadsSplitted[seq(1,length(NewReadsSplitted),6)] #get identifier. Suppose the read identifier position is costant
      close(fqFile)
      Matches<-match(IdNames,IdPass,nomatch=-10)
      MatchesNew<-which(Matches != -10)
      ReadsPassed<-BuildTotalFastq[MatchesNew]
      writeFastq(ReadsPassed,file=file.path(Directory, paste0(Label,".fq")),compress=FALSE,mode="a")
      if (FASTA==TRUE) {      
        writeFasta(ReadsPassed,file=file.path(Directory, paste0(Label,".fa")),compress=FALSE,mode="a")      
      }
    }
  }

  message("Done")

}


