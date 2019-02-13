
############################################ FastqFilterG ############################################

#' @title Filter your GridION X5 .fastq files
#' @description Use this function if you have .fastq and sequencing summary files and you want to filter your .fastq files in order to obtain only the high-quality ones. FastqFilterG can optionally return a .fasta file for the high-quality sequences. FastqFilterG can optionally return a total .fastq file (equivalent to a 'cat' shell command) and translate into .fasta
#' @param Data Path to .fastq and sequencing summary files returned from GridION X5 (can find .fastq and sequencing summary files recursively)
#' @param DataOut Where the .fastq (and, optionally, .fasta) file will be saved
#' @param FASTQTOT Logical. If TRUE, combine all the .fastq together and store in the DataOut folder. Default to FALSE
#' @param FASTA Logical. If FALSE, return only .fastq file else, if TRUE, return both .fastq and .fasta files. Default to FALSE
#' @param Cores Number of cores to use to accelerate sequencing summary files reading. Default to 1
#' @param Label Llabel to use, together with the Flow Cell identifier, to identify the experiment
#' @return High-quality .fastq file and, optionally, high quality .fasta file, total .fastq file and total .fasta file. If one or more .fastq file is "ill-formatted", FastqFilterG stops and prints the number of the "guilty" .fastq file.
#' @examples
#' #do not run
#' DataPath<-"/data/basecalled/ExperimentName/FlowCellId"
#' FastqFilterG(Data=DataPath, DataOut="Path/To/DataOut",FASTQTOT=FALSE,FASTA=FALSE)


FastqFilterG<-function(Data,DataOut,FASTQTOT=FALSE,FASTA=FALSE,Cores=1,Label) {
  
  library(ShortRead)

  Label<-as.character(Label)
  
  FastqFiles<-list.files(Data,pattern=".fastq",full.names=TRUE, recursive=TRUE)
  FastqFilesPathOrdered<-FastqFiles[order(as.numeric(gsub("[^0-9]+", "", FastqFiles)))]
  
  
  SummariesFiles<-list.files(Data,full.names=TRUE,pattern="sequencing_summary", recursive=TRUE)
  SummariesFilesOrdered<-SummariesFiles[order(as.numeric(gsub("[^0-9]+", "", SummariesFiles)))]
  
  Minimal_Read_Table_Summary<-function(i,File) {
    Table<-read.table(File[i],header=FALSE,sep="\t",skip=1)
    Flowcell_ID_Label<-unlist(strsplit(as.character(Table[,1]),"_"))[4]
    Flowcell_ID<-rep(Flowcell_ID_Label,dim(Table)[1])
    Read_Id<-as.character(Table[,2])
    Qscore<-as.numeric(Table[,12])
    Table<-cbind(Flowcell_ID,Read_Id,Qscore)
    return(Table)
  }
  
  
  cl <- makeCluster(as.numeric(Cores)) 
  clusterExport(cl, c("Minimal_Read_Table_Summary","SummariesFilesOrdered"),envir=environment())
  List<-parLapply(cl, c(1:length(SummariesFilesOrdered)),Minimal_Read_Table_Summary,SummariesFilesOrdered)
  stopCluster(cl)
  
  SummaryTable<-do.call(rbind,List)
  
  Directory<-file.path(DataOut)
  dir.create(Directory, showWarnings=FALSE)
  
  
  TablePass<-which(as.numeric(SummaryTable[,3]) >= 7)
  PassTable<-Table[TablePass,]
  IdPass<-as.character(PassTable[,2])
  Flowcell_ID_Label<-as.character(Table[1,1])
  if (Flowcell_ID_Label == "GA10000") {
    Flowcell_ID_Label<-as.character("FC1")
  }
  if (Flowcell_ID_Label == "GA20000") {
    Flowcell_ID_Label<-as.character("FC2")
  }
  if (Flowcell_ID_Label == "GA30000") {
    Flowcell_ID_Label<-as.character("FC3")
  }
  if (Flowcell_ID_Label == "GA40000") {
    Flowcell_ID_Label<-as.character("FC4")
  }
  if (Flowcell_ID_Label == "GA50000") {
    Flowcell_ID_Label<-as.character("FC5")
  }
  
  message("Starting .fastq files analysis...")
  
  for(i in 1:length(FastqFilesPathOrdered)) {
    
    fqFile<-FastqFile(FastqFilesPathOrdered[i])
    BuildTotalFastq<-tryCatch({
      readFastq(fqFile)},
      error = function(cond) {
        return(NULL)},
      warning = function(cond) {
        message(cond)
        return(NULL)}
    )
    if (is.null(BuildTotalFastq)) {
      stop(paste0("Ill-formatted .fastq file at fastq ",i, "! Can't write .fastq information. This error occurs when a .fastq file is not formatted correctly."))
    }
    else {
      if (FASTQTOT == TRUE) {
        if (i == 1) {
          message("Start writing total .fastq file...")
        }
        writeFastq(BuildTotalFastq, file.path(Directory,paste0(Label,"_",Flowcell_ID_Label,"_FastqTot.fastq")), mode="a", compress=FALSE)
        if (FASTA == TRUE) {
          if (i == 1) {
            
            message("Start writing total .fasta file...")
          }  
          
          writeFasta(BuildTotalFastq,file.path(Directory,paste0(Label,"_",Flowcell_ID_Label,"_FastaTot.fastq")), mode="a", compress=FALSE)
        }
      }
      Id<-id(BuildTotalFastq)
      CharId<-as.character(Id)
      NewReadsSplitted<-unlist(strsplit(CharId," "))
      IdNames<-NewReadsSplitted[seq(1,length(NewReadsSplitted),8)]
      close(fqFile)
      Matches<-match(IdNames,IdPass,nomatch=-10)
      MatchesNew<-which(Matches != -10)
      ReadsPassed<-BuildTotalFastq[MatchesNew]
      if (i == 1){
        message("Start writing high-quality .fastq file...")
      }
      writeFastq(ReadsPassed,file=file.path(Directory, paste0(Label,"_",Flowcell_ID_Label,"_PassFastq.fastq")),compress=FALSE,mode="a")
      if (FASTA==TRUE) {
        if (i == 1) {
          message("Start writing high-quality .fasta file...")
        }
        writeFasta(ReadsPassed,file=file.path(Directory, paste0(Label,"_",Flowcell_ID_Label,"_PassFasta.fasta")),compress=FALSE,mode="a")
      }
    }
  }
  message("Done!")
}

