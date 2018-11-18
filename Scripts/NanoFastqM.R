
############################################ NanoFastqM ############################################

#' @title Extracts .fastq informations from your MinION "Pass" .fast5 files
#' @description NanoFastqM returns a file text containing .fastq sequences (and, if you want, an additional file text containing .fasta sequences) extracted from your "Pass" .fsat5 reads.
#' @param DataPass Path to "Pass" .fast5 files folder (can find .fast5 files recursively)
#' @param DataOut Where your .fastq (and .fasta) file will be saved
#' @param Label A lable used to identify .fastq (and .fasta) file outputted
#' @param Cores Number of cores to be used: 1 by default
#' @param FASTA Logical. If FALSE (by default), return only .fastq file else, if TRUE, return both .fastq and .fasta files
#' @return .fastq file and, if you wish, .fasta file for your "Pass" .fast5 files.
#' @examples
#' #do not run
#' NanoFastqM(DataPass="Path/To/DataPass", DataOut="/Path/To/DataOutExp", Cores=6, FASTA=FALSE)
#' NanoFastqM(DataPass="Path/To/DataPass", DataOut="/Path/To/DataOutExp", Cores=6, FASTA=TRUE)


NanoFastqM<-function(DataPass,DataOut,Label,Cores=1,FASTA=FALSE) {
  
  library(rhdf5)
  library(parallel)
  library(ShortRead)
  
  Directory<-file.path(DataOut)
  dir.create(Directory,showWarnings=FALSE)
  setwd(Directory)
  
  label<-as.character(Label)
  PassFiles<-list.files(DataPass, full.names=TRUE, recursive = TRUE, pattern=".fast5")
  
  
  ### FUNCTIONS ###
  
  Read_DataSet<-function(File, PathGroup) { 
    h5errorHandling(type="suppress")
    Data1<-H5Dopen(File, PathGroup) 
    Data2<-H5Dread(Data1)
    H5Dclose(Data1)
    return(Data2) 
  }
  
  Read_Attributes<-function(PathGroup, Attribute) { 
    h5errorHandling(type="suppress")
    Data1<-H5Aopen(PathGroup, Attribute)
    Data2<-H5Aread(Data1)
    H5Aclose(Data1)
    return(Data2) 
  } 
  
  Fastq_Extraction<-function(i,File) {
    
    
    Fastq<-list()
    
    h5errorHandling(type="suppress")
    File<-H5Fopen(File[i])
    
    GroupTry<-"/Analyses/Basecall_1D_000/BaseCalled_template/Fastq"
    
    Try<-try(Read_DataSet(File,GroupTry), silent=TRUE) #exclude not-basecalled .fast5 files
    
    if (inherits(Try,"try-error")) {
      H5Fclose(File)
      Fastq[[i]]<-NA
      return(Fastq[[i]])
    }
    
    else {
      
      Group0<-"/Analyses/Basecall_1D_000/Summary/basecall_1d_template"
      Score<-H5Gopen(File,Group0)
      Quality<-Read_Attributes(Score,"mean_qscore")
      H5Gclose(Score)
      
      if (Quality >= 7) { 
        Pre_Fastq<-Read_DataSet(File,GroupTry)
        Fastq[[i]]<- strsplit(Pre_Fastq,split="\n",fixed=TRUE)[[1]]
      }
      
      else
      {
        Fastq[[i]]<-NA
      }
      H5Fclose(File)
      return(Fastq[[i]])
    }
  }
  
  ################################
  
  
  if (FASTA == FALSE) {
    message("Extracting .fastq sequences...")
  }
  else {
    message("Extracting .fastq and .fasta sequences...")
  }
  
  cl <- makeCluster(as.numeric(Cores)) 
  clusterExport(cl, c("Fastq_Extraction","PassFiles","Read_DataSet","Read_Attributes"),envir=environment())
  clusterEvalQ(cl,library(rhdf5))
  List<-parLapply(cl, c(1:length(PassFiles)),Fastq_Extraction,PassFiles)
  stopCluster(cl)
  FastqTot<-do.call(c,List)
  FastqClean<-which(is.na(FastqTot) == FALSE)
  FastqFinal<-(FastqTot[FastqClean])
  message("Writing .fastq file!")
  fileConn<-file(paste0(label,".fastq"))
  writeLines(FastqFinal,fileConn)
  close(fileConn)
  
  
  if (FASTA == TRUE) {
    message("Writing .fasta file!")
    Fastq = FastqStreamer(file.path(Directory,paste0(label,".fastq")))
    repeat {
      Fasta = yield(Fastq)
      if (length(Fasta) == 0) break
      writeFasta(Fasta, file=file.path(Directory,paste0(label,".fasta")), mode="a")
    }
  }
  message("Done!") 
}

