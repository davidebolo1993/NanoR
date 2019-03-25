
############################################ NanoFastqM ############################################

#' @title Extracts .fastq informations
#' @description NanoPrepareM prepares MinION and GridION X5 basecalled data for other functions from this package. Name of the function is inherited from NanoR previous version as MinION and GridION X5 had different default-output formats and only MinION outputted basecalled .fast5 files. From MinION release 18.12 and GridION 18.12.1 outputs are the same. This is useful expecially to filter at different treshold than the default one.
#' @param DataPass Path to passed .fast5 files folder
#' @param DataOut Where the .fastq file will be saved
#' @param Label Label to identify the MinION experiment. A folder with this name will be created in DataOut directory.
#' @param Cores Number of cores to be used: 1 by default
#' @param FASTA Logical. If TRUE translate .fastq to .fasta as well. Default to FALSE
#' @param Minquality Minimum quality to retain the .fastq sequence. Default to 7
#' @param MultiRead Logical. If TRUE, enable multiread .fast5 files support. Default to FALSE. 
#' @return .fastq file and, optionally, .fasta file for MinION passed .fast5 files
#' @examples
#' #do not run
#' DataPass<-"/path/to/fast5_pass"
#' DataOut <- "/path/to/DataOut"
#' Label<-'Exp'
#' #single-read .fast5 files
#' #Extract
#' NanoFastqM(DataPass, DataOut, Label, Cores=6) #Minquality=7
#' #Extract, convert to .fastq and filter.
#' NanoFastqM(DataPass, DataOut, Label, Cores=6, FASTA=TRUE, Minquality=10)
#' Extract from multi-read .fast5 files
#' NanoFastqM(DataPass, DataOut, Cores=6, MultiRead=TRUE)
#' At the moment, NanoFastqM do not offer gz compression, as is it very slow to do in the R environment.


NanoFastqM<-function(DataPass,DataOut,Label,Cores=1,FASTA=FALSE, Minquality=7,MultiRead=FALSE) {
  
  library(rhdf5)
  library(parallel)
  library(ShortRead)
  
  label<-as.character(Label)
  Directory<-file.path(DataOut, label)
  dir.create(Directory,showWarnings=FALSE, recursive=TRUE)
  setwd(Directory)
  
  label<-as.character(Label)

  PassFiles<-list.files(DataPass, full.names=TRUE, recursive = TRUE, pattern=".fast5")
  
  if (MultiRead ==FALSE) {

    message(length(PassFiles), " .fast5 files specified as passed")

  }

  else {

    message(length(PassFiles), " multiread .fast5 files specified as passed")

  }
  
  
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
    
        
    h5errorHandling(type="suppress")
    File<-H5Fopen(File[i])
    
    GroupTry<-"/Analyses/Basecall_1D_000/BaseCalled_template/Fastq"
    
    Try<-try(Read_DataSet(File,GroupTry), silent=TRUE) #exclude not-basecalled .fast5 files
    
    if (inherits(Try,"try-error")) {
      H5Fclose(File)
    }
    
    else {
      
      Group0<-"/Analyses/Basecall_1D_000/Summary/basecall_1d_template"
      Score<-H5Gopen(File,Group0)
      Quality<-Read_Attributes(Score,"mean_qscore")
      H5Gclose(Score)
      
      if (Quality >= as.numeric(Minquality)) { 
        Pre_Fastq<-Read_DataSet(File,GroupTry)
        Fastq<- strsplit(Pre_Fastq,split="\n",fixed=TRUE)[[1]]
        H5Fclose(File)
        return(Fastq)
      } 
      H5Fclose(File)
    }
  }

  Fastq_Extraction_SC<-function(File) { ## single core
    
        
    h5errorHandling(type="suppress")
    File<-H5Fopen(File)
    
    GroupTry<-"/Analyses/Basecall_1D_000/BaseCalled_template/Fastq"
    
    Try<-try(Read_DataSet(File,GroupTry), silent=TRUE) #exclude not-basecalled .fast5 files
    
    if (inherits(Try,"try-error")) {
      H5Fclose(File)

    }
    
    else {
      
      Group0<-"/Analyses/Basecall_1D_000/Summary/basecall_1d_template"
      Score<-H5Gopen(File,Group0)
      Quality<-Read_Attributes(Score,"mean_qscore")
      H5Gclose(Score)
      
      if (Quality >= as.numeric(Minquality)) { 
        Pre_Fastq<-Read_DataSet(File,GroupTry)
        Fastq<- strsplit(Pre_Fastq,split="\n",fixed=TRUE)[[1]]
        H5Fclose(File)
        return(Fastq)

      }
      
      H5Fclose(File)      

    }
  }

  ####################### NEW FUNCTIONS FOR MULTIREAD .fast5 FILES #########################


  Fastq_Extraction_Multiline<-function(i,File) {
    
        
    h5errorHandling(type="suppress")
    File<-H5Fopen(File[i])
    idtab<-h5ls(File, recursive=FALSE, datasetinfo=FALSE)[2]
    List<-rep(list(NA),nrow(idtab))

    for (i in 1:nrow(idtab)) {

      GroupQuality<-H5Gopen(File,paste0(idtab[i,],'/Analyses/Basecall_1D_000/Summary/basecall_1d_template'))
      QScore<-Read_Attributes(GroupQuality, 'mean_qscore')
      H5Gclose(GroupQuality)   

      if (QScore >= Minquality) {

        GroupAnalyses<-try(H5Gopen(File,paste0(idtab[i,],'/Analyses/Basecall_1D_000/BaseCalled_template')), silent=TRUE) #check if fastq exists. It is supposed to exists both for failed and passed
        if (inherits(GroupAnalyses,"try-error")) {
          H5Gclose(GroupAnalyses)
          List[[i]]<-NA
        }

        else {

          Pre_Fastq<-try(Read_DataSet(GroupAnalyses, 'Fastq'), silent=TRUE)

          if (inherits(Pre_Fastq,"try-error")) {
            H5Gclose(GroupAnalyses)     
            List[[i]]<-NA
          }

          else {

            List[[i]]<- strsplit(Pre_Fastq,split="\n",fixed=TRUE)[[1]]

          }
        }
      }

      else {

        List[[i]]<-NA

      }
    }

    H5Fclose(File)
    Res<-do.call(c,List[!is.na(List)])
    
    if (length(Res) > 0) {

      return(Res)

    }

  }



  Fastq_Extraction_Multiline_SC<-function(File) {
    
        
    h5errorHandling(type="suppress")
    File<-H5Fopen(File)
    idtab<-h5ls(File, recursive=FALSE, datasetinfo=FALSE)[2]
    List<-rep(list(NA),nrow(idtab))

    for (i in 1:nrow(idtab)) {

      GroupQuality<-H5Gopen(File,paste0(idtab[i,],'/Analyses/Basecall_1D_000/Summary/basecall_1d_template'))
      QScore<-Read_Attributes(GroupQuality, 'mean_qscore')
      H5Gclose(GroupQuality)   

      if (QScore >= Minquality) {

        GroupAnalyses<-try(H5Gopen(File,paste0(idtab[i,],'/Analyses/Basecall_1D_000/BaseCalled_template')), silent=TRUE) #check if fastq exists. It is supposed to exists both for failed and passed
        if (inherits(GroupAnalyses,"try-error")) {
          H5Gclose(GroupAnalyses)
          List[[i]]<-NA
        }

        else {

          Pre_Fastq<-try(Read_DataSet(GroupAnalyses, 'Fastq'), silent=TRUE)

          if (inherits(Pre_Fastq,"try-error")) {
            H5Gclose(GroupAnalyses)     
            List[[i]]<-NA
          }

          else {

            List[[i]]<- strsplit(Pre_Fastq,split="\n",fixed=TRUE)[[1]]

          }
        }
      }

      else {

        List[[i]]<-NA

      }
    }

    H5Fclose(File)
    Res<-do.call(c,List[!is.na(List)])
    
    if (length(Res) > 0) {

      return(Res)

    }
  }


  ################################

  message("Extracting .fastq sequences...")


  if (MultiRead == FALSE) {


    if (Cores >1) {

      cl <- makeCluster(as.numeric(Cores)) 
      clusterExport(cl, c("Fastq_Extraction","PassFiles","Read_DataSet","Read_Attributes", "Minquality"),envir=environment())
      clusterEvalQ(cl,library(rhdf5))
      List<-parLapply(cl, c(1:length(PassFiles)),Fastq_Extraction,PassFiles)
      stopCluster(cl)

    }

    else {

      List<-lapply(PassFiles,Fastq_Extraction_SC) 

    }

  }

  else { #multi read is TRUE


    if (Cores >1) {

      cl <- makeCluster(as.numeric(Cores)) 
      clusterExport(cl, c("Fastq_Extraction_Multiline","PassFiles","Read_DataSet","Read_Attributes", "Minquality"),envir=environment())
      clusterEvalQ(cl,library(rhdf5))
      List<-parLapply(cl, c(1:length(PassFiles)),Fastq_Extraction_Multiline,PassFiles)
      stopCluster(cl)

    }

    else {

      List<-lapply(PassFiles,Fastq_Extraction_Multiline_SC) 

    }

  }



  FastqTot<-do.call(c,List)
  #FastqClean<-which(is.na(FastqTot) == FALSE) #do not spend time, just skip before filtering
  #FastqFinal<-(FastqTot[FastqClean])
  message("Writing .fastq file...")


  fileConn<-file(paste0(label,".fq"))
  writeLines(FastqTot,fileConn)
  close(fileConn)  

  #stream fasta in small chunck, avoid to load all in RAM

  if (FASTA == TRUE) {
    message("Writing .fasta file...")
    Fastq = FastqStreamer(file.path(Directory,paste0(label,".fq")))
    repeat {
      Fasta = yield(Fastq)
      if (length(Fasta) == 0) break
      writeFasta(Fasta, file=file.path(Directory,paste0(label,".fa")), mode="a")
    }
    close(Fastq)
  }

  message("Done")

}
  