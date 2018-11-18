

############################################ NanoTableM ############################################


#' @title Generates an information table for your MinION .fast5 files
#' @description NanoTableM generates a table that contains useful informations for each read of the "Pass" .fast5 files folder given to NanoPrepareM function. As NanoTableM can take some time (it depends on the number of reads you obtain from your experiment), this function can be accelerated using multiple cores.
#' @param NanoPrepareMList The object of class "list" returned by NanoPrepareM function
#' @param DataOut Where the table will be saved. Don't use a directory that already contains another NanoTableM (or NanoTableG) result
#' @param Cores Number of cores to be used: 1 by default
#' @param GCC Logical. If TRUE (by default), NanoTableM computes GC content for each read
#' @return Table containing informations required by NanoStatsM and NanoCompare (and useful for users who want to analyze these data by theirselves!).
#' @examples
#' #do not run
#' NanoMTable<-NanoTableM(NanoPrepareMList=NanoMList, DataOut="/Path/To/DataOutExp",Cores=6,GCC=TRUE)


NanoTableM<-function(NanoPrepareMList,DataOut,Cores=1,GCC=TRUE) {
  
  library(parallel)
  library(rhdf5)
  library(seqinr)
  
  label<-as.character(NanoPrepareMList[[4]])
  Directory<-file.path(DataOut)
  dir.create(Directory, showWarnings = FALSE)
  setwd(Directory)
  TableInDirectory<-list.files(Directory,pattern="Information_Table.txt")
  if(length(TableInDirectory) != 0) {
    stop("Can't use a directory that already contains other NanoTableM/NanoTableG results")
  }
  
  ##### FUNCTIONS ######
  
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
  
  HDF5_File_Parsing_Table_With_GC<-function(i,File) {  ## work for R9.4 and R9.5
    
    h5errorHandling(type="suppress")
    Table<-c(Read_Id="Read_Id",Channel="Channel",Mux="Mux",Unix_Time="Unix_Time",Length="Length",Qscore="Qscore",GC_Content="GC_Content")
    
    File<-H5Fopen(File[i])
    
    Group1<-"/Analyses/Basecall_1D_000/BaseCalled_template/Fastq"
    Try<-try(Read_DataSet(File,Group1), silent=TRUE) #exclude non-basecalled .fast5 reads (will be counted as failed by NanoStats)
    if (inherits(Try,"try-error")) {
      return(Table)
    }
    
    else {
      
      Pre_Fastq<-Read_DataSet(File,Group1)
      Sequence_Fastq<-unlist(strsplit(Pre_Fastq,"\n"))[2]
      Fastq<-s2c(Sequence_Fastq)
      Table['GC_Content']<-GC(Fastq)
      
      Group2<-"/UniqueGlobalKey/channel_id"
      Chann_File<-H5Gopen(File,Group2)
      Table['Channel']<-Read_Attributes(Chann_File,"channel_number")
      H5Gclose(Chann_File)
      
      Group3<-"/Analyses/Basecall_1D_000/Summary/basecall_1d_template"
      Score_Length<-H5Gopen(File,Group3)
      Table['Qscore']<-Read_Attributes(Score_Length,"mean_qscore")
      Table['Length']<-Read_Attributes(Score_Length,"sequence_length")
      H5Gclose(Score_Length)
      
      #Group4<-"/Analyses/Basecall_1D_000/" ###No more used as sometimes it can have months in english word instead of numbers
      #Time<-H5Gopen(File,Group4)
      #if (H5Aexists(Time,"time_stamp")) {
        #Date<-Read_Attributes(Time,"time_stamp")
        #H5Gclose(Time)
        #Y_M_D<-substr(Date,1,10)
        #H_M_S<-substr(Date,12,19)
        #Time_Vector<-paste(c(Y_M_D,H_M_S), collapse=" ")
        #Table['Unix_Time']<-as.numeric(as.POSIXct(strptime(Time_Vector, "%Y-%m-%d %H:%M:%S")))
      #}
      #else {
        #H5Gclose(Time)
        #Table['Unix_Time']<-"Unix_Time"
      #}
      
      Group5<-"/Raw/Reads"
      Read_Mux<-H5Gopen(File,Group5)
      Template <- h5ls(Read_Mux,recursive=FALSE,datasetinfo=FALSE)
      H5Gclose(Read_Mux)
      Read_Path <- paste0("/Raw/Reads/",Template$name)
      
      Group6<- H5Gopen(File,Read_Path)
      
      Table['Mux']<-Read_Attributes(Group6,"start_mux")
      Table['Read_Id']<-Read_Attributes(Group6,"read_id")
      Start<-Read_Attributes(Group6,"start_time")
      AlternativeStart<-floor(Start/4000) ##actual sampling rate
      H5Gclose(Group6)
      
      ## try to compute a relative time using different parameters ##
      
     
      Group4.2<-"/UniqueGlobalKey/tracking_id"
      Time2<-H5Gopen(File,Group4.2)
      DateUnix2<-Read_Attributes(Time2,"exp_start_time")
      H5Gclose(Time2)
      if (length(unlist(strsplit(DateUnix2,"T")))==2) {### avoid problems with UNIX exp start times found on EGA samples
        Y_M_D<-substr(DateUnix2,1,10)
        H_M_S<-substr(DateUnix2,12,19)
        Time_Vector<-paste(c(Y_M_D,H_M_S), collapse=" ")
        Table['Unix_Time']<-as.numeric(as.POSIXct(strptime(Time_Vector, "%Y-%m-%d %H:%M:%S")))+AlternativeStart
      
      }
      else {
        Table['Unix_Time']<-(as.numeric(AlternativeStart)+as.numeric(DateUnix2))
      }
      
      
      H5Fclose(File)
      
      return(Table) 
    }
  }
  
  HDF5_File_Parsing_Table_Without_GC<-function(i,File) {  ## work for R9.4 and R9.5
    
    h5errorHandling(type="suppress")
    Table<-c(Read_Id="Read_Id",Channel="Channel",Mux="Mux",Unix_Time="Unix_Time",Length="Length",Qscore="Qscore",GC_Content="GC_Content")
    
    File<-H5Fopen(File[i])
    
    Group1<-"/Analyses/Basecall_1D_000/BaseCalled_template/Fastq"
    Try<-try(Read_DataSet(File,Group1), silent=TRUE) #exclude not-basecalled .fast5 files (will be counted as failed by NanoStats)
    if (inherits(Try,"try-error")) {
      return(Table)
    }
    
    else {
      
      
      Table['GC_Content']<-'GC_Content'
      
      Group2<-"/UniqueGlobalKey/channel_id"
      Chann_File<-H5Gopen(File,Group2)
      Table['Channel']<-Read_Attributes(Chann_File,"channel_number")
      H5Gclose(Chann_File)
      
      Group3<-"/Analyses/Basecall_1D_000/Summary/basecall_1d_template"
      Score_Length<-H5Gopen(File,Group3)
      Table['Qscore']<-Read_Attributes(Score_Length,"mean_qscore")
      Table['Length']<-Read_Attributes(Score_Length,"sequence_length")
      H5Gclose(Score_Length)
      
      #Group4<-"/Analyses/Basecall_1D_000/" ###No more used as sometimes it can have months in english word instead of numbers
      #Time<-H5Gopen(File,Group4)
      #if (H5Aexists(Time,"time_stamp")) {
        #Date<-Read_Attributes(Time,"time_stamp")
        #H5Gclose(Time)
        #Y_M_D<-substr(Date,1,10)
        #H_M_S<-substr(Date,12,19)
        #Time_Vector<-paste(c(Y_M_D,H_M_S), collapse=" ")
        #Table['Unix_Time']<-as.numeric(as.POSIXct(strptime(Time_Vector, "%Y-%m-%d %H:%M:%S")))
      #}
      #else {
        #H5Gclose(Time)
        #Table['Unix_Time']<-"Unix_Time"
      #}
      
      Group5<-"/Raw/Reads"
      Read_Mux<-H5Gopen(File,Group5)
      Template <- h5ls(Read_Mux,recursive=FALSE,datasetinfo=FALSE)
      H5Gclose(Read_Mux)
      Read_Path <- paste0("/Raw/Reads/",Template$name)
      
      Group6<- H5Gopen(File,Read_Path)
      
      Table['Mux']<-Read_Attributes(Group6,"start_mux")
      Table['Read_Id']<-Read_Attributes(Group6,"read_id")
      Start<-Read_Attributes(Group6,"start_time")
      AlternativeStart<-floor(Start/4000) ##actual sampling rate
      H5Gclose(Group6)
      
      ## try to compute a relative time using different parameters ##
      
      
      Group4.2<-"/UniqueGlobalKey/tracking_id"
      Time2<-H5Gopen(File,Group4.2)
      DateUnix2<-Read_Attributes(Time2,"exp_start_time")
      H5Gclose(Time2)
      if (length(unlist(strsplit(DateUnix2,"T")))==2) { ### avoid problems with UNIX exp start times found on EGA samples
        Y_M_D<-substr(DateUnix2,1,10)
        H_M_S<-substr(DateUnix2,12,19)
        Time_Vector<-paste(c(Y_M_D,H_M_S), collapse=" ")
        Table['Unix_Time']<-as.numeric(as.POSIXct(strptime(Time_Vector, "%Y-%m-%d %H:%M:%S")))+AlternativeStart
          
      }
      else {
        Table['Unix_Time']<-(as.numeric(AlternativeStart)+as.numeric(DateUnix2))
      }
        
      
      H5Fclose(File)
      
      return(Table) 
    }
  }
  
  ######################
  
  PassFiles<-NanoPrepareMList[[1]]
  
  ### CALCULATIING ###
  
  if (GCC == TRUE) {
    message("Extracting metadata from .fast5 files and calculating GC content!")
    cl <- makeCluster(as.numeric(Cores)) 
    clusterExport(cl, c("HDF5_File_Parsing_Table_With_GC","PassFiles","Read_DataSet","Read_Attributes"),envir=environment())
    clusterEvalQ(cl,c(library(rhdf5), library(seqinr)))
    List<-parLapply(cl, c(1:length(PassFiles)), HDF5_File_Parsing_Table_With_GC,PassFiles)
    stopCluster(cl)
    Table_Tot<-data.frame(matrix(unlist(List), nrow=length(List), ncol=7, byrow=TRUE),stringsAsFactors=FALSE)
    colnames(Table_Tot)<-c("Read Id", "Channel Number", "Mux Number", "Unix Time", "Length of Read", "Quality", "GC_Content")
    write.table(Table_Tot, file.path(Directory, paste0(label, "_Information_Table.txt")), col.names=T, row.names=F, quote=F, sep="\t")
    message("Information Table with GC content count created and saved at ",file.path(Directory), "!")
    message("Done!")
    return(Table_Tot)
    
  }
  if (GCC == FALSE) {
    message("Extracting metadata from .fast5 files!")
    cl <- makeCluster(as.numeric(Cores)) 
    clusterExport(cl, c("HDF5_File_Parsing_Table_Without_GC","PassFiles","Read_Attributes"),envir=environment())
    clusterEvalQ(cl,c(library(rhdf5), library(seqinr)))
    List<-parLapply(cl, c(1:length(PassFiles)), HDF5_File_Parsing_Table_Without_GC,PassFiles)
    stopCluster(cl)
    Table_Tot<-data.frame(matrix(unlist(List), nrow=length(List), ncol=7, byrow=TRUE),stringsAsFactors=FALSE)
    colnames(Table_Tot)<-c("Read Id", "Channel Number", "Mux Number", "Unix Time", "Length of Read", "Quality", "GC_Content")
    write.table(Table_Tot, file.path(Directory, paste0(label, "_Information_Table.txt")), col.names=T, row.names=F, quote=F, sep="\t")
    message("Information Table without GC content count created and saved at ",file.path(Directory), "!")
    message("Done!")
    return(Table_Tot)
    
  }
}

