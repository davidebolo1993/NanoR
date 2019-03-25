############################################ NanoTableM ############################################


#' @title Generates metadata table
#' @description NanoTableM generates a table that contains key informations for each passed .fast5 read. As NanoTableM can take some time for large sets, this function can be accelerated using multiple cores
#' @param NanoMList Object of class list returned by NanoPrepareM
#' @param DataOut Where the metadata table will be saved. Do not use a directory that already contains other results.
#' @param Cores Number of cores to accelerate metadata extraction: 1 by default
#' @param GCC Logical. If TRUE, NanoTableM computes GC content for each read. Default to FALSE
#' @return Metadata table with 7 columns
#' @examples
#' #do not run
#' DataOut <- "/path/to/DataOut"
#' # Need a list previously generated with NanoPrepareM()
#' # If List from NanoPrepareM() exists:
#' # Skip GC content calculation
#' Table<-NanoTableM(List,DataOut,Cores=6)
#' # Calculate GC content
#' Table<-NanoTableM(List,DataOut,Cores=6, GCC=TRUE) 



NanoTableM<-function(NanoMList,DataOut,Cores=1,GCC=FALSE) { #switched to FALSE
  
  library(parallel)
  library(rhdf5)
  #library(seqinr) #no more required
  
  label<-as.character(NanoMList[[4]])
  Directory<-file.path(DataOut, label)
  dir.create(Directory, showWarnings = FALSE, recursive=TRUE)

  TableInDirectory<-list.files(Directory,pattern="metadata.txt")

  if(length(TableInDirectory) != 0) {

    stop("Cannot use a directory that already contains other results")
  
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
      
      #Pre_Fastq<-Read_DataSet(File,Group1)
      Sequence_Fastq<-unlist(strsplit(Try,"\n"))[2]
      #Fastq<-s2c(Sequence_Fastq)
      Table['GC_Content']<-sum(gregexpr('[GgCc]',Sequence_Fastq)[[1]] > 0)/nchar(Sequence_Fastq)
      
      Group2<-"/UniqueGlobalKey/channel_id"
      Chann_File<-H5Gopen(File,Group2)
      Table['Channel']<-Read_Attributes(Chann_File,"channel_number")
      #Sampling <- as.numeric(Read_Attributes(Chann_File,"sampling_rate")) if this will change, uncomment
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

      #get sampling rate autonomously?

      #AlternativeStart<-floor(Start/Sampling)
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


  HDF5_File_Parsing_Table_With_GC_SC<-function(File) {  #when calling with 1 core, use lapply instead of parlapply
    
    h5errorHandling(type="suppress")
    Table<-c(Read_Id="Read_Id",Channel="Channel",Mux="Mux",Unix_Time="Unix_Time",Length="Length",Qscore="Qscore",GC_Content="GC_Content")
    
    File<-H5Fopen(File)
    
    Group1<-"/Analyses/Basecall_1D_000/BaseCalled_template/Fastq"
    Try<-try(Read_DataSet(File,Group1), silent=TRUE) #exclude non-basecalled .fast5 reads (will be counted as failed by NanoStats)
    if (inherits(Try,"try-error")) {
      return(Table)
    }
    
    else {
      
      #Pre_Fastq<-Read_DataSet(File,Group1)
      Sequence_Fastq<-unlist(strsplit(Try,"\n"))[2]
      #Fastq<-s2c(Sequence_Fastq)
      Table['GC_Content']<-sum(gregexpr('[GgCc]',Sequence_Fastq)[[1]] > 0)/nchar(Sequence_Fastq) #use regex instead of package
      
      Group2<-"/UniqueGlobalKey/channel_id"
      Chann_File<-H5Gopen(File,Group2)
      Table['Channel']<-Read_Attributes(Chann_File,"channel_number")
      #Sampling <- as.numeric(Read_Attributes(Chann_File,"sampling_rate")) if this will change, uncomment
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

      #get sampling rate autonomously?

      #AlternativeStart<-floor(Start/Sampling)
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
      #Sampling <- as.numeric(Read_Attributes(Chann_File,"sampling_rate"))
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

      #get sampling rate autonomously?
      
      #AlternativeStart<-floor(Start/Sampling) ##actual sampling rate
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


  HDF5_File_Parsing_Table_Without_GC_SC<-function(File) {  #when calling with 1 core, use lapply instead of parlapply
    
    h5errorHandling(type="suppress")
    Table<-c(Read_Id="Read_Id",Channel="Channel",Mux="Mux",Unix_Time="Unix_Time",Length="Length",Qscore="Qscore",GC_Content="GC_Content")
    
    File<-H5Fopen(File)
    
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
      #Sampling <- as.numeric(Read_Attributes(Chann_File,"sampling_rate"))
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

      #get sampling rate autonomously?
      
      #AlternativeStart<-floor(Start/Sampling) ##actual sampling rate
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

  ###################### NEW FUNCTIONS FOR MULTI-READ .fast5 ##############################


  HDF5_File_Parsing_Table_With_GC_Multiline<-function(i, File) { 
   
    h5errorHandling(type="suppress")  
    File<-H5Fopen(File[i])
    idtab<-h5ls(File, recursive=FALSE, datasetinfo=FALSE)[2]
    Table<-c(Read_Id="Read_Id",Channel="Channel",Mux="Mux",Unix_Time="Unix_Time",Length="Length",Qscore="Qscore",GC_Content="GC_Content")
    List <- rep(list(Table),nrow(idtab))
    
    for (l in 1:nrow(idtab)) {
      
      GroupAnalyses<-try(H5Gopen(File,paste0(idtab[l,],'/Analyses/Basecall_1D_000/BaseCalled_template')), silent=TRUE) #check if fastq exists. It is supposed to exists both for failed and passed
      if (inherits(GroupAnalyses,"try-error")) {
        H5Gclose(GroupAnalyses)
        List[[l]]<-Table
      }

      else {

        Pre_Fastq<-try(Read_DataSet(GroupAnalyses, 'Fastq'), silent=TRUE)

        if (inherits(Pre_Fastq,"try-error")) {
          H5Gclose(GroupAnalyses)     
          List[[l]]<-Table
        }

        else {

          Sequence_Fastq<-unlist(strsplit(Pre_Fastq,"\n"))[2]
          List[[l]]['GC_Content']<-sum(gregexpr('[GgCc]',Sequence_Fastq)[[1]] > 0)/nchar(Sequence_Fastq)
          H5Gclose(GroupAnalyses)
      
          GroupQuality<-H5Gopen(File,paste0(idtab[l,],'/Analyses/Basecall_1D_000/Summary/basecall_1d_template'))
          List[[l]]['Qscore']<-Read_Attributes(GroupQuality, 'mean_qscore')
          List[[l]]['Length']<-Read_Attributes(GroupQuality, 'sequence_length')
          H5Gclose(GroupQuality)
          
          GroupId<-H5Gopen(File,paste0(idtab[l,],'/Raw/'))
          List[[l]]['Read_Id']<-Read_Attributes(GroupId, 'read_id')
          List[[l]]['Mux']<-Read_Attributes(GroupId, 'start_mux')
          Generation<-Read_Attributes(GroupId, 'start_time')
          H5Gclose(GroupId)
          
          GroupChannel<-H5Gopen(File,paste0(idtab[l,],'/channel_id/'))
          List[[l]]['Channel']<-Read_Attributes(GroupChannel, 'channel_number')
          H5Gclose(GroupChannel)
          
          GroupStartTime<-H5Gopen(File,paste0(idtab[l,],'/tracking_id/'))
          
          Time<-Read_Attributes(GroupStartTime, 'exp_start_time')
          AlternativeStart<-floor(Generation/4000) ##actual sampling rate
          
          if (length(unlist(strsplit(Time,"T")))==2) {### avoid problems with UNIX exp start times found on EGA samples
            Y_M_D<-substr(Time,1,10)
            H_M_S<-substr(Time,12,19)
            Time_Vector<-paste(c(Y_M_D,H_M_S), collapse=" ")
            List[[l]]['Unix_Time']<-as.numeric(as.POSIXct(strptime(Time_Vector, "%Y-%m-%d %H:%M:%S")))+AlternativeStart
            
          }
          else {
            List[[l]]['Unix_Time']<-(as.numeric(AlternativeStart)+as.numeric(Time))
          }
          
          H5Gclose(GroupStartTime)        
        }
      }
    }
    
    H5Fclose(File)
    Table<-do.call(rbind,List)
    return(Table)
  }


  HDF5_File_Parsing_Table_With_GC_Multiline_SC<-function(File) { 
     
    h5errorHandling(type="suppress")      
    File<-H5Fopen(File)
    idtab<-h5ls(File, recursive=FALSE, datasetinfo=FALSE)[2]
    Table<-c(Read_Id="Read_Id",Channel="Channel",Mux="Mux",Unix_Time="Unix_Time",Length="Length",Qscore="Qscore",GC_Content="GC_Content")
    List <- rep(list(Table),nrow(idtab))
    
    for (l in 1:nrow(idtab)) {
      
      GroupAnalyses<-try(H5Gopen(File,paste0(idtab[l,],'/Analyses/Basecall_1D_000/BaseCalled_template')), silent=TRUE) #check if fastq exists. It is supposed to exists both for failed and passed
      if (inherits(GroupAnalyses,"try-error")) {
        H5Gclose(GroupAnalyses)
        List[[l]]<-Table
      }

      else {

        Pre_Fastq<-try(Read_DataSet(GroupAnalyses, 'Fastq'), silent=TRUE)

        if (inherits(Pre_Fastq,"try-error")) {
          H5Gclose(GroupAnalyses)     
          List[[l]]<-Table
        }

        else {

          Sequence_Fastq<-unlist(strsplit(Pre_Fastq,"\n"))[2]
          List[[l]]['GC_Content']<-sum(gregexpr('[GgCc]',Sequence_Fastq)[[1]] > 0)/nchar(Sequence_Fastq)
          H5Gclose(GroupAnalyses)
      
          GroupQuality<-H5Gopen(File,paste0(idtab[l,],'/Analyses/Basecall_1D_000/Summary/basecall_1d_template'))
          List[[l]]['Qscore']<-Read_Attributes(GroupQuality, 'mean_qscore')
          List[[l]]['Length']<-Read_Attributes(GroupQuality, 'sequence_length')
          H5Gclose(GroupQuality)
          
          GroupId<-H5Gopen(File,paste0(idtab[l,],'/Raw/'))
          List[[l]]['Read_Id']<-Read_Attributes(GroupId, 'read_id')
          List[[l]]['Mux']<-Read_Attributes(GroupId, 'start_mux')
          Generation<-Read_Attributes(GroupId, 'start_time')
          H5Gclose(GroupId)
          
          GroupChannel<-H5Gopen(File,paste0(idtab[l,],'/channel_id/'))
          List[[l]]['Channel']<-Read_Attributes(GroupChannel, 'channel_number')
          H5Gclose(GroupChannel)
          
          GroupStartTime<-H5Gopen(File,paste0(idtab[l,],'/tracking_id/'))
          
          Time<-Read_Attributes(GroupStartTime, 'exp_start_time')
          AlternativeStart<-floor(Generation/4000) ##actual sampling rate
          
          if (length(unlist(strsplit(Time,"T")))==2) {### avoid problems with UNIX exp start times found on EGA samples
            Y_M_D<-substr(Time,1,10)
            H_M_S<-substr(Time,12,19)
            Time_Vector<-paste(c(Y_M_D,H_M_S), collapse=" ")
            List[[l]]['Unix_Time']<-as.numeric(as.POSIXct(strptime(Time_Vector, "%Y-%m-%d %H:%M:%S")))+AlternativeStart
            
          }
          else {
            List[[l]]['Unix_Time']<-(as.numeric(AlternativeStart)+as.numeric(Time))
          }
          
          H5Gclose(GroupStartTime)        
        }
      }
    }
    
    H5Fclose(File)
    Table<-do.call(rbind,List)
    return(Table)
  }


  HDF5_File_Parsing_Table_Without_GC_Multiline<-function(i,File) {  ## work for R9.4 and R9.5
    
    h5errorHandling(type="suppress")      
    File<-H5Fopen(File[i])
    idtab<-h5ls(File, recursive=FALSE, datasetinfo=FALSE)[2]
    Table<-c(Read_Id="Read_Id",Channel="Channel",Mux="Mux",Unix_Time="Unix_Time",Length="Length",Qscore="Qscore",GC_Content="GC_Content")
    List <- rep(list(Table),nrow(idtab))
    
    for (l in 1:nrow(idtab)) {
      
      GroupAnalyses<-try(H5Gopen(File,paste0(idtab[l,],'/Analyses/Basecall_1D_000/BaseCalled_template')), silent=TRUE) #check if fastq exists. It is supposed to exists both for failed and passed
      if (inherits(GroupAnalyses,"try-error")) {
        H5Gclose(GroupAnalyses)
        List[[l]]<-Table
      }

      else {

        Pre_Fastq<-try(Read_DataSet(GroupAnalyses, 'Fastq'), silent=TRUE)

        if (inherits(Pre_Fastq,"try-error")) {
          H5Gclose(GroupAnalyses)     
          List[[l]]<-Table
        }

        else {

          #Sequence_Fastq<-unlist(strsplit(Pre_Fastq,"\n"))[2]
          List[[l]]['GC_Content']<-'GC_Content'
          H5Gclose(GroupAnalyses)
      
          GroupQuality<-H5Gopen(File,paste0(idtab[l,],'/Analyses/Basecall_1D_000/Summary/basecall_1d_template'))
          List[[l]]['Qscore']<-Read_Attributes(GroupQuality, 'mean_qscore')
          List[[l]]['Length']<-Read_Attributes(GroupQuality, 'sequence_length')
          H5Gclose(GroupQuality)
          
          GroupId<-H5Gopen(File,paste0(idtab[l,],'/Raw/'))
          List[[l]]['Read_Id']<-Read_Attributes(GroupId, 'read_id')
          List[[l]]['Mux']<-Read_Attributes(GroupId, 'start_mux')
          Generation<-Read_Attributes(GroupId, 'start_time')
          H5Gclose(GroupId)
          
          GroupChannel<-H5Gopen(File,paste0(idtab[l,],'/channel_id/'))
          List[[l]]['Channel']<-Read_Attributes(GroupChannel, 'channel_number')
          H5Gclose(GroupChannel)
          
          GroupStartTime<-H5Gopen(File,paste0(idtab[l,],'/tracking_id/'))
          
          Time<-Read_Attributes(GroupStartTime, 'exp_start_time')
          AlternativeStart<-floor(Generation/4000) ##actual sampling rate
          
          if (length(unlist(strsplit(Time,"T")))==2) {### avoid problems with UNIX exp start times found on EGA samples
            Y_M_D<-substr(Time,1,10)
            H_M_S<-substr(Time,12,19)
            Time_Vector<-paste(c(Y_M_D,H_M_S), collapse=" ")
            List[[l]]['Unix_Time']<-as.numeric(as.POSIXct(strptime(Time_Vector, "%Y-%m-%d %H:%M:%S")))+AlternativeStart
            
          }
          else {
            List[[l]]['Unix_Time']<-(as.numeric(AlternativeStart)+as.numeric(Time))
          }
          
          H5Gclose(GroupStartTime)        
        }
      }
    }
    
    H5Fclose(File)
    Table<-do.call(rbind,List)
    return(Table)
  }

  HDF5_File_Parsing_Table_Without_GC_Multiline_SC<-function(File) {  ## work for R9.4 and R9.5
    
    h5errorHandling(type="suppress")      
    File<-H5Fopen(File)
    idtab<-h5ls(File, recursive=FALSE, datasetinfo=FALSE)[2]
    Table<-c(Read_Id="Read_Id",Channel="Channel",Mux="Mux",Unix_Time="Unix_Time",Length="Length",Qscore="Qscore",GC_Content="GC_Content")
    List <- rep(list(Table),nrow(idtab))
    
    for (l in 1:nrow(idtab)) {
      
      GroupAnalyses<-try(H5Gopen(File,paste0(idtab[l,],'/Analyses/Basecall_1D_000/BaseCalled_template')), silent=TRUE) #check if fastq exists. It is supposed to exists both for failed and passed
      if (inherits(GroupAnalyses,"try-error")) {
        H5Gclose(GroupAnalyses)
        List[[l]]<-Table
      }

      else {

        Pre_Fastq<-try(Read_DataSet(GroupAnalyses, 'Fastq'), silent=TRUE)

        if (inherits(Pre_Fastq,"try-error")) {
          H5Gclose(GroupAnalyses)     
          List[[l]]<-Table
        }

        else {

          #Sequence_Fastq<-unlist(strsplit(Pre_Fastq,"\n"))[2]
          List[[l]]['GC_Content']<-'GC_Content'
          H5Gclose(GroupAnalyses)
      
          GroupQuality<-H5Gopen(File,paste0(idtab[l,],'/Analyses/Basecall_1D_000/Summary/basecall_1d_template'))
          List[[l]]['Qscore']<-Read_Attributes(GroupQuality, 'mean_qscore')
          List[[l]]['Length']<-Read_Attributes(GroupQuality, 'sequence_length')
          H5Gclose(GroupQuality)
          
          GroupId<-H5Gopen(File,paste0(idtab[l,],'/Raw/'))
          List[[l]]['Read_Id']<-Read_Attributes(GroupId, 'read_id')
          List[[l]]['Mux']<-Read_Attributes(GroupId, 'start_mux')
          Generation<-Read_Attributes(GroupId, 'start_time')
          H5Gclose(GroupId)
          
          GroupChannel<-H5Gopen(File,paste0(idtab[l,],'/channel_id/'))
          List[[l]]['Channel']<-Read_Attributes(GroupChannel, 'channel_number')
          H5Gclose(GroupChannel)
          
          GroupStartTime<-H5Gopen(File,paste0(idtab[l,],'/tracking_id/'))
          
          Time<-Read_Attributes(GroupStartTime, 'exp_start_time')
          AlternativeStart<-floor(Generation/4000) ##actual sampling rate
          
          if (length(unlist(strsplit(Time,"T")))==2) {### avoid problems with UNIX exp start times found on EGA samples
            Y_M_D<-substr(Time,1,10)
            H_M_S<-substr(Time,12,19)
            Time_Vector<-paste(c(Y_M_D,H_M_S), collapse=" ")
            List[[l]]['Unix_Time']<-as.numeric(as.POSIXct(strptime(Time_Vector, "%Y-%m-%d %H:%M:%S")))+AlternativeStart
            
          }
          else {
            List[[l]]['Unix_Time']<-(as.numeric(AlternativeStart)+as.numeric(Time))
          }
          
          H5Gclose(GroupStartTime)        
        }
      }
    }
    
    H5Fclose(File)
    Table<-do.call(rbind,List)
    return(Table)
  }


  #### PROCESS #####
  
  PassFiles<-NanoMList[[1]] #need to re-assign to another value

  Multi<-NanoMList[[5]]
  
  ### CALCULATIING ###

  if (Multi == FALSE) { #classic behaviour, no multi-line .fast5

  
    if (GCC == TRUE) {

      message("Extracting metadata and calculating GC content from .fast5 files...")

      if (Cores > 1) {

        cl <- makeCluster(as.numeric(Cores)) 
        clusterExport(cl, c("HDF5_File_Parsing_Table_With_GC","PassFiles","Read_DataSet","Read_Attributes"),envir=environment())
        clusterEvalQ(cl,library(rhdf5))
        List<-parLapply(cl, c(1:length(PassFiles)), HDF5_File_Parsing_Table_With_GC,PassFiles)
        stopCluster(cl)
      
      }

      else {

        List<-lapply(PassFiles,HDF5_File_Parsing_Table_With_GC_SC)

      }
    
    }

    else {

      message("Extracting metadata from .fast5 files...")


      if (Cores >1) {

        cl <- makeCluster(as.numeric(Cores)) 
        clusterExport(cl, c("HDF5_File_Parsing_Table_Without_GC","PassFiles","Read_Attributes"),envir=environment())
        clusterEvalQ(cl,library(rhdf5))
        List<-parLapply(cl, c(1:length(PassFiles)), HDF5_File_Parsing_Table_Without_GC,PassFiles)
        stopCluster(cl)

      }

      else {

        List<-lapply(PassFiles,HDF5_File_Parsing_Table_Without_GC_SC)

      }

    }

  }


  else { ## multiline.fast5


    if (GCC == TRUE) {

      message("Extracting metadata and calculating GC content from multi-read .fast5 files...")

      if (Cores > 1) {

        cl <- makeCluster(as.numeric(Cores)) 
        clusterExport(cl, c("HDF5_File_Parsing_Table_With_GC_Multiline","PassFiles","Read_DataSet","Read_Attributes"),envir=environment())
        clusterEvalQ(cl,library(rhdf5))
        List<-parLapply(cl, c(1:length(PassFiles)), HDF5_File_Parsing_Table_With_GC_Multiline,PassFiles)
        stopCluster(cl)
        
      }

      else {

        List<-lapply(PassFiles,HDF5_File_Parsing_Table_With_GC_Multiline_SC)

      }
    
    }

    else {

      message("Extracting metadata from multi-read .fast5 files...")


      if (Cores >1) {

        cl <- makeCluster(as.numeric(Cores)) 
        clusterExport(cl, c("HDF5_File_Parsing_Table_Without_GC_Multiline","PassFiles","Read_Attributes", "Read_DataSet"),envir=environment())
        clusterEvalQ(cl,library(rhdf5))
        List<-parLapply(cl, c(1:length(PassFiles)), HDF5_File_Parsing_Table_Without_GC_Multiline,PassFiles)
        stopCluster(cl)

      }

      else {

        List<-lapply(PassFiles,HDF5_File_Parsing_Table_Without_GC_Multiline_SC)

      }

    }

  }    

  Table_Tot<-do.call(rbind,List)
  colnames(Table_Tot)<-c("Read Id", "Channel Number", "Mux Number", "Unix Time", "Length of Read", "Quality", "GC Content")
  write.table(Table_Tot, file.path(Directory, "metadata.txt"), col.names=T, row.names=F, quote=F, sep="\t")
  message("Done")
  return(Table_Tot)

}
