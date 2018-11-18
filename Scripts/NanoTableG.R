
############################################ NanoTableG ############################################


#' @title Generates an information table for your GridION X5 .fast5 files
#' @description NanoTableG generates a table that contains useful informations for each read identified by NanoPrepareG function. When analyzing GridION X5 basecalled .fast5 files, metadata extraction can be accelerated using multiple cores. Please note that sometimes .fastq files returned by GridION X5 can be ill-formatted: in that case, NanoTableG will stop and you can run this function again after set GCC to FALSE. This problem is avoided if basecalled .fast5 files are used.
#' @param NanoPrepareGList The object of class "list" returned by NanoPrepareG function;
#' @param DataOut Where the table will be saved. Don't use a directory that already contains NanoTableG (or NanoTableM) result;
#' @param GCC Logical. If TRUE (by default), NanoTableG computes GC content for each read
#' @param Cores Number of cores to be used to accelerate metadata extraction from basecalled .fast5 files. Doesn't affect time when dealing with .fastq and sequencing summary files.
#' @return Table containing informations required by NanoStatsG and NanoCompare (and useful for users who want to analyze these data by theirselves!).
#' @examples
#' #do not run
#' #when working with sequencing summary files and .fastq files
#' NanoGTable<-NanoTableG(NanoPrepareGList=NanoGList, DataOut="/Path/To/DataOutEx", GCC=TRUE) #set GCC to "FALSE" on error
#' #when working with basecalled .fast5 files
#' NanoGTable<-NanoTableG(NanoPrepareGList=NanoGList, DataOut="/Path/To/DataOutEx", GCC=TRUE, Cores=6) #accelerate using multiple cores 




NanoTableG<-function(NanoPrepareGList,DataOut,Cores=1,GCC=TRUE) {
  
  Directory<-file.path(DataOut)
  dir.create(Directory,showWarnings=FALSE)
  Label<-NanoPrepareGList[[6]]
  
  if (is.na(NanoPrepareGList[[1]][1]) == FALSE & is.na(NanoPrepareGList[[3]][1])) {
    
    
    Table<-NanoPrepareGList[[2]]
    
    Flowcell_ID_Label<-unique(as.character(Table[,1]))
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
    Read_Id<-as.character(Table[,2])
    Channel<-as.numeric(Table[,3])
    Mux<-rep("Mux", nrow(Table))
    Length<-as.numeric(Table[,5])
    Qscore<-as.numeric(Table[,6])
    Relative_Time<-as.numeric(Table[,4])
    
    if (GCC == TRUE) {
      
      
      library(ShortRead)
      library(seqinr)
      message("Extracting .fastq sequences and calculating GC content...")
      
      
      FastqFilesPath<-NanoPrepareGList[[1]]
      
      Content<-function(CharRead) {
        Gc<-GC(s2c(CharRead))
        return(Gc)
      }
      
      
      Gc_Con<-function(Element) {
        fqFile<-FastqFile(Element)
        Fastq <- tryCatch({
          readFastq(fqFile)},
          error = function(cond) {
            return(NULL)},
          warning = function(cond) {
            message(cond)
            return(NULL)}
        )
        if (is.null(Fastq)) {
          stop(paste0("Ill-formatted .fastq file: fastq", i-1, ". Analysis aborted. Retry NanoTableG setting GCC to FALSE!"))
        }
        else {
          CharRead<-as.character(sread(Fastq))
          close(fqFile)
          Gc<-lapply(CharRead,Content)
        }
        return(Gc)
      }
      
      List<-lapply(FastqFilesPath,Gc_Con)
      
      
      GC_Content<-unlist(List)
      
      if (length(GC_Content) != dim(Table)[1]) { ##lack of some fastq sequences
        LackOfReads<-dim(Table)[1]-length(GC_Content)
        message(LackOfReads," missing GC content counts: creating ",LackOfReads," fake GC content count...")
        GC_To_Add<-sample(GC_Content,LackOfReads)
        GC_Content<-c(GC_Content,GC_To_Add)
      }
      
      Table_Tot<-cbind(Read_Id,Channel,Mux,Relative_Time,Length,Qscore,GC_Content)
      colnames(Table_Tot)<-c("Read Id","Channel Number","Mux Number","Relative Time","Length of Read","Quality", "GC content")
      write.table(Table_Tot, file.path(DataOut, paste0(Label,"_",Flowcell_ID_Label,"_Information_Table.txt")), sep="\t", quote=FALSE,col.names=T, row.names=FALSE) 
      message("Information Table with GC content count created and saved at ",file.path(Directory), "!")
      message("Done!")
    }
    
    else {
      
      
      GC_Content<-rep("GC_Content",dim(Table)[1])
      Table_Tot<-cbind(Read_Id,Channel,Mux,Relative_Time,Length,Qscore,GC_Content)
      colnames(Table_Tot)<-c("Read Id","Channel Number","Mux Number","Relative Time","Length of Read","Quality", "GC content")
      write.table(Table_Tot, file.path(DataOut, paste0(Label,"_",Flowcell_ID_Label,"_Information_Table.txt")), sep="\t", quote=FALSE,col.names=T, row.names=FALSE) 
      message("Information Table without GC content count created and saved at ",file.path(Directory), "!")
      message("Done!")
    }
    
    return(Table_Tot)
  }
  
  else {
    
    Fast5Data<-NanoPrepareGList[[3]]
    
    library(rhdf5)
    library(seqinr)
    
    
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
      Table<-c(FlowCell_Id="FlowCell_Id",Read_Id="Read_Id",Channel="Channel",Mux="Mux",Unix_Time="Unix_Time",Length="Length",Qscore="Qscore",GC_Content="GC_Content")
      
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
        FlowcellId<-Read_Attributes(Time2, "device_id")
        Table['FlowCell_Id']<-FlowcellId
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
      Table<-c(FlowCell_Id="FlowCell_Id",Read_Id="Read_Id",Channel="Channel",Mux="Mux",Unix_Time="Unix_Time",Length="Length",Qscore="Qscore",GC_Content="GC_Content")
      
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
        FlowcellId<-Read_Attributes(Time2, "device_id")
        Table['FlowCell_Id']<-FlowcellId
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
    
    if (GCC == TRUE) {
      message("Extracting .fast5 files metadata and calculating GC content!")
      cl <- makeCluster(as.numeric(Cores)) 
      clusterExport(cl, c("HDF5_File_Parsing_Table_With_GC","Fast5Data","Read_DataSet","Read_Attributes"),envir=environment())
      clusterEvalQ(cl,c(library(rhdf5), library(seqinr)))
      List<-parLapply(cl, c(1:length(Fast5Data)), HDF5_File_Parsing_Table_With_GC,Fast5Data)
      stopCluster(cl)
      Table_Tot<-data.frame(matrix(unlist(List), nrow=length(List), ncol=8, byrow=TRUE),stringsAsFactors=FALSE)
      Flowcell_ID_Label<-unique(as.character(Table_Tot[,1]))
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
      Table_Tot<-Table_Tot[,2:8]
      colnames(Table_Tot)<-c("Read Id", "Channel Number", "Mux Number", "Unix Time", "Length of Read", "Quality", "GC_Content")
      write.table(Table_Tot, file.path(Directory, paste0(Label,"_",Flowcell_ID_Label, "_Information_Table.txt")), col.names=T, row.names=F, quote=F, sep="\t")
      message("Information Table with GC content count created and saved at ",file.path(Directory), "!")
      message("Done!")
      return(Table_Tot)
      
    }
    if (GCC == FALSE) {
      message("Extracting .fast5 files metadata!")
      cl <- makeCluster(as.numeric(Cores)) 
      clusterExport(cl, c("HDF5_File_Parsing_Table_Without_GC","Fast5Data","Read_Attributes"),envir=environment())
      clusterEvalQ(cl,c(library(rhdf5), library(seqinr)))
      List<-parLapply(cl, c(1:length(Fast5Data)), HDF5_File_Parsing_Table_Without_GC,Fast5Data)
      stopCluster(cl)
      Table_Tot<-data.frame(matrix(unlist(List), nrow=length(List), ncol=8, byrow=TRUE),stringsAsFactors=FALSE)
      Flowcell_ID_Label<-unique(as.character(Table_Tot[,1]))
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
      Table_Tot<-Table_Tot[,2:8]
      colnames(Table_Tot)<-c("Read Id", "Channel Number", "Mux Number", "Unix Time", "Length of Read", "Quality", "GC_Content")
      write.table(Table_Tot, file.path(Directory, paste0(Label,"_",Flowcell_ID_Label, "_Information_Table.txt")), col.names=T, row.names=F, quote=F, sep="\t")
      message("Information Table without GC content count created and saved at ",file.path(Directory), "!")
      message("Done!")
      return(Table_Tot)
      
    }
  }
}

