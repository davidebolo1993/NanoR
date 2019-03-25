
############################################ NanoTableG ############################################


#' @title Generates metadata table
#' @description NanoTableG filters the sequencing summary table retaining only the most-useful statistics and optionall extract GC content from .fastq files
#' @param NanoGList Object of class list returned by NanoPrepareG
#' @param DataOut Where the metadata table will be saved. Do not use a directory that already contains other results.
#' @param GCC  Logical. If TRUE, NanoTableG computes GC content for each sequence in .fastq files. Default to FALSE. Calculating GCC can be slow in R.
#' @return Metadata table with 7 columns
#' @examples
#' #do not run
#' DataOut <- "/path/to/DataOut"
#' # Need a list previously generated with NanoPrepareG()
#' # If List from NanoPrepareM() exists:
#' # Skip GC content calculation
#' Table<-NanoTableG(List, DataOut)
#' # Calculate GC content
#' Table<-NanoTableG(List, DataOut, GCC=TRUE)


NanoTableG<-function(NanoGList,DataOut,GCC=FALSE) { #using multiple cores to deal with .fastq files is not a good idea
  

  Label<-NanoGList[[3]]    
  #Flowcell_ID_Label<-as.character(NanoGList[[2]][1,1])

  #FlowCellId

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
  #else {
    #Flowcell_ID_Label<-as.character("FC5")
  #}


  Directory<-file.path(DataOut,Label)
  dir.create(Directory,showWarnings=FALSE, recursive=TRUE)

  TableInDirectory<-list.files(Directory,pattern="metadata.txt")

  if(length(TableInDirectory) != 0) {

    stop("Cannot use a directory that already contains other results")
  
  }

  Read_Id<-as.character(NanoGList[[2]][,1])
  Channel<-as.numeric(NanoGList[[2]][,2])
  Mux<-NanoGList[[2]][,3]
  Length<-as.numeric(NanoGList[[2]][,5])
  Qscore<-as.numeric(NanoGList[[2]][,6])
  Relative_Time<-as.numeric(NanoGList[[2]][,4])

  if (GCC == TRUE) {
      
    library(ShortRead)
    #library(seqinr)
    message("Calculating GC content...") ## very hard to speed up Fastq parsing in R ... :( 
    # 1 hour for 550 fastq file, 8000 sequences each. Definitely slow.
    
    
    FastqFilesPath<-NanoGList[[1]]
    
    GCC<-function(seq) {      
      GC<-sum(gregexpr('[GgCc]',seq)[[1]] > 0)/nchar(seq) #faster than using library
      return(GC)
    }    
    
    Gc_Con<-function(Element) {
      fqFile<-FastqFile(Element) ## deal with possible error. Not possible to remove single ill-formatted .fastq as far as I know.
      Fastq <- tryCatch({
        readFastq(fqFile)},
        error = function(cond) {
          return(NULL)},
        warning = function(cond) {
          message(cond)
          return(NULL)}
      )
      if (is.null(Fastq)) {
        warning("Ill-formatted .fastq file: ", Element, " .Skipped")
      }
      else {
        CharRead<-as.character(sread(Fastq))
        close(fqFile)
        GC<-lapply(CharRead,GCC)
      }
      return(GC)
    }

    List<-lapply(FastqFilesPath,Gc_Con)      
    GC_Content<-unlist(List)
    GCL<-length(GC_Content)
    NT<-nrow(NanoGList[[2]])
    
    if (GCL != NT) { ##lack of some fastq sequences . Get GC content from what we have in the future. Do not sample what we have
      LackOfReads<-NT-GCL
      GC_To_Add<-rep(NA, LackOfReads)
      GC_Content<-c(GC_Content,GC_To_Add)
    }
    
    Table_Tot<-cbind(Read_Id,Channel,Mux,Relative_Time,Length,Qscore,GC_Content)
    colnames(Table_Tot)<-c("Read Id","Channel Number","Mux Number","Relative Time","Length of Read","Quality", "GC content")
    write.table(Table_Tot, file.path(Directory, 'metadata.txt'), sep="\t", quote=FALSE, col.names=T, row.names=FALSE) 
  
  }
  
  else {
        
    GC_Content<-rep("GC_Content",nrow(NanoGList[[2]]))
    Table_Tot<-cbind(Read_Id,Channel,Mux,Relative_Time,Length,Qscore,GC_Content)
    colnames(Table_Tot)<-c("Read Id","Channel Number","Mux Number","Relative Time","Length of Read", "Quality", "GC content")
    write.table(Table_Tot, file.path(Directory, 'metadata.txt'), sep="\t", quote=FALSE, col.names=T, row.names=FALSE) 
  }

  message("Done")
  return(Table_Tot)
  
}

