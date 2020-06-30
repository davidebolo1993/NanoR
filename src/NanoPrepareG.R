############################################ NanoPrepareG ############################################



#' @title NanoPrepareG
#' @description Organize MinION/GridION sequencing summary and FASTQ
#' @param DataSummary Path to sequencing summary
#' @param DataFastq Path to passed FASTQ
#' @details FASTQ are found recursively
#' @return Object of class list
#' @examples
#' #do not run
#' DataSummary<-'path/to/sequencing_summary'
#' DataFastq<-'path/to/fastq_pass'
#' List<-NanoPrepareG(DataSummary,DataFastq)



NanoPrepareG<-function(DataSummary,DataFastq) {

  library(data.table)

  List<-list()
  FastqFiles<-list.files(DataFastq,pattern=".fastq",full.names=TRUE, recursive=TRUE)
  message('Passed FASTQ: ', length(FastqFiles))


  Read_Table_Summary<-function(File) {

    Table<-fread(File,header=TRUE,sep="\t")
    RealativeTimeToAdd<-(as.numeric(Table$template_start)+as.numeric(Table$template_duration))
    Read_Id<-as.character(Table$read_id)
    Channel<-as.numeric(Table$channel)
    Mux<-as.numeric(Table$mux)
    Length<-as.numeric(Table$sequence_length_template)
    Qscore<-as.numeric(Table$mean_qscore_template)
    Table<-cbind(Read_Id,Channel,Mux,RealativeTimeToAdd,Length,Qscore)
    return(Table)
  }

  message('Reading ', file.path(DataSummary), ' ...')

  SummaryTable<-Read_Table_Summary(DataSummary)

  colnames(SummaryTable)<-c("Read Id","Channel Number","Mux Number","Relative Time","Length of Read","Quality")

  List[['fastq']] <-FastqFiles
  List[['summary']]<-SummaryTable

  message("Done")

  return(List)

}


