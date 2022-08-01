#' @title report
#' @description generate a TSV file with several informations summarizing the sequencing run 
#' @param summary path to a sequencing summary file from a single ONT run
#' @param out path to TSV file
#' @return TSV file 
#' @examples
#' #do not run
#' summary<-file.path("/path/to/sequencing_summary.txt")
#' outtsv<-file.path("/path/to/summary.tsv")
#' NanoR::report(summary=summary,out=outtsv)


report<-function(summary,out) {

	summary<-normalizePath(file.path(summary))
	out<-file.path(out)

	message("[",Sys.time(),"]"," reading sequencing summary file")

	tab<-fread(summary,sep="\t",header=TRUE,showProgress=FALSE)

	run<-unique(tab$run_id)

	if (length(run) != 1) {

		stop("[",Sys.time(),"]"," sequencing summary must contain only one run id")

	}

	message("[",Sys.time(),"]"," generating report")

	read_type="all"
	min_length<-min(tab$sequence_length_template)
	max_length<-max(tab$sequence_length_template)
	mean_length<-mean(tab$sequence_length_template)
	len.sorted <- rev(sort(tab$sequence_length_template))
	n50_length<-len.sorted[cumsum(as.numeric(len.sorted)) >= sum(len.sorted)*0.5][1]
	q_length<-as.numeric(quantile(tab$sequence_length_template))
	q1_length<-q_length[2]
	median_length<-q_length[3]
	q3_length<-q_length[4]
	longest_read<-tab$read_id[which.max(tab$sequence_length_template)]
	min_qscore<-min(tab$mean_qscore_template)
	max_qscore<-max(tab$mean_qscore_template)
	mean_qscore<-mean(tab$mean_qscore_template)
	q_qual<-as.numeric(quantile(tab$mean_qscore_template))
	highest_quality_read<-tab$read_id[which.max(tab$mean_qscore_template)]
	q1_qscore<-q_qual[2]
	median_qscore<-q_qual[3]
	q3_qscore<-q_qual[4]
	reads<-nrow(tab)
	bases<-sum(tab$sequence_length_template)
	pass_ratio<-nrow(subset(tab, (tab$passes_filtering == TRUE)))/nrow(tab)
	fail_ratio<-1-pass_ratio

	values_all<-data.table(cbind(run,read_type,min_length,mean_length,max_length,q1_length,median_length,q3_length,n50_length,longest_read,min_qscore,mean_qscore,max_qscore,q1_qscore,median_qscore,q3_qscore,highest_quality_read,reads,bases,pass_ratio,fail_ratio))

	#pass
	read_type="pass"
	tab_s<-subset(tab, (tab$passes_filtering == TRUE))
	min_length<-min(tab_s$sequence_length_template)
	max_length<-max(tab_s$sequence_length_template)
	mean_length<-mean(tab_s$sequence_length_template)
	len.sorted <- rev(sort(tab_s$sequence_length_template))
	n50_length<-len.sorted[cumsum(as.numeric(len.sorted)) >= sum(len.sorted)*0.5][1]
	q_length<-as.numeric(quantile(tab_s$sequence_length_template))
	q1_length<-q_length[2]
	median_length<-q_length[3]
	q3_length<-q_length[4]
	longest_read<-tab_s$read_id[which.max(tab_s$sequence_length_template)]
	min_qscore<-min(tab_s$mean_qscore_template)
	max_qscore<-max(tab_s$mean_qscore_template)
	mean_qscore<-mean(tab_s$mean_qscore_template)
	q_qual<-as.numeric(quantile(tab_s$mean_qscore_template))
	highest_quality_read<-tab_s$read_id[which.max(tab_s$mean_qscore_template)]
	q1_qscore<-q_qual[2]
	median_qscore<-q_qual[3]
	q3_qscore<-q_qual[4]
	reads<-nrow(tab_s)
	bases<-sum(tab_s$sequence_length_template)
	pass_ratio<-1
	fail_ratio<-0

	values_pass<-data.table(cbind(run,read_type,min_length,mean_length,max_length,q1_length,median_length,q3_length,n50_length,longest_read,min_qscore,mean_qscore,max_qscore,q1_qscore,median_qscore,q3_qscore,highest_quality_read,reads,bases,pass_ratio,fail_ratio))
	values<-data.table(rbind(values_all,values_pass))

	message("[",Sys.time(),"]"," writing report to file")

	fwrite(values,file = out, sep = "\t", quote=FALSE,row.names = FALSE,col.names = TRUE)
	rm(list=ls())
	invisible(gc())

	message("[",Sys.time(),"]"," done")


}