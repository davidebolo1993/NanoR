#' @title alfredstats
#' @description generate an interactive HTML report with statistics from alfred qc (zgrep "^ME <qc.tsv.gz> | cut -f 2- | datamash transpose)
#' @param qcdata path to alfred qc.tsv.gz 
#' @param out path to HTML report
#' @return HTML file
#' @examples
#' #do not run
#' qcdata<-file.path("/path/to/qc.tsv.gz")
#' outhtml<-file.path("/path/to/qc.html")
#' NanoR::alfredstats(qcdata=qcdata, out=outhtml)

alfredstats<-function(qcdata,out) {

	qcdata<-normalizePath(file.path(qcdata))
	out<-file.path(out)

	message("[",Sys.time(),"]"," reading qc data")

	#read, re-organize
	df<-fread(qcdata)

	#error-rate
	keep<-c("#MatchedBases", "MatchRate",
			"#MismatchedBases", "MismatchRate",
			"#DeletionsCigarD", "DeletionRate",
			"#InsertionsCigarI", "InsertionRate",
			"ErrorRate")

	error_df<- df[which(df$Sample %in% keep),]
	error_df[match(keep, error_df$Sample),]

	error_df<-transpose(error_df)
	names(error_df) <- as.character(unlist(error_df[1,]))
	error_df<-error_df[-1,]
	TotalErrors<-as.integer(error_df[["#MismatchedBases"]]) + as.integer(error_df[["#DeletionsCigarD"]]) + as.integer(error_df[["#InsertionsCigarI"]])
	error_df<-cbind(error_df,TotalErrors)
	colnames(error_df)<-c(colnames(error_df)[c(1:length(colnames(error_df))-1)], c("#TotalErrors"))
	cols<-colnames(error_df)
	error_df<-transpose(error_df)
	error_df$V2<-error_df$V1
	error_df$V1<-cols

	#counts
	error_counts_df<-error_df[grepl("#", error_df$V1),]
	error_counts_df$V1<-factor(error_counts_df$V1, levels=as.character(error_counts_df$V1))
	#rates
	error_rates_df<-error_df[!grepl("#", error_df$V1),]
	error_rates_df$V1<-factor(error_rates_df$V1, levels=as.character(error_rates_df$V1))


	#aligned reads

	keep<-c("#Unmapped", "UnmappedFraction",
			"#Mapped", "MappedFraction",
			"#MappedForward", "MappedForwardFraction",
			"#MappedReverse", "MappedReverseFraction",
			"#SecondaryAlignments", "SecondaryAlignmentFraction",
			"#SupplementaryAlignments", "SupplementaryAlignmentFraction",
			"#SplicedAlignments", "SplicedAlignmentFraction"            
			)

	alignment_df<-df[which(df$Sample %in% keep),]
	colnames(alignment_df)<-c("V1", "V2")

	#counts
	alignment_counts_df<-alignment_df[grepl("#", alignment_df$V1),]
	alignment_counts_df$V1<-factor(alignment_counts_df$V1, levels=as.character(alignment_counts_df$V1))
	#rates
	alignment_rates_df<-alignment_df[!grepl("#", alignment_df$V1),]
	alignment_rates_df$V1<-factor(alignment_rates_df$V1, levels=as.character(alignment_rates_df$V1))

	f <- list(
		size = 10,
		color = "#7f7f7f"
	)


	chart_type <- list(
		type = "dropdown",
		xanchor = 'center',
		yanchor = "bottom",
		buttons = list(
		list(label = "counts",
			method = "restyle",
			args = list("visible", list(TRUE,FALSE))),
		list(label = "rates",
			method = "restyle",
			args = list("visible", list(FALSE,TRUE)))
		)
	)

	message("[",Sys.time(),"]"," plotting")

	p1<-plot_ly() %>% add_trace(x=error_counts_df$V1,y=as.numeric(error_counts_df$V2),type="bar", color = I("darkblue"), name = "#bp", visible=TRUE) %>%
					add_trace(x=error_rates_df$V1,y=as.numeric(error_rates_df$V2),type="bar", color = I("darkred"), name = ":bp", visible=FALSE) %>%
					layout(showlegend=FALSE,updatemenus = list(chart_type))

	p2<-plot_ly() %>% add_trace(x=alignment_counts_df$V1,y=as.numeric(alignment_counts_df$V2),type="bar", color = I("darkblue"), name = "#reads", visible=TRUE) %>%
					add_trace(x=alignment_rates_df$V1,y=as.numeric(alignment_rates_df$V2),type="bar", color = I("darkred"), name = ":reads", visible=FALSE) %>%
					layout(showlegend=FALSE,updatemenus = list(chart_type), title="Alfred statistics")

	
	fig<- subplot(p1,p2,nrows = 2,titleX=TRUE, titleY=TRUE,margin=.05)

	message("[",Sys.time(),"]"," storing plot to file")

	htmlwidgets::saveWidget(fig, out)

	rm(list=ls())
	invisible(gc())

	message("[",Sys.time(),"]"," done")

}
