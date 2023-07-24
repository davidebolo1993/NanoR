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
	error_df<-error_df[match(keep, error_df$Sample),]

	error_df<-transpose(error_df)
	names(error_df) <- as.character(unlist(error_df[1,]))
	error_df<-error_df[-1,]
	TotalErrors<-as.integer(error_df[["#MismatchedBases"]]) + as.integer(error_df[["#DeletionsCigarD"]]) + as.integer(error_df[["#InsertionsCigarI"]])
	error_df<-cbind(error_df,TotalErrors)
	colnames(error_df)<-c(colnames(error_df)[c(1:length(colnames(error_df))-1)], c("#TotalErrors"))
	cols<-colnames(error_df)
	error_df<-transpose(error_df)
  colnames(error_df)<-colnames(df)[c(2:ncol(df))]
  error_df$V0<-cols

	#counts
	error_counts_df<-error_df[grepl("#", error_df$V0),]
	error_counts_df$V0<-factor(error_counts_df$V0, levels=as.character(error_counts_df$V0))
	#rates
	error_rates_df<-error_df[!grepl("#", error_df$V0),]
	error_rates_df$V0<-factor(error_rates_df$V0, levels=as.character(error_rates_df$V0))


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
	colnames(alignment_df)[1]<-c("V0")
	alignment_df<-alignment_df %>% dplyr::relocate(V0, .after = last_col())

	#counts
	alignment_counts_df<-alignment_df[grepl("#", alignment_df$V0),]
	alignment_counts_df$V0<-factor(alignment_counts_df$V0, levels=as.character(alignment_counts_df$V0))
	#rates
	alignment_rates_df<-alignment_df[!grepl("#", alignment_df$V0),]
	alignment_rates_df$V0<-factor(alignment_rates_df$V0, levels=as.character(alignment_rates_df$V0))

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
			args = list("visible", sapply(rep(c(TRUE,FALSE),each=(ncol(error_counts_df)-1)),list))),
		list(label = "rates",
			method = "restyle",
			args = list("visible", sapply(rep(c(FALSE,TRUE),each=(ncol(error_counts_df)-1)),list)))
		)
	)

	message("[",Sys.time(),"]"," plotting")

	p1<-plot_ly()

	for (n in c(1:(ncol(error_counts_df)-1))) {

	  p1<-p1 %>% add_trace(x=error_counts_df$V0, y=as.numeric(unlist(error_counts_df[,..n])), name = colnames(error_counts_df[,..n]), legendgroup= colnames(error_counts_df[,..n]),type="bar", visible=TRUE)

	}

	for (n in c(1:(ncol(error_rates_df)-1))) {

	  p1<-p1 %>% add_trace(x=error_rates_df$V0, y=as.numeric(unlist(error_rates_df[,..n])), name = colnames(error_rates_df[,..n]), legendgroup= colnames(error_rates_df[,..n]), type="bar", visible=FALSE)

	}

	p1<-p1%>%layout(showlegend=FALSE,updatemenus = list(chart_type),  barmode = 'group')


	p2<-plot_ly()

	for (n in c(1:(ncol(alignment_counts_df)-1))) {

	  p2<-p2 %>% add_trace(x=alignment_counts_df$V0, y=as.numeric(unlist(alignment_counts_df[,..n])), name = colnames(alignment_counts_df[,..n]), legendgroup= colnames(alignment_counts_df[,..n]), type="bar", visible=TRUE)

	}

	for (n in c(1:(ncol(alignment_rates_df)-1))) {

	  p2<-p2 %>% add_trace(x=alignment_rates_df$V0, y=as.numeric(unlist(alignment_rates_df[,..n])), name = colnames(alignment_rates_df[,..n]), legendgroup= colnames(alignment_rates_df[,..n]), type="bar", visible=FALSE)

	}

	p2<-p2%>%layout(showlegend=TRUE,updatemenus = list(chart_type),  barmode = 'group')

	fig<- subplot(p1,p2,nrows = 2,titleX=TRUE, titleY=TRUE,margin=.05)

	message("[",Sys.time(),"]"," storing plot to file")

	htmlwidgets::saveWidget(fig, out)

	rm(list=ls())
	invisible(gc())

	message("[",Sys.time(),"]"," done")

}
