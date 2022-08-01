#' @title lenqualdist
#' @description generate an interactive HTML report with length and quality distributions for a subset of reads
#' @param summary path to a sequencing summary file from a single ONT run
#' @param fraction fraction of reads to sample from the original sequencing summary file. Default to 0.01
#' @param out path to HTML file
#' @return HTML file
#' @examples
#' #do not run
#' summary<-file.path("/path/to/sequencing_summary.txt")
#' outhtml<-file.path("/path/to/lengthqualitydist.html")
#' NanoR::lenqualdist(summary=summary,fraction=.01,out=outhtml)


lenqualdist<-function(summary,fraction=.01,out) {

	summary<-normalizePath(file.path(summary))
	out<-file.path(out)

	message("[",Sys.time(),"]"," reading sequencing summary file")

	tab<-fread(summary,sep="\t",header=TRUE,showProgress=FALSE)

	run<-unique(tab$run_id)

	if (length(run) != 1) {

		stop("[",Sys.time(),"]"," sequencing summary must contain only one run id")

	}

	message("[",Sys.time(),"]"," calculating stats")

	chart_type <- list(
		type = "dropdown",
		xanchor = 'center',
		yanchor = "bottom",
		buttons = list(
		list(label = "all",
			method = "restyle",
			args = list("visible", list(TRUE,FALSE))),
		list(label = "pass",
			method = "restyle",
			args = list("visible", list(FALSE,TRUE)))
		)
	)

	f <- list(
		size = 10,
		color = "#7f7f7f"
	)

	tab_s<-tab %>% sample_frac(fraction)

	message("[",Sys.time(),"]"," plotting")

	p1<-plot_ly(x = tab_s$sequence_length_template, type = "histogram",color = I("darkred"), name="length histogram", visible=TRUE) %>% 
		add_trace(x= subset(tab_s, (tab_s$passes_filtering==TRUE))$sequence_length_template,type = "histogram", color = I("darkred"),name="length histogram", visible=FALSE)%>%
		layout(xaxis=list(titlefont=f,title="Length"), yaxis=list(titlefont=f, title="Count"), updatemenus = list(chart_type))

	p2<-plot_ly(x = tab_s$sequence_length_template, type = "box", color = I("darkred"), name="length boxplot", visible=TRUE) %>% 
		add_trace(x= subset(tab_s, (tab_s$passes_filtering==TRUE))$sequence_length_template,type = "box", color = I("darkred"),name="length boxplot", visible=FALSE)%>%
		layout(xaxis=list(titlefont=f,title="Length"), yaxis=list(showline=FALSE, showticklabels=FALSE, showgrid=TRUE,ticks="",titlefont=f),updatemenus = list(chart_type))


	p3<-plot_ly(x = tab_s$mean_qscore_template, type = "histogram",color = I("darkblue"), name="qscore histogram", visible=TRUE) %>% 
		add_trace(x= subset(tab_s, (tab_s$passes_filtering==TRUE))$mean_qscore_template,type = "histogram", color = I("darkblue"),name="qscore histogram", visible=FALSE)%>%
		layout(xaxis=list(titlefont=f,title="Qscore"), yaxis=list(titlefont=f), updatemenus = list(chart_type))


	p4<-plot_ly(x = tab_s$mean_qscore_template, type = "box", color = I("darkblue"), name="qscore boxplot", visible=TRUE) %>% 
		add_trace(x= subset(tab_s, (tab_s$passes_filtering==TRUE))$mean_qscore_template,type = "box", color = I("darkblue"),name="qscore boxplot", visible=FALSE)%>%
		layout(xaxis=list(titlefont=f,title="Qscore"), yaxis=list(showline=FALSE, showticklabels=FALSE, showgrid=TRUE,ticks="",titlefont=f),updatemenus = list(chart_type), title="Reads length and quality distribution (subsample)")

	fig<-subplot(p1,p3,p2,p4,nrows=2,titleX=TRUE, titleY=TRUE,shareX=TRUE)

	message("[",Sys.time(),"]"," storing plot to file")

	htmlwidgets::saveWidget(fig, out)

	rm(list=ls())
	invisible(gc())

	message("[",Sys.time(),"]"," done")

}