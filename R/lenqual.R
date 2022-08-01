#' @title lenqual
#' @description generate an interactive HTML report with reads length and quality over time
#' @param summary path to a sequencing summary file from a single ONT run
#' @param time hours to split the sequencing run into. Default to 1
#' @param out path to HTML file
#' @return HTML file
#' @examples
#' #do not run
#' summary<-file.path("/path/to/sequencing_summary.txt")
#' outhtml<-file.path("/path/to/lenqual.html")
#' NanoR::lenqual(summary=summary,time=1,out=outhtml)


lenqual<-function(summary,time=1,out) {

	summary<-normalizePath(file.path(summary))
	out<-file.path(out)

	message("[",Sys.time(),"]"," reading sequencing summary file")

	tab<-fread(summary,sep="\t",header=TRUE,showProgress=FALSE)

	run<-unique(tab$run_id)

	if (length(run) != 1) {

		stop("[",Sys.time(),"]"," sequencing summary must contain only one run id")

	}

	message("[",Sys.time(),"]"," calculating stats")

	unixtime<-tab$template_start+tab$template_duration
	runtime<-round(as.numeric(difftime(as.POSIXct(max(unixtime),origin="1/1/1970"), as.POSIXct(min(unixtime), origin="1/1/1970"), units="hours")))
	rescaled<-rescale(unixtime, to=c(0,runtime))
	tab$template_unix<-rescaled
	bins<-seq(from=0,to=max(rescaled),by=time)

	len<-qual<-data.frame(matrix(0, ncol = 7, nrow = length(bins)-1),stringsAsFactors=FALSE)
	colnames(len)<-colnames(qual)<-c("x", "y", "y1", "y2", "y3", "y4", "y5")

	for (i in c(1:(length(bins)-1))) {

		from<-bins[i]
		to<-bins[i+1]
		subtab<-subset(tab, (tab$template_unix > from & tab$template_unix <= to))

		len[i,]$x<-qual[i,]$x<-paste(from, to, sep="-")

		len[i,]$y<-mean(subtab$sequence_length_template)
		qual[i,]$y<-mean(subtab$mean_qscore_template)

		len[i,]$y1<-min(subtab$sequence_length_template)
		qual[i,]$y1<-min(subtab$mean_qscore_template)

		len[i,]$y2<-max(subtab$sequence_length_template)
		qual[i,]$y2<-max(subtab$mean_qscore_template)

		len[i,]$y3<-mean(subset(subtab, (subtab$passes_filtering == TRUE))$sequence_length_template)
		qual[i,]$y3<-mean(subset(subtab, (subtab$passes_filtering == TRUE))$mean_qscore_template)

		len[i,]$y4<-min(subset(subtab, (subtab$passes_filtering == TRUE))$sequence_length_template)
		qual[i,]$y4<-min(subset(subtab, (subtab$passes_filtering == TRUE))$mean_qscore_template)

		len[i,]$y5<-max(subset(subtab, (subtab$passes_filtering == TRUE))$sequence_length_template)
		qual[i,]$y5<-max(subset(subtab, (subtab$passes_filtering == TRUE))$mean_qscore_template)

	}

	len$x<-qual$x<-factor(len$x,levels=as.character(len$x))


	f <- list(
		size = 10,
		color = "#7f7f7f"
	)



	chart_type <- list(
		type = "dropdown",
		xanchor = 'center',
		yanchor = "bottom",
		buttons = list(
		list(label = "all",
			method = "restyle",
			args = list("visible", list(TRUE, TRUE,TRUE,FALSE,FALSE,FALSE))),
		list(label = "pass",
			method = "restyle",
			args = list("visible", list(FALSE,FALSE,FALSE,TRUE,TRUE,TRUE)))
		)
	)

	message("[",Sys.time(),"]"," plotting")

	p1 <- plot_ly(len, x = ~x, y = ~y1, type = 'scatter', mode = 'lines', fill = 'tozeroy',fillcolor='rgba(0,100,80,0.1)',name = 'Min')
	p1 <- p1 %>% add_trace(y = ~y2, type = 'scatter', mode = 'lines', fill = 'tozeroy',fillcolor='rgba(0,100,80,0.1)', name = 'Max')
	p1 <- p1 %>% add_trace(y = ~y, type = 'scatter', mode = 'lines', fill = 'tozeroy', fillcolor='rgba(0,100,80,0.1)', name = 'Average')
	p1 <- p1 %>% add_trace(y = ~y4, type = 'scatter', mode = 'lines', name = 'Min', fill = 'tozeroy',fillcolor='rgba(0,100,80,0.1)', visible=FALSE)
	p1 <- p1 %>% add_trace(y = ~y5, type = 'scatter', mode = 'lines', fill = 'tozeroy', name = 'Max', fillcolor='rgba(0,100,80,0.1)', visible=FALSE)
	p1 <- p1 %>% add_trace(y = ~y3, type = 'scatter', mode = 'lines', name = 'Average',fill = 'tozeroy',fillcolor='rgba(0,100,80,0.1)', visible=FALSE)


	p1 <- p1 %>% layout(xaxis = list(title = "Sequencing run-time (h)",
	                    showgrid = TRUE,
	                    showticklabels = TRUE,
	                    titlefont=f),
	                    yaxis = list(title = "Length",
	                    showgrid = TRUE,
	                    showline = FALSE,
	                    showticklabels = TRUE,
	                    type = "log",titlefont=f),updatemenus = list(chart_type))


	p2 <- plot_ly(qual, x = ~x, y = ~y1, type = 'scatter', mode = 'lines', fill = 'tozeroy',fillcolor='rgba(0,100,80,0.1)', name = 'Min')
	p2 <- p2 %>% add_trace(y = ~y2, type = 'scatter', mode = 'lines', fill = 'tozeroy', fillcolor='rgba(0,100,80,0.1)', name = 'Max')
	p2 <- p2 %>% add_trace(y = ~y, type = 'scatter', mode = 'lines', fill = 'tozeroy', fillcolor='rgba(0,100,80,0.1)',name = 'Average')
	p2 <- p2 %>% add_trace(y = ~y4, type = 'scatter', mode = 'lines', name = 'Min', fill = 'tozeroy', fillcolor='rgba(0,100,80,0.1)', visible=FALSE)
	p2 <- p2 %>% add_trace(y = ~y5, type = 'scatter', mode = 'lines', fill = 'tozeroy', fillcolor='rgba(0,100,80,0.1)', name = 'Max', visible=FALSE)
	p2 <- p2 %>% add_trace(y = ~y3, type = 'scatter', mode = 'lines', name = 'Average',fill = 'tozeroy', fillcolor='rgba(0,100,80,0.1)',visible=FALSE)


	p2 <- p2 %>% layout(xaxis = list(title = "Sequencing run-time (h)",
	                    showgrid = TRUE,
	                    showticklabels = TRUE,
	                    titlefont=f),
	                    yaxis = list(title = "Qscore",
	                    showgrid = TRUE,
	                    showline = FALSE,
	                    showticklabels = TRUE,
	                    titlefont=f),updatemenus = list(chart_type), title="Reads length and qscore over time")


	s <- subplot(p1,p2,nrows = 2,shareX = TRUE,titleX=TRUE,titleY=TRUE)
	fig <- layout(s, hovermode = "x unified")
	
	message("[",Sys.time(),"]"," storing plot to file")

	htmlwidgets::saveWidget(fig, out)
	
	rm(list=ls())
	invisible(gc())

	message("[",Sys.time(),"]"," done")



}
