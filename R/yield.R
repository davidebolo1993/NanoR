#' @title yield
#' @description generate an interactive HTML report with number of reads and bases generated during the sequencing run
#' @param summary path to a sequencing summary file from a single ONT run
#' @param time hours to split the sequencing run into. Default to 1
#' @param out path to HTML file
#' @return HTML file
#' @examples
#' #do not run
#' summary<-file.path("/path/to/sequencing_summary.txt")
#' outhtml<-file.path("/path/to/yield.html")
#' NanoR::yield(summary=summary,time=1,out=outhtml)


yield<-function(summary,time=1,out) {

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

	yield<-yield_pass<-data.frame(matrix(0, ncol = 5, nrow = length(bins)-1),stringsAsFactors=FALSE)
	colnames(yield)<-colnames(yield_pass)<-c("x", "y", "y2", "y3", "y4")

	for (i in c(1:(length(bins)-1))) {

		from<-bins[i]
		to<-bins[i+1]
		subtab<-subset(tab, (tab$template_unix > from & tab$template_unix <= to))

		yield[i,]$x<-yield_pass[i,]$x<-paste(from, to, sep="-")
		yield[i,]$y<-nrow(subtab)
		yield_pass[i,]$y<-nrow(subset(subtab, (subtab$passes_filtering == TRUE)))
		yield[i,]$y2<-sum(subtab$sequence_length_template)
		yield_pass[i,]$y2<-sum(subset(subtab, (subtab$passes_filtering == TRUE))$sequence_length_template)

	}


	yield$y3<-cumsum(yield$y)
	yield_pass$y3<-cumsum(yield_pass$y)
	yield$y4<-cumsum(yield$y2)
	yield_pass$y4<-cumsum(yield_pass$y2)


	#plot #read and #bp, all and passed reads

	yield$x<-yield_pass$x<-factor(yield$x,levels=as.character(yield$x))


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
			args = list("visible", list(TRUE, TRUE,FALSE,FALSE))),
		list(label = "pass",
			method = "restyle",
			args = list("visible", list(FALSE,FALSE,TRUE,TRUE)))
		)
	)

	message("[",Sys.time(),"]"," plotting")

	p1 <- plot_ly() %>% add_trace(x=yield$x,y=yield$y,type="scatter",mode="lines+markers",yaxis="y",name="#reads", visible=TRUE) %>%
						add_trace(x=yield$x,y=yield$y2,type="scatter",mode="lines+markers",yaxis="y2",name="#bp",visible=TRUE) %>%
						add_trace(x=yield_pass$x,y=yield_pass$y,type="scatter",mode="lines+markers",yaxis="y",name="#reads", visible=FALSE) %>%
						add_trace(x=yield_pass$x,y=yield_pass$y2,type="scatter",mode="lines+markers",yaxis="y2",name="#bp",visible=FALSE) %>%
						layout(yaxis=list(side="left",title="#reads",titlefont=f),
							yaxis2=list(side="right",overlaying="y",title="#bases",showgrid=FALSE,titlefont=f),
							xaxis=list(title="Sequencing run-time (h)",titlefont=f,tickangle=45),
							showlegend=TRUE,updatemenus = list(chart_type))

	p2 <- plot_ly() %>% add_trace(x=yield$x,y=yield$y3,type="scatter",mode="lines+markers",yaxis="y",name="#reads", visible=TRUE) %>%
						add_trace(x=yield$x,y=yield$y4,type="scatter",mode="lines+markers",yaxis="y2",name="#bp",visible=TRUE) %>%
						add_trace(x=yield_pass$x,y=yield_pass$y3,type="scatter",mode="lines+markers",yaxis="y",name="#reads", visible=FALSE) %>%
						add_trace(x=yield_pass$x,y=yield_pass$y4,type="scatter",mode="lines+markers",yaxis="y2",name="#bp",visible=FALSE) %>%
						layout(yaxis=list(side="left",title="#reads",titlefont=f),
							yaxis2=list(side="right",overlaying="y3",title="#bases",showgrid=FALSE,titlefont=f),
							xaxis=list(title="Sequencing run-time (h)",titlefont=f,tickangle=45),
							showlegend=TRUE,updatemenus = list(chart_type),title = '#reads and #bp (over time and cumulative)')

	s <- subplot(p1,p2,nrows = 2,shareX = TRUE,titleX=TRUE,titleY=TRUE)
	fig <- layout(s, hovermode = "x unified")

	message("[",Sys.time(),"]"," storing plot to file")

	htmlwidgets::saveWidget(fig, out)

	rm(list=ls())
	invisible(gc())

	message("[",Sys.time(),"]"," done")


}