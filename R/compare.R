#' @title compare
#' @description generate an interactive HTML report comparing length and quality of reads from multiple sequencing summary files
#' @param summaries a character vector containing paths to multiple sequencing summaries (each from a different ONT sequencing run)
#' @param time compare run statistics using this time interval (hours). Default to 10
#' @param out path to HTML report
#' @return HTML file
#' @examples
#' #do not run
#' exp1<-file.path("/path/to/exp1_sequencing_summary.txt")
#' exp2<-file.path("/path/to/exp2_sequencing_summary.txt")
#' outhtml<-file.path("/path/to/comparison.html")
#' NanoR::compare(summaries=c(exp1,exp2), fraction=.001, out=outhtml)


compare<-function(summaries,time=10,out){


	out<-file.path(out)
	bins<-seq(from=0,to=100,by=1)
	taball<-taball_pass<-list()
	compare_<-seq(from=0,to=100,by=time)
	clist<-list()

	for (i in c(1:(length(compare_)-1))) {

		h<-seq(from=compare_[i], to=compare_[i+1], by=1)
		
		for (n in h) {

			clist[[as.character(n)]]<-paste(h[1], h[length(h)],sep="-")
		
		}
	}


	for (summary in summaries) {

		summary<-normalizePath(file.path(summary))

		message("[",Sys.time(),"]"," reading sequencing summary file ", summary)

		tab<-fread(summary,sep="\t",header=TRUE,showProgress=FALSE)
		run<-unique(tab$run_id)

		if (length(run) != 1) {

			stop("[",Sys.time(),"]"," sequencing summary must contain only one run id")

		}

		summ<-summ_pass<-data.frame(matrix(0, ncol = 5, nrow = length(bins)-1),stringsAsFactors=FALSE)
		colnames(summ)<-colnames(summ_pass)<-c("x", "y", "y2", "y3", "y4")

		unixtime<-tab$template_start+tab$template_duration
		runtime<-round(as.numeric(difftime(as.POSIXct(max(unixtime),origin="1/1/1970"), as.POSIXct(min(unixtime), origin="1/1/1970"), units="hours")))
		rescaled<-rescale(unixtime, to=c(0,runtime))
		tab$template_unix<-rescaled

		message("[",Sys.time(),"]"," calculating stats")

		for (i in c(1:(length(bins)-1))) {

			from<-bins[i]
			to<-bins[i+1]
			subtab<-subset(tab, (tab$template_unix > from & tab$template_unix <= to))

			range_<-clist[[as.character(from)]]

			summ[i,]$x<-summ_pass[i,]$x<-range_
			
			summ[i,]$y<-nrow(subtab)
			summ_pass[i,]$y<-nrow(subset(subtab, (subtab$passes_filtering == TRUE)))

			summ[i,]$y2<-sum(subtab$sequence_length_template)
			summ_pass[i,]$y2<-sum(subset(subtab, (subtab$passes_filtering == TRUE))$sequence_length_template)

			summ[i,]$y3<-mean(subtab$sequence_length_template)
			summ_pass[i,]$y3<-mean(subset(subtab, (subtab$passes_filtering == TRUE))$sequence_length_template)

			summ[i,]$y4<-mean(subtab$mean_qscore_template)
			summ_pass[i,]$y4<-mean(subset(subtab, (subtab$passes_filtering == TRUE))$mean_qscore_template)

		}

		summ<-data.frame(summ)
		summ_pass<-data.frame(summ_pass)
		summ$z<-summ_pass$z<-run
		summ<-summ[complete.cases(summ),]
		summ_pass<-summ_pass[complete.cases(summ_pass),]
		taball[[run]]<-summ
		taball_pass[[run]]<-summ_pass

	}
	
	taball<-do.call(rbind, taball)
	taball_pass<-do.call(rbind,taball_pass)

	taball$x<-factor(taball$x, levels=unique(taball$x))
	taball_pass$x<-factor(taball_pass$x, levels=unique(taball_pass$x))


	listA<-as.list(rep(c(TRUE,FALSE),each=length(summaries)))
	listB<-as.list(rep(c(FALSE,TRUE),each=length(summaries)))

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
			args = list("visible", listA)),
		list(label = "pass",
			method = "restyle",
			args = list("visible", listB))
		)
	)


	p1 <- plot_ly() %>% add_trace(x = taball$x, y = taball$y, color = taball$z, type = "box", visible=TRUE,legendgroup = taball$z,showlegend=F) %>%
		add_trace(x=taball_pass$x,y=taball_pass$y,color = taball_pass$z, type = "box", visible=FALSE,legendgroup = taball_pass$z,showlegend=F) %>% 
		layout(boxmode = "group",legend = list(orientation = 'h'),updatemenus = list(chart_type), xaxis=list(title="Sequencing run-time (h)",titlefont=f,tickangle=45), 
			yaxis=list(side="left",title="#reads",titlefont=f))

	p2 <- plot_ly() %>% add_trace(x = taball$x, y = taball$y2, color = taball$z, type = "box", visible=TRUE,legendgroup = taball$z,showlegend=F) %>%
		add_trace(x=taball_pass$x,y=taball_pass$y2,color = taball_pass$z, type = "box", visible=FALSE,legendgroup=taball_pass$z,showlegend=F) %>% 
		layout(boxmode = "group",legend = list(orientation = 'h'),updatemenus = list(chart_type), xaxis=list(title="Sequencing run-time (h)",titlefont=f,tickangle=45),
			yaxis=list(side="left",title="#bases",titlefont=f))

	p3 <- plot_ly() %>% add_trace(x = taball$x, y = taball$y3, color = taball$z, type = "box", visible=TRUE,legendgroup = taball$z,showlegend=F) %>%
		add_trace(x=taball_pass$x,y=taball_pass$y3,color = taball_pass$z, type = "box", visible=FALSE,legendgroup=taball_pass$z,showlegend=F) %>% 
		layout(boxmode = "group",legend = list(orientation = 'h'),updatemenus = list(chart_type), xaxis=list(title="Sequencing run-time (h)",titlefont=f,tickangle=45),
			yaxis=list(side="left",title="Length",titlefont=f))

	p4<- plot_ly() %>% add_trace(x = taball$x, y = taball$y4, color = taball$z, type = "box", visible=TRUE,legendgroup = taball$z) %>%
		add_trace(x=taball_pass$x,y=taball_pass$y4,color = taball_pass$z, type = "box", visible=FALSE, legendgroup=taball_pass$z) %>% 
		layout(boxmode = "group",legend = list(orientation = 'h'),updatemenus = list(chart_type), xaxis=list(title="Sequencing run-time (h)",titlefont=f,tickangle=45),
			yaxis=list(side="left",title="Qscore",titlefont=f))

	s <- subplot(p1,p2,p3,p4,nrows = 2,shareX = TRUE,titleX=TRUE,titleY=TRUE)
	fig <- layout(s, hovermode = "x unified")


	message("[",Sys.time(),"]"," storing plot to file")

	htmlwidgets::saveWidget(fig, out)

	rm(list=ls())
	invisible(gc())

	message("[",Sys.time(),"]"," done")


}