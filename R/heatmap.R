#' @title heatmap
#' @description generate an interactive HTML report with channels activity over time and channels activity with respect to their disposition on the flowcell
#' @param summary path to a sequencing summary file from a single ONT run
#' @param time hours to split the sequencing run into. Default to 1
#' @param platfrom flowcell used for sequencing (minion, promethion. Default to minion)
#' @param out path to HTML file
#' @return HTML file
#' @examples
#' #do not run
#' summary<-file.path("/path/to/sequencing_summary.txt")
#' outhtml<-file.path("/path/to/heatmap.html")
#' NanoR::heatmap(summary=summary time=1,out=outhtml)


heatmap<-function(summary,time=1,platform="minion",out) {



	#the following functions were provided by Barbara Ottolini, fron Oxford Nanopore Technologies, pointing at https://github.com/sagrudd/nanopoRe/blob/master/R/FlowcellLayout.R
	#these should provide the channels layout in MinION and PromethION flowcells

	getPromethIONChannelMap <- function() {

	chunk <- function(i) {
		
		m <- matrix(seq_len(250), ncol=10, byrow=TRUE)
		m + i
	
	}
	
	layout <- do.call(cbind, lapply(seq(from=0, to=2750, by=250), chunk))
	channelMap <- as.data.frame(cbind(channel = as.vector(t(layout)), which(layout == as.vector(layout), arr.ind = TRUE)))
	return(channelMap)
	
	}


	#getFlongleChannelMap <- function() {

		#layout <- matrix(c(seq(1, 12), 0, seq(13, 24), 0, seq(25, 114), 0, seq(115, 126), 0), ncol = 13, byrow = TRUE)
		#layout <- layout[rev(seq(10)), ]
		#channelMap <- as.data.frame(cbind(channel = as.vector(t(layout)), which(layout == as.vector(layout), arr.ind = TRUE)))
		#return(channelMap)
	
	#}


	getMinIONChannelMap <- function() {

		# build the map for R9.4.1 flowcell, as a long-form dataframe
		blockCalc <- function(i) {
			
			m <- matrix(seq(i, i + 63, by = 1), ncol = 8, byrow = TRUE)
			cbind(m[seq(5, 8, by = 1), ], m[seq(4), rev(seq(8))])
		
		}
		
		layout <- do.call(rbind, lapply(c(1, 449, 385, 321, 257, 193, 129, 65), blockCalc))
		
		# transpose the layout for cleaner presentation ...
		#layout <- t(layout)
		channelMap <- as.data.frame(cbind(channel = as.vector(t(layout)), which( layout == as.vector(layout), arr.ind = TRUE)))
		return(channelMap)
	
	}

	
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

	if (platform == "minion") {

		n_channels<-512
		layout<-getMinIONChannelMap()

	#} else if (platform  == "flongle") { did not have the chance to test this.

		#n_channels<-126
		#layout<-getFlongleChannelMap()

	#}
	} else if (platform == "promethion")  {

		n_channels<-3000
		layout<-getPromethIONChannelMap()

	} else { #unknown platform

		stop("[",Sys.time(),"]"," wrong platform specified. Choose either minion or promethion")


	}

	channels_activity_overtime<-channels_activity_overtime_pass<-matrix(0,ncol=n_channels, nrow=length(bins)-1)
	channels_activity<-channels_activity_pass<-rep(0,n_channels)
	
	for (i in c(1:n_channels)) {
		
		subtab<-tab[channel==i]
		channels_activity[i]<-sum(subtab$sequence_length_template)
		channels_activity_pass[i]<-sum(subtab[passes_filtering == TRUE]$sequence_length_template)
		
		for (l in c(1:(length(bins)-1))) {
		
		from<-bins[l]
		to<-bins[l+1]
		subsubtab<-subtab[template_unix > from & template_unix <= to]
		channels_activity_overtime[l,i]<-sum(subsubtab$sequence_length_template)
		channels_activity_overtime_pass[l,i]<-sum(subsubtab[passes_filtering == TRUE]$sequence_length_template)
		
		}
		
	}

	channels_activity_labels<-matrix("0",nrow=max(layout$row),ncol=max(layout$col))
	channels_activity_map<-channels_activity_map_pass<-matrix(0,nrow=max(layout$row),ncol=max(layout$col))

	for (m in c(1:nrow(layout))) {

		r<-layout$row[m]
		c<-layout$col[m]
		label<-layout$channel[m]
		channels_activity_labels[r,c]<-as.character(label)
		channels_activity_map[r,c]<-channels_activity[m]
		channels_activity_map_pass[r,c]<-channels_activity_pass[m]

	}


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
			args = list("visible", list(TRUE,FALSE))),
		list(label = "pass",
			method = "restyle",
			args = list("visible", list(FALSE,TRUE)))
		)
	)

	message("[",Sys.time(),"]"," plotting")

	p1<-plot_ly() %>% add_trace(x=as.character(c(1:n_channels)),y=as.character(bins),z = channels_activity_overtime,type="heatmap", xgap=.4,ygap=.4,colors="OrRd",colorbar=list(x=1.02,y=1),hoverinfo="x+y+z") %>%
					add_trace(x=as.character(c(1:n_channels)),y=as.character(bins),z = channels_activity_overtime_pass,type="heatmap", xgap=.4,ygap=.4,colors="OrRd",colorbar=list(x=1.02,y=1),hoverinfo="x+y+z",visible=FALSE) %>%
					layout(yaxis=list(title="Sequencing run-time (h)",titlefont=f),
							xaxis=list(title="#channels",titlefont=f),
							showlegend=TRUE,updatemenus = list(chart_type))


	p2<-plot_ly() %>% add_trace(z = channels_activity_map, type = "heatmap", colors="OrRd", customdata=apply(channels_activity_labels,1,as.list),hovertemplate="%{customdata}<extra></extra>#bases: %{z}<extra></extra>", xgap=.8, ygap=.8,colorbar=list(x=1.02,y=.450)) %>%
					add_trace(z = channels_activity_map_pass, type = "heatmap", colors="OrRd", customdata=apply(channels_activity_labels,1,as.list),hovertemplate="%{customdata}<extra></extra>#bases: %{z}<extra></extra>", xgap=.8, ygap=.8,visible=FALSE,colorbar=list(x=1.02,y=.450)) %>%
					layout(yaxis=list(zeroline=FALSE, showline=FALSE, showticklabels=FALSE, showgrid=FALSE,ticks=""),
							xaxis=list(zeroline=FALSE, showline=FALSE, showticklabels=FALSE, showgrid=FALSE,ticks=""),
							showlegend=TRUE,updatemenus = list(chart_type), title="#bp per channel (over time and space)")

	fig<- subplot(p1,p2,nrows = 2,titleX=TRUE, titleY=TRUE,margin=.05)
	
	message("[",Sys.time(),"]"," storing plot to file")

	htmlwidgets::saveWidget(fig, out)

	rm(list=ls())
	invisible(gc())

	message("[",Sys.time(),"]"," done")


}
