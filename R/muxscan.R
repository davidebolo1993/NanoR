#' @title muxscan
#' @description generate an interactive HTML report with muxes status over time 
#' @param muxdata path to a mux scan data file from a single ONT run
#' @param out path to HTML report
#' @return HTML file
#' @examples
#' #do not run
#' muxdata<-file.path("/path/to/mux_scan_data.csv")
#' outhtml<-file.path("/path/to/muxscan.html")
#' NanoR::muxscan(muxdata=muxdata, out=outhtml)

muxscan<-function(muxdata,out) {

	muxdata<-normalizePath(file.path(muxdata))
	out<-file.path(out)

	message("[",Sys.time(),"]"," reading mux scan data")

	csv<-fread(muxdata, sep=",", header=TRUE, select=c("mux_scan_assessment", "repeat"), showProgress=FALSE)
	colnames(csv)<-c("A", "B")

	message("[",Sys.time(),"]"," calculating stats")

	csv<-as.data.frame(csv %>% group_by(B) %>% count(B, A))
	csv<-as.data.frame(complete(csv,B,A))
	csv[is.na(csv)] <- 0
	csv<-unstack(form=n ~ A, x=csv)
	csv<-setDT(csv, keep.rownames = TRUE)[]
	csv$rn<-factor(csv$rn, levels=as.character(csv$rn))


	f <- list(
		size = 10,
		color = "#7f7f7f"
	)

	message("[",Sys.time(),"]"," plotting")

	fig <- plot_ly(csv, x = ~rn, y = ~single_pore, type = "bar", name = "single_pore")
	fig <- fig %>% add_trace(y = ~saturated, name = "saturated")
	fig <- fig %>% add_trace(y = ~multiple, name = "multiple")
	fig <- fig %>% add_trace(y = ~other, name = "other")
	fig <- fig %>% add_trace(y = ~zero, name = "zero")
	fig <- fig %>% add_trace(y = ~unavailable, name = "unavailable")
	fig <- fig %>% layout(yaxis = list(title = "Count", titlefont=f), xaxis = list(title = "#mux scan", titlefont=f), barmode = "stack", title="Mux scan over time")

	message("[",Sys.time(),"]"," storing plot to file")

	htmlwidgets::saveWidget(fig,out)

	rm(list=ls())
	invisible(gc())

	message("[",Sys.time(),"]"," done")


}