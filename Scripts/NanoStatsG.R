
############################################ NanoStatsG ############################################


#' @title Plots statistics for your GridION X5 .fast5 files 
#' @description NanoStatsG plots statistics for sequences with quality greater (or equal) to 7 returned from GridION X5 and outputs 4 tables required by NanoCompare.
#' @param NanoPrepareGList The object of class "list" returned by NanoPrepareG function
#' @param NanoGTable The table returned by NanoTableG function
#' @param DataOut Where NanoStatsG results will be saved. Use the same directory specified for NanoTableG function and be sure that it doesn't already contain NanoStatsM/NanoStatsG results
#' @return Plots: \cr 
#' - Cumulative reads and cumulative base pairs (Cumulative_Reads_&_Cumulative_Basepairs.pdf); \cr 
#' - Reads number, base pairs number, reads length (min, max, mean), reads quality (min, max, mean) calculated every half an hour (Reads_Basepairs_Length_Quality.pdf); \cr 
#' - Reads length versus reads quality (Length_versus_Quality.pdf); \cr 
#' - "Pass", "Fail" and "Skip" reads number, percentage and reads GC content count (if "GCC = TRUE" is provided to NanoTableG function) (Pass_Fail_Skip_and_GC_Content.pdf / Pass_Fail_Skip_NO_GC_Content.pdf); \cr 
#' - Channels activity (and muxes activity if analyzing basecalled .fast5 files) (Channels_Activity.pdf / Channels_and_Muxes_Activity.pdf) with respect to their REAL disposition in the flowcell (more at: https://community.nanoporetech.com/technical_documents/hardware/v/hwtd_5000_v1_revh_03may2016/flow-cell-chip). Not-working channels and muxes are grey-colored.
#' @examples
#' #do not run
#' #knows how to deal with different inputs type autonomously 
#' NanoStatsG(NanoPrepareGList=NanoGList, NanoGTable=NanoGTable, DataOut="/Path/To/DataOutEx")



NanoStatsG<-function(NanoPrepareGList,NanoGTable,DataOut) {
  
  library(reshape2)
  library(scales)
  library(ggplot2)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
  
  
  Increment <- function(x)
  {
    return(x+8)
  }
  
  Increment2 <- function(x)
  {
    return(x+32)
  }
  
  Increment3<-function(x)
  {
    return(x+128)
  }
  
  rotate<-function(x) 
  { 
    t(apply(x, 1, rev))
  }
  
  define_region <- function(row,col)
  {
    viewport(layout.pos.row = row, layout.pos.col = col)
  }
  
  
  
  Directory<-file.path(DataOut)
  dir.create(Directory, showWarnings = FALSE)
  setwd(Directory)
  TableInDirectory<-list.files(Directory,pattern="Information_Table")

  if(length(TableInDirectory) == 0) {
    stop("Use the same directory specified for NanoTableG function")
  }
  
  label<-as.character(unlist(strsplit(TableInDirectory,"_"))[2])


  CumulativeInDirectory<-list.files(Directory,pattern=("Cumulative_Basepairs.pdf"))
  
  if(length(CumulativeInDirectory) != 0) {
    stop("Can't use a directory that already contains NanoStatsM/NanoStatsG results")
  }
  
  options(scipen=9999)
  Label<-NanoPrepareGList[[6]]


  
  if (is.na(NanoPrepareGList[[1]][1]) == FALSE & is.na(NanoPrepareGList[[3]][1])) {
    
    
    
    
    Really_Pass_File<-which(as.numeric(NanoGTable[,6]) >= 7)
    Fail_File_Length<-length(which(as.numeric(NanoGTable[,6]) < 7))
    NanoTable2<-NanoGTable[Really_Pass_File,]
    List.Files.HDF5_Fail_Length<-as.numeric(NanoPrepareGList[[4]])+Fail_File_Length
    List.Files.HDF5_Pass.length<-dim(NanoGTable)[1]-Fail_File_Length
    List.Files.HDF5_Skip_Length<-as.numeric(NanoPrepareGList[[5]])
    
    
    write.table(NanoTable2, file.path(Directory, paste0(Label,"_",label, "_Information_Table.txt")), col.names=T, row.names=F, quote=F, sep="\t")
    
    
    Table_HDF5<-NanoTable2[,1:6]
    row.names(Table_HDF5)<-c()
    
    Time_2<-as.numeric(Table_HDF5[,4])
    Run_Duration<-round(as.numeric(difftime(as.POSIXct(max(Time_2),origin="1/1/1970"), as.POSIXct(min(Time_2), origin="1/1/1970"), units="hours")))
    Time_Rescaled <- scales::rescale(Time_2, to=c(0,Run_Duration))
    
    Table_HDF5_Def<-cbind(Table_HDF5,Time_Rescaled)
    colnames(Table_HDF5_Def)<-c("Flowcell ID","Read Id","Channel Number","Relative Time","Length of Read","Quality","Relative Experimental Time")
    
    
    Relative_Time<-as.numeric(Table_HDF5_Def[,7])
    Relative_Time_Per_Hours<-seq(from=min(round(Relative_Time)), to=max(round(Relative_Time)), by=0.5)
    Template_Length<-as.numeric(Table_HDF5_Def[,5])
    Quality_Score<-as.numeric(Table_HDF5_Def[,6])
    
    
    message("Analyzing...")
    
    Reads_Per_Hour<-c()
    Base_Pairs_Per_Hour<-c()
    Max_Length_Per_Hour<-c()
    Mean_Length_Per_Hour<-c()
    Min_Length_Per_Hour<-c()
    Min_Quality_Score_Per_Hour<-c()
    Mean_Quality_Score_Per_Hour<-c()
    Max_Quality_Score_Per_Hour<-c()
    
    
    
    for (ii in 1:(length(Relative_Time_Per_Hours))) {
      
      if (ii < length(Relative_Time_Per_Hours)) {
        Index_Hours<-which(Relative_Time >= Relative_Time_Per_Hours[ii] & Relative_Time < Relative_Time_Per_Hours[ii+1])
        if (length(Index_Hours) == 0) {
          Reads_Per_Hour[ii]<-0
          Base_Pairs_Per_Hour[ii]<-0
          Mean_Length_Per_Hour[ii]<-0
          Max_Length_Per_Hour[ii]<-0
          Min_Length_Per_Hour[ii]<-0
          Mean_Quality_Score_Per_Hour[ii]<-0
          Min_Quality_Score_Per_Hour[ii]<-0
          Max_Quality_Score_Per_Hour[ii]<-0 
        }
        else
        {
          Reads_Per_Hour[ii]<-length(Index_Hours)
          Base_Pairs_Per_Hour[ii]<-sum(Template_Length[Index_Hours])
          Mean_Length_Per_Hour[ii]<-mean(Template_Length[Index_Hours])
          Max_Length_Per_Hour[ii]<-max(Template_Length[Index_Hours])
          Min_Length_Per_Hour[ii]<-min(Template_Length[Index_Hours])
          Mean_Quality_Score_Per_Hour[ii]<-mean(Quality_Score[Index_Hours])
          Min_Quality_Score_Per_Hour[ii]<-min(Quality_Score[Index_Hours])
          Max_Quality_Score_Per_Hour[ii]<-max(Quality_Score[Index_Hours])
        }
      }
      else {
        Index_Hours<-which(Relative_Time == Relative_Time_Per_Hours[ii])
        
        if (length(Index_Hours) == 0) {
          Reads_Per_Hour[ii]<-0
          Base_Pairs_Per_Hour[ii]<-0
          Mean_Length_Per_Hour[ii]<-0
          Max_Length_Per_Hour[ii]<-0
          Min_Length_Per_Hour[ii]<-0
          Mean_Quality_Score_Per_Hour[ii]<-0
          Min_Quality_Score_Per_Hour[ii]<-0
          Max_Quality_Score_Per_Hour[ii]<-0   
        }
        else {
          Reads_Per_Hour[ii]<-length(Index_Hours)
          Base_Pairs_Per_Hour[ii]<-sum(Template_Length[Index_Hours])
          Mean_Length_Per_Hour[ii]<-mean(Template_Length[Index_Hours])
          Max_Length_Per_Hour[ii]<-max(Template_Length[Index_Hours])
          Min_Length_Per_Hour[ii]<-min(Template_Length[Index_Hours])
          Mean_Quality_Score_Per_Hour[ii]<-mean(Quality_Score[Index_Hours])
          Min_Quality_Score_Per_Hour[ii]<-min(Quality_Score[Index_Hours])
          Max_Quality_Score_Per_Hour[ii]<-max(Quality_Score[Index_Hours])
        }
      }
    }
    
    
    Cumulative_Reads<-cumsum(Reads_Per_Hour)
    Cumulative_Basepairs<-cumsum(Base_Pairs_Per_Hour)
    
    Channel_Vector<-as.numeric(Table_HDF5_Def[,2])
    Channels_Number<-c(1:512)
    
    
    Base_Pairs_Per_Channel<-c()
    
    Table_HDF5_Reordered<-c()
    
    for (iii in 1:length(Channels_Number)) {
      Ind_Chann<-which(Channel_Vector == Channels_Number[iii])
      if (length(Ind_Chann) == 0) {
        next
      }
      else {
        Table_HDF5_Re<-Table_HDF5_Def[Ind_Chann,]
        Table_HDF5_Reordered<-rbind(Table_HDF5_Reordered,Table_HDF5_Re)
      }
      Base_Pairs_Per_Channel[iii]<-sum(Template_Length[Ind_Chann])
    }
    
    rownames(Table_HDF5_Reordered)<-c()
    
    
    message("Plotting...")
    
    x<-Relative_Time_Per_Hours
    y0.1<-Cumulative_Reads
    data0.1<-data.frame('x'=x,'y'=y0.1)
    data0.1$group<-"Cumulative Reads"
    
    
    Cumulative_Reads_Plot<-ggplot(data0.1, aes(x=x, y=y, group=group, fill=group)) +
      geom_line(size=.5) + 
      scale_x_continuous(name="Experimental Time", breaks=(seq(0,Run_Duration,2)))+
      scale_y_continuous(name="Number Of Reads")+
      geom_ribbon(data=subset(data0.1,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0,show.legend=FALSE) +
      scale_fill_manual(name='', values=c("Cumulative Reads" = "coral3"))+
      theme_bw(base_size = 12)+
      theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=14),legend.background = (element_rect(size=0.5, linetype="solid",colour ="black")), legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=12,face="italic"),axis.title.x = element_text(size=11, face="italic"),axis.title.y = element_text(size=11, face="italic"))+
      theme(legend.position="bottom")+
      ggtitle("Cumulative Reads")
    
    y0.2<-Cumulative_Basepairs
    data0.2<-data.frame('x'=x,'y'=y0.2)
    data0.2$group<-"Cumulative Base Pairs"
    
    Cumulative_Base_Pairs_Plot<-ggplot(data0.2, aes(x=x, y=y, group=group, fill=group)) +
      geom_line(size=.5) + 
      scale_x_continuous(name="Experimental Time", breaks=(seq(0,Run_Duration,2)))+
      scale_y_continuous(name="Number Of Base Pairs")+
      geom_ribbon(data=subset(data0.2,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0,show.legend=FALSE) +
      scale_fill_manual(name='', values=c("Cumulative Base Pairs" = "darkcyan"))+
      theme_bw(base_size = 12)+
      theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=14),legend.background = (element_rect(size=0.5, linetype="solid",colour ="black")), legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=12,face="italic"),axis.title.x = element_text(size=11, face="italic"),axis.title.y = element_text(size=11, face="italic"))+
      theme(legend.position="bottom")+
      ggtitle("Cumulative Base Pairs")
    
    Cumulative_Plot<-arrangeGrob(grid.arrange(Cumulative_Reads_Plot,Cumulative_Base_Pairs_Plot, nrow=2, ncol=1))
    
    ggsave(paste0(Label,"_",label,"_Cumulative_Reads_&_Cumulative_Basepairs.pdf"), device="pdf", Cumulative_Plot, height=10,width=15)
    
    #PLOT PER-HOUR READS/BPs/QUALITY/LENGTH
    
    y1<-Reads_Per_Hour
    data1<-data.frame('x'=x,'y'=y1)
    data1$group<-"Reads Per Hour"
    
    Reads_Per_Hour_Plot<-ggplot(data1, aes(x=x, y=y, group=group, fill=group)) +
      geom_line(size=.5) + 
      scale_x_continuous(name="Experimental Time", breaks=(seq(0,Run_Duration,2)))+
      scale_y_continuous(name="Number Of Reads")+
      geom_ribbon(data=subset(data1,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0,show.legend=FALSE) +
      scale_fill_manual(name='', values=c("Reads Per Hour" = "coral3"))+
      theme_bw(base_size = 12)+
      theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=14),legend.background = (element_rect(size=0.5, linetype="solid",colour ="black")), legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=12),axis.title.x = element_text(size=11, face="italic"),axis.title.y = element_text(size=11, face="italic"))+
      theme(legend.position="bottom")+
      ggtitle("Number Of Reads")
    
    
    y2<-Base_Pairs_Per_Hour
    data2<-data.frame('x'=x,'y'=y2)
    data2$group<-"Base Pairs Per Hour"
    
    
    Base_Pairs_Per_Hour_Plot<-ggplot(data2, aes(x=x, y=y, group=group, fill=group)) +
      geom_line(size=.5) + 
      scale_x_continuous(name="Experimental Time", breaks=(seq(0,Run_Duration,2)))+
      scale_y_continuous(name="Number Of Base Pairs")+
      geom_ribbon(data=subset(data2,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0,show.legend=FALSE) +
      scale_fill_manual(name='', values=c("Base Pairs Per Hour" = "darkcyan"))+
      theme_bw(base_size = 12)+
      theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=14),legend.background = (element_rect(size=0.5, linetype="solid",colour ="black")), legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=12),axis.title.x = element_text(size=11, face="italic"),axis.title.y = element_text(size=11, face="italic"))+
      theme(legend.position="bottom")+
      ggtitle("Number Of Base Pairs")
    
    
    y3.0<-log10(Mean_Length_Per_Hour)
    data3.0.0<-data.frame('x'=x,'y'=Mean_Length_Per_Hour)
    data3.0<-data.frame('x'=x,'y'=y3.0)
    data3.0$group<-"Mean Length"
    
    y3.1<-log10(Max_Length_Per_Hour)
    data3.1<-data.frame('x'=x,'y'=y3.1)
    data3.1$group<-"Max Length"
    
    y3.2<-log10(Min_Length_Per_Hour)
    data3.2<-data.frame('x'=x,'y'=y3.2)
    data3.2$group<-"Min Length"
    
    data3<-rbind(data3.0, data3.1, data3.2)
    
    Length_Per_Hour_Plot<-ggplot(data3, aes(x=x, y=y, group=group, fill=group)) +
      geom_line(size=.5) + 
      scale_x_continuous(name="Experimental Time", breaks=(seq(0,Run_Duration,2)))+
      scale_y_continuous(name="Length", breaks=c(2,3,4,5,6), labels=c("100","1000","10000","100000","1000000"))+
      geom_ribbon(data=subset(data3,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0,alpha=.7) +
      scale_fill_manual(name='', values=c("Mean Length" = "darkolivegreen1", "Min Length" = "darkolivegreen", "Max Length" = "cornsilk"))+
      theme_bw(base_size = 12)+
      theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=14),legend.background = (element_rect(size=0.5, linetype="solid",colour ="black")), legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=12),axis.title.x = element_text(size=11, face="italic"),axis.title.y = element_text(size=11, face="italic"))+
      theme(legend.position="bottom")+
      ggtitle("Length")
    
    y4.0<-Mean_Quality_Score_Per_Hour
    data4.0<-data.frame('x'=x,'y'=y4.0)
    data4.0$group<-"Mean Quality"
    
    y4.1<-Min_Quality_Score_Per_Hour
    data4.1<-data.frame('x'=x,'y'=y4.1)
    data4.1$group<-"Min Quality"
    
    y4.2<-Max_Quality_Score_Per_Hour
    data4.2<-data.frame('x'=x,'y'=y4.2)
    data4.2$group<-"Max Quality"
    
    data4<-rbind(data4.0,data4.1,data4.2)
    
    Quality_Score_Per_Hour_Plot<-ggplot(data4, aes(x=x, y=y, group=group, fill=group)) +
      geom_line(size=.5) + 
      scale_x_continuous(name="Experimental Time", breaks=(seq(0,Run_Duration,2)))+
      scale_y_continuous(name="Quality")+
      geom_ribbon(data=subset(data4,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0,alpha=.7) +
      scale_fill_manual(name='', values=c("Mean Quality" = "chocolate1", "Min Quality" = "chocolate", "Max Quality" = "cornsilk"))+
      theme_bw(base_size = 12)+
      theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=14),legend.background = (element_rect(size=0.5, linetype="solid",colour ="black")), legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=12),axis.title.x = element_text(size=11, face="italic"),axis.title.y = element_text(size=11, face="italic"))+
      theme(legend.position="bottom")+
      ggtitle("Quality Score")
    
    
    Others_Plot<-arrangeGrob(grid.arrange(Reads_Per_Hour_Plot,Base_Pairs_Per_Hour_Plot,Length_Per_Hour_Plot,Quality_Score_Per_Hour_Plot, nrow=2, ncol=2))
    
    ggsave(paste0(Label,"_",label,"_Reads_Basepairs_Length_Quality.pdf"), device="pdf", Others_Plot, height=10,width=15)
    
    
    #PASS/FAIL
    
    blank_theme <- theme_minimal()+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_blank()
      )
    
    Data_Pass_Fail_Percentage <- data.frame(
      group = c("Pass", "Fail/Skip"),
      value = c((List.Files.HDF5_Pass.length/(List.Files.HDF5_Pass.length+List.Files.HDF5_Fail_Length+List.Files.HDF5_Skip_Length)), ((List.Files.HDF5_Fail_Length+List.Files.HDF5_Skip_Length)/(List.Files.HDF5_Pass.length+List.Files.HDF5_Fail_Length+List.Files.HDF5_Skip_Length)))
    )
    
    
    Data_Pass_Fail_Percentage_Plot<-ggplot(Data_Pass_Fail_Percentage, aes(x="", y=value, fill=group))+
      geom_bar(width = 1, size = 1, color = "white", stat = "identity")+
      geom_text(aes(label = scales::percent(value)), position = position_stack(vjust = 0.5))+
      coord_polar("y", start=0) +
      scale_fill_brewer(palette="Dark2") +
      blank_theme+
      theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=14),legend.title=element_blank(),legend.position="bottom")+
      ggtitle("Passed and Failed/Skipped Percentage")
    
    
    Data_Pass_Fail_Tot <- data.frame(
      group = c("Pass", "Fail/Skip"),
      value = c(List.Files.HDF5_Pass.length, (List.Files.HDF5_Fail_Length+List.Files.HDF5_Skip_Length))
    )
    
    Data_Pass_Fail_Tot_Plot<-ggplot(Data_Pass_Fail_Tot, aes(x="", y=value, fill=group))+
      geom_bar(width = 1, size = 1, color = "white", stat = "identity")+
      geom_text(aes(label = value), position = position_stack(vjust = 0.5))+
      coord_polar("y", start=0) +
      scale_fill_brewer(palette="Accent") +
      blank_theme+
      theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=14),legend.title=element_blank(), legend.position="bottom")+
      ggtitle("Passed and Failed/Skipped Number")
    
    
    ###################### LENGTH VS QUALITY ######################
    
    Tot_Length<-data.frame(Template_Length)
    Tot_Quality<-data.frame(Quality_Score)
    
    
    ScatterTheme <- list(labs(x="Length",y="Quality"),theme_bw(), theme(legend.position=c(1,0),legend.justification=c(1,0), legend.background=element_blank(),legend.direction="horizontal", legend.title=element_text(face="bold.italic"),axis.title=element_text(face="italic")))
    
    
    hist_top_mean_length<-ggplot(Tot_Length, aes(x=Template_Length))+theme_bw()+ geom_histogram(aes(y = ..count../1000),col="darkolivegreen", fill="darkolivegreen1",boundary = min(Tot_Length$Template_Length))+labs(x="",y="Count (x 1e3)")+scale_x_continuous(limits=c(min(Tot_Length$Template_Length),max(Tot_Length$Template_Length)))+theme(axis.title=element_text(face="italic"))
    hist_right_mean_quality<-ggplot(Tot_Quality, aes(x=Quality_Score))+ theme_bw()+ geom_histogram(aes(y = ..count../1000),col="chocolate", fill="chocolate1",boundary = min(Tot_Quality$Quality_Score))+ labs(x="",y="Count (x 1e3)")+coord_flip()+scale_x_continuous(limits=c(min(Tot_Quality$Quality_Score),max(Tot_Quality$Quality_Score)))+theme(axis.title=element_text(face="italic"))
    empty <- ggplot()+geom_point(aes(1,1), colour="white")+
      theme(axis.ticks=element_blank(), 
            panel.background=element_blank(), 
            axis.text.x=element_blank(), axis.text.y=element_blank(),           
            axis.title.x=element_blank(), axis.title.y=element_blank())
    scatter <- ggplot(data.frame(cbind(Tot_Length$Template_Length,Tot_Quality$Quality_Score)),aes(x=Tot_Length$Template_Length, y=Tot_Quality$Quality_Score))+geom_point(col="grey27", size=.09, alpha=.4)+scale_x_continuous(limits=c(min(Tot_Length$Template_Length),max(Tot_Length$Template_Length)))+scale_y_continuous(limits=c(min(Tot_Quality$Quality_Score),max(Tot_Quality$Quality_Score)))+stat_density2d(aes(col=..level.., alpha=..level..)) + scale_color_continuous(low="darkblue",high="darkred") +geom_smooth(method=lm,linetype=2, size=.5,col="black",se=F) + guides(alpha="none",col=guide_legend(title="Density"))+ ScatterTheme
    
    
    
    
    Length_VS_Quality_Plot<-grid.arrange(hist_top_mean_length, empty, scatter, hist_right_mean_quality, ncol=2, nrow=2, widths=c(4,1), heights=c(1, 4))
    
    
    ggsave(paste0(Label,"_",label,"_Length_versus_Quality.pdf"), device="pdf", Length_VS_Quality_Plot, height=10,width=15)
    
    
    
    ### skip muxes rapresentation
    
    
    #PLOT CORRELATION MATRIXES (CHANNEL AND MUXES)
    
    m1<-matrix(Base_Pairs_Per_Channel[1:32], ncol=8, nrow=4, byrow=TRUE)
    m2<-matrix(Base_Pairs_Per_Channel[449:480], ncol=8, nrow=4, byrow=TRUE)
    m3<-matrix(Base_Pairs_Per_Channel[385:416], ncol=8, nrow=4, byrow=TRUE)
    m4<-matrix(Base_Pairs_Per_Channel[321:352], ncol=8, nrow=4, byrow=TRUE)
    m5<-matrix(Base_Pairs_Per_Channel[257:288], ncol=8, nrow=4, byrow=TRUE)
    m6<-matrix(Base_Pairs_Per_Channel[193:224], ncol=8, nrow=4, byrow=TRUE)
    m7<-matrix(Base_Pairs_Per_Channel[129:160], ncol=8, nrow=4, byrow=TRUE)
    m8<-matrix(Base_Pairs_Per_Channel[65:96], ncol=8, nrow=4, byrow=TRUE)
    mdef3<-rbind(m1,m2,m3,m4,m5,m6,m7,m8)
    m9<-rotate(matrix(Base_Pairs_Per_Channel[33:64], ncol=8, nrow=4, byrow=TRUE))
    m10<-rotate(matrix(Base_Pairs_Per_Channel[481:512], ncol=8, nrow=4, byrow=TRUE))
    m11<-rotate(matrix(Base_Pairs_Per_Channel[417:448], ncol=8, nrow=4, byrow=TRUE))
    m12<-rotate(matrix(Base_Pairs_Per_Channel[353:384], ncol=8, nrow=4, byrow=TRUE))
    m13<-rotate(matrix(Base_Pairs_Per_Channel[289:320], ncol=8, nrow=4, byrow=TRUE))
    m14<-rotate(matrix(Base_Pairs_Per_Channel[225:256], ncol=8, nrow=4, byrow=TRUE))
    m15<-rotate(matrix(Base_Pairs_Per_Channel[161:192], ncol=8, nrow=4, byrow=TRUE))
    m16<-rotate(matrix(Base_Pairs_Per_Channel[97:128], ncol=8, nrow=4, byrow=TRUE))
    mdef4<-rbind(m9,m10,m11,m12,m13,m14,m15,m16)
    Matrixbpchannel<-cbind(mdef3,mdef4)
    
    
    #PLOTTING "FALSE" CORRELATION MATRIXES
    
    Palette <- colorRampPalette(brewer.pal(9, "Reds"))
    
    adjMatrixbpchannel<-melt(rotate(t(Matrixbpchannel)))
    Plot_Channel_Activity<-ggplot(data=adjMatrixbpchannel, aes(x=Var1, y=Var2)) +
      geom_tile(aes(fill=value), color="white", size=2)+
      scale_fill_gradientn(colours= Palette(4), na.value="grey70",limits=c(min(adjMatrixbpchannel[,3], na.rm=TRUE), max(adjMatrixbpchannel[,3], na.rm=TRUE))) +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="bottom",
            legend.title=element_text(face="bold.italic"),
            panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),plot.background=element_blank(),legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5,face="bold.italic",size=16),legend.text.align = 0)+
      guides(fill = guide_colorbar(barwidth= 61,barheight=.5,title="Base Pairs Productivity", title.position="top",title.hjust=0.5))+
      ggtitle("Channels Activity")
    
    ggsave(paste0(Label,"_",label,"_Channels_Activity.pdf"), device="pdf",Plot_Channel_Activity, height=10, width=18)
    
    
    ###PLOT GC CONTENT#######
    
    
    
    if (NanoTable2[1,7] == "GC_Content") {
      
      pdf(paste0(Label,"_",label,"_Pass_Fail_Skip_NO_GC_Content.pdf"), height=10, width=15, onefile=TRUE)
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(1, 2)))
      print(Data_Pass_Fail_Percentage_Plot, vp=define_region(1, 1))
      print(Data_Pass_Fail_Tot_Plot, vp = define_region(1, 2))
      dev.off()
    }
    
    else {
      
      GC_Content_To_Plot<-as.numeric(NanoTable2[,7])
      Hist_GC_Content<-ggplot(data.frame(GC_Content_To_Plot), aes(x=GC_Content_To_Plot))+theme_bw()+ geom_histogram(col="midnightblue", fill="cyan4")+labs(x="",y="Count")+ggtitle("GC Content")+theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=16), axis.title.y=element_text(face="italic"))
      pdf(paste0(Label,"_",label,"_Pass_Fail_Skip_and_GC_Content.pdf"), height=10, width=15, onefile=TRUE)
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(2, 2)))
      print(Data_Pass_Fail_Percentage_Plot, vp=define_region(1, 1))
      print(Data_Pass_Fail_Tot_Plot, vp = define_region(1, 2))
      print(Hist_GC_Content, vp = define_region(2, 1:2))
      dev.off()
    }
    
    #save data for other analyses
    
    write.table(data1, file.path(Directory, paste0(Label,"_",label, "_ReadsProduced.txt")),col.names=T, sep="\t")
    write.table(data2, file.path(Directory, paste0(Label,"_",label, "_BasePairsProduced.txt")),col.names=T, sep="\t")
    write.table(data3.0.0, file.path(Directory, paste0(Label,"_",label, "_MeanLength.txt")),col.names=T, sep="\t")
    write.table(data4.0, file.path(Directory, paste0(Label,"_",label, "_MeanQuality.txt")),col.names=T, sep="\t")
    
    message("Done!")
    
  }
  
  else {
    
    
    NanoGTable<-NanoGTable[,2:8]
    FilesAnalyzed<-which(as.character(NanoGTable[,1]) != "Read_Id")
    Length_Not_Analyzed<-length(which(as.character(NanoGTable[,1]) == "Read_Id"))
    NanoTable<-NanoGTable[FilesAnalyzed,]
    Really_Pass_File<-which(as.numeric(NanoTable[,6]) >= 7)
    Fail_File_Length<-length(which(as.numeric(NanoTable[,6]) < 7))
    NanoTable<-NanoTable[Really_Pass_File,]
    List.Files.HDF5_Fail_Length<-as.numeric(NanoPrepareGList[[4]])+Length_Not_Analyzed+Fail_File_Length
    List.Files.HDF5_Pass.length<-length(NanoPrepareGList[[3]])-Length_Not_Analyzed-Fail_File_Length
    List.Files.HDF5_Skip_Length<-as.numeric(NanoPrepareGList[[5]])
    
    write.table(NanoTable, file.path(Directory, paste0(Label,"_",label, "_Information_Table.txt")), col.names=T, row.names=F, quote=F, sep="\t")
    
    Table_HDF5<-NanoTable[,1:6]
    
    Time_2<-as.numeric(Table_HDF5[,4])
    
    Run_Duration<-round(as.numeric(difftime(as.POSIXct(max(Time_2),origin="1/1/1970"), as.POSIXct(min(Time_2), origin="1/1/1970"), units="hours")))
    
    Time_Rescaled <- scales::rescale(Time_2, to=c(0,Run_Duration))
    Table_HDF5_Def<-cbind(Table_HDF5,Time_Rescaled)
    
    colnames(Table_HDF5_Def)<-c("Read", "Channel Number", "Mux Number", "Unix Time", "Length of Read", "Quality", "Relative Experimental Time")
    
    
    Relative_Time<-as.numeric(Table_HDF5_Def[,7])
    Relative_Time_Per_Hours<-seq(from=min(round(Relative_Time)), to=max(round(Relative_Time)), by=0.5)
    Template_Length<-as.numeric(Table_HDF5_Def[,5])
    Quality_Score<-as.numeric(Table_HDF5_Def[,6])
    
    
    message("Analyzing...")
    
    Reads_Per_Hour<-c()
    Base_Pairs_Per_Hour<-c()
    Max_Length_Per_Hour<-c()
    Mean_Length_Per_Hour<-c()
    Min_Length_Per_Hour<-c()
    Min_Quality_Score_Per_Hour<-c()
    Mean_Quality_Score_Per_Hour<-c()
    Max_Quality_Score_Per_Hour<-c()
    
    
    for (ii in 1:(length(Relative_Time_Per_Hours))) {
      
      if (ii < length(Relative_Time_Per_Hours)) {
        Index_Hours<-which(Relative_Time >= Relative_Time_Per_Hours[ii] & Relative_Time < Relative_Time_Per_Hours[ii+1])
        if (length(Index_Hours) == 0) {
          Reads_Per_Hour[ii]<-0
          Base_Pairs_Per_Hour[ii]<-0
          Mean_Length_Per_Hour[ii]<-0
          Max_Length_Per_Hour[ii]<-0
          Min_Length_Per_Hour[ii]<-0
          Mean_Quality_Score_Per_Hour[ii]<-0
          Min_Quality_Score_Per_Hour[ii]<-0
          Max_Quality_Score_Per_Hour[ii]<-0 
        }
        else
        {
          Reads_Per_Hour[ii]<-length(Index_Hours)
          Base_Pairs_Per_Hour[ii]<-sum(Template_Length[Index_Hours])
          Mean_Length_Per_Hour[ii]<-mean(Template_Length[Index_Hours])
          Max_Length_Per_Hour[ii]<-max(Template_Length[Index_Hours])
          Min_Length_Per_Hour[ii]<-min(Template_Length[Index_Hours])
          Mean_Quality_Score_Per_Hour[ii]<-mean(Quality_Score[Index_Hours])
          Min_Quality_Score_Per_Hour[ii]<-min(Quality_Score[Index_Hours])
          Max_Quality_Score_Per_Hour[ii]<-max(Quality_Score[Index_Hours])
        }
      }
      else {
        Index_Hours<-which(Relative_Time == Relative_Time_Per_Hours[ii])
        
        if (length(Index_Hours) == 0) {
          Reads_Per_Hour[ii]<-0
          Base_Pairs_Per_Hour[ii]<-0
          Mean_Length_Per_Hour[ii]<-0
          Max_Length_Per_Hour[ii]<-0
          Min_Length_Per_Hour[ii]<-0
          Mean_Quality_Score_Per_Hour[ii]<-0
          Min_Quality_Score_Per_Hour[ii]<-0
          Max_Quality_Score_Per_Hour[ii]<-0   
        }
        else {
          Reads_Per_Hour[ii]<-length(Index_Hours)
          Base_Pairs_Per_Hour[ii]<-sum(Template_Length[Index_Hours])
          Mean_Length_Per_Hour[ii]<-mean(Template_Length[Index_Hours])
          Max_Length_Per_Hour[ii]<-max(Template_Length[Index_Hours])
          Min_Length_Per_Hour[ii]<-min(Template_Length[Index_Hours])
          Mean_Quality_Score_Per_Hour[ii]<-mean(Quality_Score[Index_Hours])
          Min_Quality_Score_Per_Hour[ii]<-min(Quality_Score[Index_Hours])
          Max_Quality_Score_Per_Hour[ii]<-max(Quality_Score[Index_Hours])
        }
      }
    }
    
    
    Cumulative_Reads<-cumsum(Reads_Per_Hour)
    Cumulative_Basepairs<-cumsum(Base_Pairs_Per_Hour)
    
    Channel_Vector<-as.numeric(Table_HDF5_Def[,2])
    Mux_Vector<-as.numeric(Table_HDF5_Def[,3])
    Channels_Number<-c(1:512)
    
    
    Base_Pairs_Per_Channel<-c()
    
    Table_HDF5_Reordered<-c()
    
    for (iii in 1:length(Channels_Number)) {
      Ind_Chann<-which(Channel_Vector == Channels_Number[iii])
      Mux_Associated<-sort(Mux_Vector[Ind_Chann], index.return=TRUE)$ix
      if (length(Ind_Chann) == 0) {
        next
      }
      else {
        Table_HDF5_Re<-Table_HDF5_Def[Ind_Chann,][Mux_Associated,]
        Table_HDF5_Reordered<-rbind(Table_HDF5_Reordered,Table_HDF5_Re)
      }
      Base_Pairs_Per_Channel[iii]<-sum(Template_Length[Ind_Chann])
    }
    
    rownames(Table_HDF5_Reordered)<-c()
    
    
    #PLOT CUMULATIVE READS/BP
    
    
    message("Plotting...")
    
    x<-Relative_Time_Per_Hours
    y0.1<-Cumulative_Reads
    data0.1<-data.frame('x'=x,'y'=y0.1)
    data0.1$group<-"Cumulative Reads"
    
    
    Cumulative_Reads_Plot<-ggplot(data0.1, aes(x=x, y=y, group=group, fill=group)) +
      geom_line(size=.5) + 
      scale_x_continuous(name="Experimental Time", breaks=(seq(0,Run_Duration,2)))+
      scale_y_continuous(name="Number Of Reads")+
      geom_ribbon(data=subset(data0.1,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0,show.legend=FALSE) +
      scale_fill_manual(name='', values=c("Cumulative Reads" = "coral3"))+
      theme_bw(base_size = 12)+
      theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=14),legend.background = (element_rect(size=0.5, linetype="solid",colour ="black")), legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=12,face="italic"),axis.title.x = element_text(size=11, face="italic"),axis.title.y = element_text(size=11, face="italic"))+
      theme(legend.position="bottom")+
      ggtitle("Cumulative Reads")
    
    y0.2<-Cumulative_Basepairs
    data0.2<-data.frame('x'=x,'y'=y0.2)
    data0.2$group<-"Cumulative Base Pairs"
    
    Cumulative_Base_Pairs_Plot<-ggplot(data0.2, aes(x=x, y=y, group=group, fill=group)) +
      geom_line(size=.5) + 
      scale_x_continuous(name="Experimental Time", breaks=(seq(0,Run_Duration,2)))+
      scale_y_continuous(name="Number Of Base Pairs")+
      geom_ribbon(data=subset(data0.2,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0,show.legend=FALSE) +
      scale_fill_manual(name='', values=c("Cumulative Base Pairs" = "darkcyan"))+
      theme_bw(base_size = 12)+
      theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=14),legend.background = (element_rect(size=0.5, linetype="solid",colour ="black")), legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=12,face="italic"),axis.title.x = element_text(size=11, face="italic"),axis.title.y = element_text(size=11, face="italic"))+
      theme(legend.position="bottom")+
      ggtitle("Cumulative Base Pairs")
    
    Cumulative_Plot<-arrangeGrob(grid.arrange(Cumulative_Reads_Plot,Cumulative_Base_Pairs_Plot, nrow=2, ncol=1))
    
    ggsave(paste0(Label,"_",label,"_Cumulative_Reads_&_Cumulative_Basepairs.pdf"), device="pdf", Cumulative_Plot, height=10,width=15)
    
    #PLOT PER-HOUR READS/BPs/QUALITY/LENGTH
    
    
    y1<-Reads_Per_Hour
    data1<-data.frame('x'=x,'y'=y1)
    data1$group<-"Reads Per Hour"
    
    Reads_Per_Hour_Plot<-ggplot(data1, aes(x=x, y=y, group=group, fill=group)) +
      geom_line(size=.5) + 
      scale_x_continuous(name="Experimental Time", breaks=(seq(0,Run_Duration,2)))+
      scale_y_continuous(name="Number Of Reads")+
      geom_ribbon(data=subset(data1,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0,show.legend=FALSE) +
      scale_fill_manual(name='', values=c("Reads Per Hour" = "coral3"))+
      theme_bw(base_size = 12)+
      theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=14),legend.background = (element_rect(size=0.5, linetype="solid",colour ="black")), legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=12),axis.title.x = element_text(size=11, face="italic"),axis.title.y = element_text(size=11, face="italic"))+
      theme(legend.position="bottom")+
      ggtitle("Number Of Reads")
    
    
    y2<-Base_Pairs_Per_Hour
    data2<-data.frame('x'=x,'y'=y2)
    data2$group<-"Base Pairs Per Hour"
    
    
    Base_Pairs_Per_Hour_Plot<-ggplot(data2, aes(x=x, y=y, group=group, fill=group)) +
      geom_line(size=.5) + 
      scale_x_continuous(name="Experimental Time", breaks=(seq(0,Run_Duration,2)))+
      scale_y_continuous(name="Number Of Base Pairs")+
      geom_ribbon(data=subset(data2,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0,show.legend=FALSE) +
      scale_fill_manual(name='', values=c("Base Pairs Per Hour" = "darkcyan"))+
      theme_bw(base_size = 12)+
      theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=14),legend.background = (element_rect(size=0.5, linetype="solid",colour ="black")), legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=12),axis.title.x = element_text(size=11, face="italic"),axis.title.y = element_text(size=11, face="italic"))+
      theme(legend.position="bottom")+
      ggtitle("Number Of Base Pairs")
    
    
    y3.0<-log10(Mean_Length_Per_Hour)
    data3.0.0<-data.frame('x'=x,'y'=Mean_Length_Per_Hour)
    data3.0<-data.frame('x'=x,'y'=y3.0)
    data3.0$group<-"Mean Length"
    
    y3.1<-log10(Max_Length_Per_Hour)
    data3.1<-data.frame('x'=x,'y'=y3.1)
    data3.1$group<-"Max Length"
    
    y3.2<-log10(Min_Length_Per_Hour)
    data3.2<-data.frame('x'=x,'y'=y3.2)
    data3.2$group<-"Min Length"
    
    data3<-rbind(data3.0, data3.1, data3.2)
    
    Length_Per_Hour_Plot<-ggplot(data3, aes(x=x, y=y, group=group, fill=group)) +
      geom_line(size=.5) + 
      scale_x_continuous(name="Experimental Time", breaks=(seq(0,Run_Duration,2)))+
      scale_y_continuous(name="Length", breaks=c(2,3,4,5,6), labels=c("100","1000","10000","100000","1000000"))+
      geom_ribbon(data=subset(data3,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0, alpha=.7) +
      scale_fill_manual(name='', values=c("Mean Length" = "darkolivegreen1", "Min Length" = "darkolivegreen", "Max Length" = "cornsilk"))+
      theme_bw(base_size = 12)+
      theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=14),legend.background = (element_rect(size=0.5, linetype="solid",colour ="black")), legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=12),axis.title.x = element_text(size=11, face="italic"),axis.title.y = element_text(size=11, face="italic"))+
      theme(legend.position="bottom")+
      ggtitle("Length")
    
    y4.0<-Mean_Quality_Score_Per_Hour
    data4.0<-data.frame('x'=x,'y'=y4.0)
    data4.0$group<-"Mean Quality"
    
    y4.1<-Min_Quality_Score_Per_Hour
    data4.1<-data.frame('x'=x,'y'=y4.1)
    data4.1$group<-"Min Quality"
    
    y4.2<-Max_Quality_Score_Per_Hour
    data4.2<-data.frame('x'=x,'y'=y4.2)
    data4.2$group<-"Max Quality"
    
    data4<-rbind(data4.0,data4.1,data4.2)
    
    Quality_Score_Per_Hour_Plot<-ggplot(data4, aes(x=x, y=y, group=group, fill=group)) +
      geom_line(size=.5) + 
      scale_x_continuous(name="Experimental Time", breaks=(seq(0,Run_Duration,2)))+
      scale_y_continuous(name="Quality")+
      geom_ribbon(data=subset(data4,x>=0 & x<=Run_Duration),aes(x=x,ymax=y),ymin=0,alpha=.7) +
      scale_fill_manual(name='', values=c("Mean Quality" = "chocolate1", "Min Quality" = "chocolate", "Max Quality" = "cornsilk"))+
      theme_bw(base_size = 12)+
      theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=14),legend.background = (element_rect(size=0.5, linetype="solid",colour ="black")), legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=12),axis.title.x = element_text(size=11, face="italic"),axis.title.y = element_text(size=11, face="italic"))+
      theme(legend.position="bottom")+
      ggtitle("Quality Score")
    
    
    Others_Plot<-arrangeGrob(grid.arrange(Reads_Per_Hour_Plot,Base_Pairs_Per_Hour_Plot,Length_Per_Hour_Plot,Quality_Score_Per_Hour_Plot, nrow=2, ncol=2))
    
    ggsave(paste0(Label,"_",label,"_Reads_Basepairs_Length_Quality.pdf"), device="pdf", Others_Plot, height=10,width=15)
    
    #PASS/FAIL
    
    blank_theme <- theme_minimal()+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_blank()
      )
    
    Data_Pass_Fail_Percentage <- data.frame(
      group = c("Pass", "Fail/Skip"),
      value = c((List.Files.HDF5_Pass.length/(List.Files.HDF5_Pass.length+List.Files.HDF5_Fail_Length+List.Files.HDF5_Skip_Length)), ((List.Files.HDF5_Fail_Length+List.Files.HDF5_Skip_Length)/(List.Files.HDF5_Pass.length+List.Files.HDF5_Fail_Length+List.Files.HDF5_Skip_Length)))
    )
    
    
    Data_Pass_Fail_Percentage_Plot<-ggplot(Data_Pass_Fail_Percentage, aes(x="", y=value, fill=group))+
      geom_bar(width = 1, size = 1, color = "white", stat = "identity")+
      geom_text(aes(label = scales::percent(value)), position = position_stack(vjust = 0.5))+
      coord_polar("y", start=0) +
      scale_fill_brewer(palette="Dark2") +
      blank_theme+
      theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=14),legend.title=element_blank(),legend.position="bottom")+
      ggtitle("Passed and Failed/Skipped Percentage")
    
    
    Data_Pass_Fail_Tot <- data.frame(
      group = c("Pass", "Fail/Skip"),
      value = c(List.Files.HDF5_Pass.length, (List.Files.HDF5_Fail_Length+List.Files.HDF5_Skip_Length))
    )
    
    Data_Pass_Fail_Tot_Plot<-ggplot(Data_Pass_Fail_Tot, aes(x="", y=value, fill=group))+
      geom_bar(width = 1, size = 1, color = "white", stat = "identity")+
      geom_text(aes(label = value), position = position_stack(vjust = 0.5))+
      coord_polar("y", start=0) +
      scale_fill_brewer(palette="Accent") +
      blank_theme+
      theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=14),legend.title=element_blank(), legend.position="bottom")+
      ggtitle("Passed and Failed/Skipped Number")
    
    #WILL BE SAVED WITH GC CONTENT
    
    
    ###################### LENGTH_VS_QUALITY ######################
    
    Tot_Length<-data.frame(Template_Length)
    Tot_Quality<-data.frame(Quality_Score)
    
    
    ScatterTheme <- list(labs(x="Length",y="Quality"),theme_bw(), theme(legend.position=c(1,0),legend.justification=c(1,0), legend.background=element_blank(),legend.direction="horizontal", legend.title=element_text(face="bold.italic"),axis.title=element_text(face="italic")))
    
    
    hist_top_mean_length<-ggplot(Tot_Length, aes(x=Template_Length))+theme_bw()+ geom_histogram(aes(y = ..count../1000),col="darkolivegreen", fill="darkolivegreen1",boundary = min(Tot_Length$Template_Length))+labs(x="",y="Count (x 1e3)")+scale_x_continuous(limits=c(min(Tot_Length$Template_Length),max(Tot_Length$Template_Length)))+theme(axis.title=element_text(face="italic"))
    hist_right_mean_quality<-ggplot(Tot_Quality, aes(x=Quality_Score))+ theme_bw()+ geom_histogram(aes(y = ..count../1000),col="chocolate", fill="chocolate1",boundary = min(Tot_Quality$Quality_Score))+ labs(x="",y="Count (x 1e3)")+coord_flip()+scale_x_continuous(limits=c(min(Tot_Quality$Quality_Score),max(Tot_Quality$Quality_Score)))+theme(axis.title=element_text(face="italic"))
    empty <- ggplot()+geom_point(aes(1,1), colour="white")+
      theme(axis.ticks=element_blank(), 
            panel.background=element_blank(), 
            axis.text.x=element_blank(), axis.text.y=element_blank(),           
            axis.title.x=element_blank(), axis.title.y=element_blank())
    scatter <- ggplot(data.frame(cbind(Tot_Length$Template_Length,Tot_Quality$Quality_Score)),aes(x=Tot_Length$Template_Length, y=Tot_Quality$Quality_Score))+geom_point(col="grey27", size=.09, alpha=.4)+scale_x_continuous(limits=c(min(Tot_Length$Template_Length),max(Tot_Length$Template_Length)))+scale_y_continuous(limits=c(min(Tot_Quality$Quality_Score),max(Tot_Quality$Quality_Score)))+stat_density2d(aes(col=..level.., alpha=..level..)) + scale_color_continuous(low="darkblue",high="darkred") +geom_smooth(method=lm,linetype=2, size=.5,col="black",se=F) + guides(alpha="none",col=guide_legend(title="Density"))+ ScatterTheme
    
    
    
    
    Length_VS_Quality_Plot<-grid.arrange(hist_top_mean_length, empty, scatter, hist_right_mean_quality, ncol=2, nrow=2, widths=c(4,1), heights=c(1, 4))
    
    
    ggsave(paste0(Label,"_",label,"_Length_versus_Quality.pdf"), device="pdf", Length_VS_Quality_Plot, height=10,width=15)
    
    
    # MUX PRODUCTIVITY PER CHANNEL
    
    Mux_Numbers<-c(1:4)
    
    Chan<-as.numeric(Table_HDF5_Reordered[,2])
    Mu<-as.numeric(Table_HDF5_Reordered[,3])
    Le<-as.numeric(Table_HDF5_Reordered[,5])
    
    List_Of_Mux<-list()
    
    for (iii in 1:length(Channels_Number)) {
      Ind_Chann<-which(Chan == Channels_Number[iii])
      Mux_Associated_Number<-sort(Mu[Ind_Chann])
      Table_Mux<-c()
      for (lll in 1:length(Mux_Numbers)) {
        Ind_Mux<-which(Mux_Associated_Number == Mux_Numbers[lll])
        Chan_Mux<-Chan[Ind_Chann][Ind_Mux]
        if (length(Chan_Mux) == 0) {
          Mux<-NA
          Lenght_Per_Mux<-NA
          Table_Mu<-cbind(NA, NA, NA)
          Table_Mux<-rbind(Table_Mux,Table_Mu)
        }
        else {
          Mux<-Mu[Ind_Chann][Ind_Mux]
          Lenght_Per_Mux<-sum(Le[Ind_Chann][Ind_Mux])
          Table_Mu<-cbind(unique(Chan_Mux), unique(Mux),Lenght_Per_Mux)
          Table_Mux<-rbind(Table_Mux,Table_Mu)
        }
      }
      List_Of_Mux[[iii]]<-Table_Mux
    }
    
    Table_Mux_Def<-do.call(rbind,List_Of_Mux)
    
    
    colnames(Table_Mux_Def)<-c("Channel Number", "Mux Number", "Total Reads Produced Per Mux")
    
    #PLOT CORRELATION MATRIXES (CHANNEL AND MUXES)
    
    m1<-matrix(Base_Pairs_Per_Channel[1:32], ncol=8, nrow=4, byrow=TRUE)
    m2<-matrix(Base_Pairs_Per_Channel[449:480], ncol=8, nrow=4, byrow=TRUE)
    m3<-matrix(Base_Pairs_Per_Channel[385:416], ncol=8, nrow=4, byrow=TRUE)
    m4<-matrix(Base_Pairs_Per_Channel[321:352], ncol=8, nrow=4, byrow=TRUE)
    m5<-matrix(Base_Pairs_Per_Channel[257:288], ncol=8, nrow=4, byrow=TRUE)
    m6<-matrix(Base_Pairs_Per_Channel[193:224], ncol=8, nrow=4, byrow=TRUE)
    m7<-matrix(Base_Pairs_Per_Channel[129:160], ncol=8, nrow=4, byrow=TRUE)
    m8<-matrix(Base_Pairs_Per_Channel[65:96], ncol=8, nrow=4, byrow=TRUE)
    mdef3<-rbind(m1,m2,m3,m4,m5,m6,m7,m8)
    m9<-rotate(matrix(Base_Pairs_Per_Channel[33:64], ncol=8, nrow=4, byrow=TRUE))
    m10<-rotate(matrix(Base_Pairs_Per_Channel[481:512], ncol=8, nrow=4, byrow=TRUE))
    m11<-rotate(matrix(Base_Pairs_Per_Channel[417:448], ncol=8, nrow=4, byrow=TRUE))
    m12<-rotate(matrix(Base_Pairs_Per_Channel[353:384], ncol=8, nrow=4, byrow=TRUE))
    m13<-rotate(matrix(Base_Pairs_Per_Channel[289:320], ncol=8, nrow=4, byrow=TRUE))
    m14<-rotate(matrix(Base_Pairs_Per_Channel[225:256], ncol=8, nrow=4, byrow=TRUE))
    m15<-rotate(matrix(Base_Pairs_Per_Channel[161:192], ncol=8, nrow=4, byrow=TRUE))
    m16<-rotate(matrix(Base_Pairs_Per_Channel[97:128], ncol=8, nrow=4, byrow=TRUE))
    mdef4<-rbind(m9,m10,m11,m12,m13,m14,m15,m16)
    Matrixbpchannel<-cbind(mdef3,mdef4)
    
    BasePairs_Per_Mux<-as.numeric(Table_Mux_Def[,3])
    
    First_Eight_Disposition<-c(3,4,1,2,6,5,8,7)
    Second_Eight_Disposition<-Increment(First_Eight_Disposition)
    Third_Eight_Disposition<-Increment(Second_Eight_Disposition)
    Fouth_Eight_Disposition<-Increment(Third_Eight_Disposition)
    First_Line<-c(First_Eight_Disposition,Second_Eight_Disposition,Third_Eight_Disposition,Fouth_Eight_Disposition)
    Second_Line<-Increment2(First_Line)
    Third_Line<-Increment2(Second_Line)
    Fourth_Line<-Increment2(Third_Line)
    First_Block<-c(First_Line,Second_Line,Third_Line,Fourth_Line)
    Second_Block<-Increment3(First_Block)
    Third_Block<-Increment3(Second_Block)
    Fourth_Block<-Increment3(Third_Block)
    Fifth_Block<-Increment3(Fourth_Block)
    Sixth_Block<-Increment3(Fifth_Block)
    Seventh_Block<-Increment3(Sixth_Block)
    Eight_Block<-Increment3(Seventh_Block)
    Ninth_Block<-Increment3(Eight_Block)
    Tenth_Block<-Increment3(Ninth_Block)
    Eleventh_Block<-Increment3(Tenth_Block)
    Twelfth_Block<-Increment3(Eleventh_Block)
    Thirtheenth_Block<-Increment3(Twelfth_Block)
    Fourtheenth_Block<-Increment3(Thirtheenth_Block)
    Fiftheenth_Block<-Increment3(Fourtheenth_Block)
    Sixtheenth_Block<-Increment3(Fiftheenth_Block)
    
    M1<-matrix(BasePairs_Per_Mux[First_Block], ncol= 32, nrow= 4, byrow=TRUE)
    M2<-matrix(BasePairs_Per_Mux[Fiftheenth_Block], ncol= 32, nrow= 4, byrow=TRUE)
    M3<-matrix(BasePairs_Per_Mux[Thirtheenth_Block], ncol= 32, nrow= 4, byrow=TRUE)
    M4<-matrix(BasePairs_Per_Mux[Eleventh_Block], ncol= 32, nrow= 4, byrow=TRUE)
    M5<-matrix(BasePairs_Per_Mux[Ninth_Block], ncol= 32, nrow= 4, byrow=TRUE)
    M6<-matrix(BasePairs_Per_Mux[Seventh_Block], ncol= 32, nrow= 4, byrow=TRUE)
    M7<-matrix(BasePairs_Per_Mux[Fifth_Block], ncol= 32, nrow= 4, byrow=TRUE)
    M8<-matrix(BasePairs_Per_Mux[Third_Block], ncol= 32, nrow= 4, byrow=TRUE)
    Mdef3<-rbind(M1,M2,M3,M4,M5,M6,M7,M8)
    M9<-rotate(matrix(BasePairs_Per_Mux[Second_Block], ncol= 32, nrow= 4, byrow=TRUE))
    M10<-rotate(matrix(BasePairs_Per_Mux[Sixtheenth_Block],ncol= 32, nrow= 4, byrow=TRUE))
    M11<-rotate(matrix(BasePairs_Per_Mux[Fourtheenth_Block], ncol= 32, nrow= 4, byrow=TRUE))
    M12<-rotate(matrix(BasePairs_Per_Mux[Twelfth_Block], ncol= 32, nrow= 4, byrow=TRUE))
    M13<-rotate(matrix(BasePairs_Per_Mux[Tenth_Block], ncol= 32, nrow= 4, byrow=TRUE))
    M14<-rotate(matrix(BasePairs_Per_Mux[Eight_Block], ncol= 32, nrow= 4, byrow=TRUE))
    M15<-rotate(matrix(BasePairs_Per_Mux[Sixth_Block], ncol= 32, nrow= 4, byrow=TRUE))
    M16<-rotate(matrix(BasePairs_Per_Mux[Fourth_Block], ncol= 32, nrow= 4, byrow=TRUE))
    Mdef4<-rbind(M9,M10,M11,M12,M13,M14,M15,M16)
    MatrixMuxActivity<-cbind(Mdef3,Mdef4)
    
    
    #PLOTTING "FALSE" CORRELATION MATRIXES
    
    Palette <- colorRampPalette(brewer.pal(9, "Reds"))
    
    adjMatrixbpchannel<-melt(rotate(t(Matrixbpchannel)))
    Plot_Channel_Activity<-ggplot(data=adjMatrixbpchannel, aes(x=Var1, y=Var2)) +
      geom_tile(aes(fill=value), color="white", size=2)+
      scale_fill_gradientn(colours= Palette(4), na.value="grey70",limits=c(min(adjMatrixbpchannel[,3], na.rm=TRUE), max(adjMatrixbpchannel[,3], na.rm=TRUE))) +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="bottom",
            legend.title=element_text(face="bold.italic"),
            panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),plot.background=element_blank(),legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5,face="bold.italic",size=16),legend.text.align = 0)+
      guides(fill = guide_colorbar(barwidth= 61,barheight=.5,title="Base Pairs Productivity", title.position="top",title.hjust=0.5))+
      ggtitle("Channels Activity")
    
    adjMatrixMuxActivity<-melt(rotate(t(MatrixMuxActivity)))
    Plot_Mux_Activity<-ggplot(data=adjMatrixMuxActivity, aes(x=Var1, y=Var2)) +
      geom_tile(aes(fill=value), color="white", size=2)+
      scale_fill_gradientn(colours=Palette(4), na.value="grey70",limits=c(min(adjMatrixMuxActivity[,3], na.rm=TRUE), max(adjMatrixMuxActivity[,3], na.rm=TRUE))) +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="bottom",
            legend.title=element_text(face="bold.italic"),
            panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),plot.background=element_blank(),legend.text=element_text(size=10),
            plot.title = element_text(hjust = 0.5,face="bold.italic",size=16))+
      guides(fill = guide_colorbar(barwidth= 61, barheight=.5,title="Base Pairs Productivity", title.position="top",title.hjust=0.5))+
      ggtitle("Muxes Activity")
    
    
    Plot_Tot<-grid.arrange(arrangeGrob(Plot_Channel_Activity,Plot_Mux_Activity,nrow=2, ncol=1, widths=15, heights=c(12,12)))
    
    ggsave(paste0(Label,"_",label,"_Channels_and_Muxes_Activity.pdf"), device="pdf",Plot_Tot, height=10, width=18)
    
    
    ###PLOT GC CONTENT#######
    
    if (NanoGTable[1,7] == "GC_Content") {
      
      pdf(paste0(Label,"_",label,"_Pass_Fail_Skip_NO_GC_Content.pdf"), height=10, width=15, onefile=TRUE)
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(1, 2)))
      print(Data_Pass_Fail_Percentage_Plot, vp=define_region(1, 1))
      print(Data_Pass_Fail_Tot_Plot, vp = define_region(1, 2))
      dev.off()
    }
    
    else {
      
      GC_Content_To_Plot<-as.numeric(NanoTable[,7])
      Hist_GC_Content<-ggplot(data.frame(GC_Content_To_Plot), aes(x=GC_Content_To_Plot))+theme_bw()+ geom_histogram(col="midnightblue", fill="cyan4")+labs(x="",y="Count")+ggtitle("GC Content")+theme(plot.title = element_text(hjust = 0.5,face="bold.italic",size=16), axis.title.y=element_text(face="italic"))
      pdf(paste0(Label,"_",label,"_Pass_Fail_Skip_and_GC_Content.pdf"), height=10, width=15, onefile=TRUE)
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(2, 2)))
      print(Data_Pass_Fail_Percentage_Plot, vp=define_region(1, 1))
      print(Data_Pass_Fail_Tot_Plot, vp = define_region(1, 2))
      print(Hist_GC_Content, vp = define_region(2, 1:2))
      dev.off()
    }
    
    #save data for other analyses
    
    write.table(data1, file.path(Directory, paste0(Label,"_",label, "_ReadsProduced.txt")),col.names=T, sep="\t")
    write.table(data2, file.path(Directory, paste0(Label,"_",label, "_BasePairsProduced.txt")),col.names=T, sep="\t")
    write.table(data3.0.0, file.path(Directory, paste0(Label,"_",label, "_MeanLength.txt")),col.names=T, sep="\t")
    write.table(data4.0, file.path(Directory, paste0(Label,"_",label, "_MeanQuality.txt")),col.names=T, sep="\t")
    
    message("Done!")
    
  } 
}

