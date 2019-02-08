############################################ NanoCompare ############################################


#' @title Compare ONT experiments 
#' @description NanoCompare plots comparison statistics for MinION and GridION X5 experiments analyzed with the other functions from this package.
#' @param DataIn  Character vector containing paths to folders contaning analyzed MinION and/or GridION X5 experiments
#' @param DataOut Where NanoCompare will save the results
#' @param Labels Character vector containing ordered labels used to identify experiments in DataIn
#' @return Plots: \cr 
#' - Violins plot comparing reads number, base pairs number, reads length and reads quality calculated every 10 hours of experimental run (Violins.pdf); \cr
#' - Histograms plot comparing Length, Quality and GCC (optionally) count distributions (Histograms.pdf);
#' @examples
#' #do not run
#' DataIn<-c("Path/To/AnalyzedFolder1","Path/To/AnalyzedFolder2","Path/To/AnalyzedFolder3",...) #path to the NanoR-analyzed data
#' Labels<-c("Label1","Label2","Label3") #labels used
#' NanoCompare(DataIn=DataIn,DataOut="Path/To/DataOut",Labels=Labels,GCC=TRUE) #compare




NanoCompare<-function(DataIn,DataOut,Labels,GCC=TRUE) {
  
  library(scales)
  library(ggplot2)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
  
  Directory<-file.path(DataOut)
  dir.create(Directory, showWarnings = FALSE)
  setwd(Directory)
  
  g_legend<-function(a.gplot)
  {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  } 
  
  split.vars<-DataIn
  split.vars2<-Labels
  
  if (length(split.vars) != length(split.vars2)) {
    stop("There is no correspondence in number of folders and number of labels. For each experiment one label must be provided!")
  }
  
  Name_Of_Runs_To_Compare<-split.vars2
  
  Table_List<-list()
  Length_Data<-c()
  
  for (n in 1:length(split.vars)) {
    Data_Reads<-read.table(file.path(split.vars[n], paste0(Name_Of_Runs_To_Compare[n],"_ReadsProduced.txt")), header=T, sep="\t")
    Length_Data[n]<-length(Data_Reads[,1])
    Table_List[[n]]<-Data_Reads
  }
  
  Major_Index<-which.max(Length_Data)
  
  if (Length_Data[Major_Index] %% 10 != 0) {
    Group2<-c(rep(head(seq(min(Table_List[[Major_Index]]$x), max(Table_List[[Major_Index]]$x), 10), length(seq(min(Table_List[[Major_Index]]$x), max(Table_List[[Major_Index]]$x), 10))-1), each=20),rep(tail(seq(min(Table_List[[Major_Index]]$x), max(Table_List[[Major_Index]]$x), 10),1),2*(max(Table_List[[Major_Index]]$x)-tail(seq(min(Table_List[[Major_Index]]$x), max(Table_List[[Major_Index]]$x), 10),1))+1))
  }
  
  if (Length_Data[Major_Index] %% 10 == 0) {
    Group2<-c(rep(head(seq(min(Table_List[[Major_Index]]$x), max(Table_List[[Major_Index]]$x), 10), length(seq(min(Table_List[[Major_Index]]$x), max(Table_List[[Major_Index]]$x), 10))-1), each=20),max(Table_List[[Major_Index]]$x)) 
  }
  
  
  ##########
  New_Group_2<-c()
  
  for (i in 1:length(Group2)) {
    if (Group2[i] == 0) {
      New_Group_2[i]<-as.character("0-10 hours")
    }
    if (Group2[i] == 10) {
      New_Group_2[i]<-as.character("10-20 hours")
    }
    if (Group2[i] == 20) {
      New_Group_2[i]<-as.character("20-30 hours")
    }
    if (Group2[i] == 30) {
      New_Group_2[i]<-as.character("30-40 hours")
    }
    if (Group2[i] == 40) {
      New_Group_2[i]<-as.character("40-50 hours")
    }
    if (Group2[i] == 50) {
      New_Group_2[i]<-as.character("50-60 hours")
    }
  }
  
  ############VIOLINS COMPARISON#########
  
  List_Data_Reads<-list()
  List_Data_BasePairs<-list()
  List_Data_Length<-list()
  List_Data_Quality<-list()
  
  
  options(scipen=9999)


  message("Creating comparison violins...")
  
  for (nn in 1:length(split.vars)) {
    Data_Reads1<-Table_List[[nn]]
    if (length(Data_Reads1[,1]) != length(Group2)) {  
      x<-c(Data_Reads1$x,seq(max(Data_Reads1$x)+.5, max(Table_List[[Major_Index]]$x),0.5))
      y<-c(Data_Reads1$y, rep(0, length(seq(max(Data_Reads1$x)+.5, max(Table_List[[Major_Index]]$x),0.5))))
      group<-rep("Reads", length(y))
      value<-rep(Name_Of_Runs_To_Compare[nn], length(y))
      group2<-factor(as.character(New_Group_2), levels=c("0-10 hours","10-20 hours","20-30 hours","30-40 hours","40-50 hours"))
      Data_Reads1<-data.frame(x,y,group,value,group2)
    }
    else {    
      Data_Reads1$group<-"Reads"
      Data_Reads1$value<-Name_Of_Runs_To_Compare[nn]
      Data_Reads1$group2<-factor(as.character(New_Group_2), levels=c("0-10 hours","10-20 hours","20-30 hours","30-40 hours","40-50 hours"))
    }
    List_Data_Reads[[nn]]<-Data_Reads1
    if (length(List_Data_Reads) == length(split.vars)) {
      Reads_To_Plot<-do.call(rbind, List_Data_Reads)
      Violins_Reads<-ggplot(Reads_To_Plot, aes(x=group,y=y, fill=value)) + 
        geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
        labs(x="", y = " # Reads")+ 
        scale_fill_brewer(palette= "Set3") + 
        scale_x_discrete(breaks=NULL)+
        theme_minimal(base_size=12)+
        theme(panel.grid.minor=element_blank(),legend.position="none",axis.text.x=element_blank(),axis.title.y=element_text(size=10, face="bold.italic"),axis.text.y=element_text(size=8), strip.text.x=element_text(face="bold.italic"))+
        guides(alpha = FALSE)+
        facet_wrap(~ group2, scale="free_y", nrow=1)
    }
    Data_BasePairs1<-read.table(file.path(split.vars[nn], paste0(Name_Of_Runs_To_Compare[nn],"_BasePairsProduced.txt")), header=T, sep="\t")
    if (length(Data_BasePairs1[,1]) != length(Group2)) {  
      x<-c(Data_BasePairs1$x,seq(max(Data_BasePairs1$x)+.5, max(Table_List[[Major_Index]]$x),0.5))
      y<-c(Data_BasePairs1$y, rep(0, length(seq(max(Data_BasePairs1$x)+.5, max(Table_List[[Major_Index]]$x),0.5))))
      group<-rep("Base Pairs", length(y))
      value<-rep(Name_Of_Runs_To_Compare[nn], length(y))
      group2<-factor(as.character(New_Group_2), levels=c("0-10 hours","10-20 hours","20-30 hours","30-40 hours","40-50 hours"))
      Data_BasePairs1<-data.frame(x,y,group,value,group2)
    }
    else {
      Data_BasePairs1$group<-"Base Pairs"
      Data_BasePairs1$value<-Name_Of_Runs_To_Compare[nn]
      Data_BasePairs1$group2<-factor(as.character(New_Group_2), levels=c("0-10 hours","10-20 hours","20-30 hours","30-40 hours","40-50 hours"))
    }
    List_Data_BasePairs[[nn]]<-Data_BasePairs1
    if (length(List_Data_BasePairs) == length(split.vars)) {
      BasePairs_To_Plot<-do.call(rbind, List_Data_BasePairs)
      Violins_BasePairs<-ggplot(BasePairs_To_Plot, aes(x=group,y=y, fill=value)) + 
        geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
        labs(x="", y = "# Base Pairs")+ 
        scale_fill_brewer(palette="Set3") + 
        scale_x_discrete(breaks=NULL)+
        theme_minimal(base_size=12)+
        theme(panel.grid.minor=element_blank(),legend.position="none",axis.text.x=element_blank(),axis.title.y=element_text(size=10, face="bold.italic"),axis.text.y=element_text(size=8),strip.text=element_blank())+
        guides(alpha = FALSE)+
        facet_wrap(~ group2, scale="free_y", nrow=1)
    }
    Data_Mean_Length1<-read.table(file.path(split.vars[nn], paste0(Name_Of_Runs_To_Compare[nn],"_MeanLength.txt")), header=T, sep="\t")
    if (length(Data_Mean_Length1[,1]) != length(Group2)) {
      x<-c(Data_Mean_Length1$x,seq(max(Data_Mean_Length1$x)+.5, max(Table_List[[Major_Index]]$x),0.5))
      y<-c(Data_Mean_Length1$y, rep(0, length(seq(max(Data_Mean_Length1$x)+.5, max(Table_List[[Major_Index]]$x),0.5))))
      group<-rep("Length", length(y))
      value<-rep(Name_Of_Runs_To_Compare[nn], length(y))
      group2<-factor(as.character(New_Group_2), levels=c("0-10 hours","10-20 hours","20-30 hours","30-40 hours","40-50 hours"))
      Data_Mean_Length1<-data.frame(x,y,group,value,group2)
    }
    else {
      Data_Mean_Length1$group<-"Length"
      Data_Mean_Length1$value<-Name_Of_Runs_To_Compare[nn]
      Data_Mean_Length1$group2<-factor(as.character(New_Group_2), levels=c("0-10 hours","10-20 hours","20-30 hours","30-40 hours","40-50 hours"))
    }
    List_Data_Length[[nn]]<-Data_Mean_Length1
    if (length(List_Data_Length) == length(split.vars)) {
      Length_To_Plot<-do.call(rbind, List_Data_Length)
      Violins_Length<-ggplot(Length_To_Plot, aes(x=group,y=y, fill=value)) + 
        geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
        labs(x="", y = "Length (bp)")+ 
        scale_fill_brewer(palette="Set3") + 
        scale_x_discrete(breaks=NULL)+
        theme_minimal(base_size=12)+
        theme(panel.grid.minor=element_blank(),legend.position="none",axis.text.x=element_blank(),axis.title.y=element_text(size=10, face="bold.italic"),axis.text.y=element_text(size=8), strip.text=element_blank())+
        guides(alpha = FALSE)+
        facet_wrap(~ group2, scale="free_y", nrow=1)
    }
    Data_Mean_Quality1<-read.table(file.path(split.vars[nn], paste0(Name_Of_Runs_To_Compare[nn],"_MeanQuality.txt")), header=T, sep="\t")
    if (length(Data_Mean_Quality1[,1]) != length(Group2)) {
      x<-c(Data_Mean_Quality1$x,seq(max(Data_Mean_Quality1$x)+.5, max(Table_List[[Major_Index]]$x),0.5))
      y<-c(Data_Mean_Quality1$y, rep(0, length(seq(max(Data_Mean_Quality1$x)+.5, max(Table_List[[Major_Index]]$x),0.5))))
      group<-rep("Quality", length(y))
      value<-rep(Name_Of_Runs_To_Compare[nn], length(y))
      group2<-factor(as.character(New_Group_2), levels=c("0-10 hours","10-20 hours","20-30 hours","30-40 hours","40-50 hours"))
      Data_Mean_Quality1<-data.frame(x,y,group,value,group2)
    }
    else {
      Data_Mean_Quality1$group<-"Quality"
      Data_Mean_Quality1$value<-Name_Of_Runs_To_Compare[nn]
      Data_Mean_Quality1$group2<-factor(as.character(New_Group_2), levels=c("0-10 hours","10-20 hours","20-30 hours","30-40 hours","40-50 hours"))
    }
    List_Data_Quality[[nn]]<-Data_Mean_Quality1
    if (length(List_Data_Quality) == length(split.vars)) {
      Quality_To_Plot<-do.call(rbind, List_Data_Quality)
      Violins_Quality<-ggplot(Quality_To_Plot, aes(x=group,y=y, fill=value)) + 
        geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
        labs(x="", y = "Quality (Phred)")+ 
        scale_fill_brewer(palette="Set3") + 
        scale_x_discrete(breaks=NULL)+
        theme_minimal(base_size=12)+
        theme(panel.grid.minor=element_blank(),legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=8),axis.text.x=element_blank(),axis.title.y=element_text(size=10, face="bold.italic"),axis.text.y=element_text(size=8), strip.text=element_blank())+
        guides(alpha = FALSE,fill = guide_legend(nrow = 1,label.position = "bottom"))+
        facet_wrap(~ group2, scale="free_y", nrow=1)
    }
  }
  mylegend<-g_legend(Violins_Quality)
  
  
  All_Violins_Together<-grid.arrange(Violins_Reads + theme(legend.position="none"), Violins_BasePairs + theme(legend.position="none"), Violins_Length + theme(legend.position="none"), Violins_Quality + theme(legend.position="none"), mylegend, nrow=5,heights=c(2.3,2.3,2.3,2.3,.8))


  
  ggsave("Comparison_Violins.pdf", device="pdf",All_Violins_Together, height=10, width=15)

  message("Done!")
  
  ######HISTOGRAM COMPARISON
  
  List_Reads_Length<-list()
  List_Reads_Quality<-list()
  List_Reads_GC<-list()
  
  if (GCC == TRUE) {

  	message("Creating comparison histograms with GC content...")
    
    
    for (nnn in 1:length(split.vars)) {
      Info_Table<-read.table(file.path(split.vars[nnn], paste0(Name_Of_Runs_To_Compare[nnn],"_Information_Table.txt")), header=T, sep="\t")
      Read_Length<-as.numeric(log10(Info_Table[,5]))
      value<-rep(Name_Of_Runs_To_Compare[nnn], length(Read_Length))
      group<-rep("Reads Length",length(Read_Length))
      Read_Length_DF<-data.frame(x=Read_Length, value, group)
      List_Reads_Length[[nnn]]<-Read_Length_DF
      if (length(List_Reads_Length) == length(split.vars)) {
        Table_Reads<-do.call(rbind,List_Reads_Length)
        Reads_Plot<-ggplot(Table_Reads, aes(x=x, fill=value)) + theme_bw() +geom_histogram(col="black")+scale_fill_brewer(palette="Set3") + labs(x="Length (bp)",y="Count")+theme(legend.position="none",axis.title.y=element_text(face="bold.italic"),axis.title.x=element_text(face="bold.italic"),strip.text.x = element_blank())+facet_wrap(~ value, scale="free_y",nrow=1)
      }
      Read_Quality<-as.numeric(Info_Table[,6])
      value<-rep(Name_Of_Runs_To_Compare[nnn], length(Read_Quality))
      group<-rep("Reads Quality",length(Read_Quality))
      Read_Quality_DF<-data.frame(x=Read_Quality, value, group)
      List_Reads_Quality[[nnn]]<-Read_Quality_DF
      if (length(List_Reads_Quality) == length(split.vars)) {
        Table_Quality<-do.call(rbind,List_Reads_Quality)
        Quality_Plot<-ggplot(Table_Quality, aes(x=x, fill=value)) + theme_bw()+ xlim(0,20)+geom_histogram(col="black")+scale_fill_brewer(palette="Set3") + labs(x="Quality (Phred)",y="Count")+theme(legend.position="none",axis.title.y=element_text(face="bold.italic"),axis.title.x=element_text(face="bold.italic"),strip.text.x = element_blank())+facet_wrap(~value, scale="free_y",nrow=1)
      }
      Read_GC_Content<-as.numeric(Info_Table[,7])
      value<-rep(Name_Of_Runs_To_Compare[nnn], length(Read_GC_Content))
      group<-rep("Reads GC Content",length(Read_GC_Content))
      Read_GC_DF<-data.frame(x=Read_GC_Content, value, group)
      List_Reads_GC[[nnn]]<-Read_GC_DF
      if(length(List_Reads_GC) == length(split.vars)) {
        GC_Table<-do.call(rbind,c(List_Reads_GC))
        GC_Plot<-ggplot(GC_Table, aes(x=x, fill=value)) + theme_bw()+ xlim(0,1)+geom_histogram(col="black")+scale_fill_brewer(palette="Set3") + labs(x="GCC",y="Count")+theme(legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=8),axis.title.y=element_text(face="bold.italic"),axis.title.x=element_text(face="bold.italic"),strip.text.x = element_blank())+facet_wrap(~value, scale="free_y",nrow=1)+guides(fill = guide_legend(nrow = 1,label.position = "bottom"))
      }
    }
    
    mylegend2<-g_legend(GC_Plot)
    
    Reads_Length_Quality_GC<-grid.arrange(Reads_Plot,Quality_Plot,GC_Plot+theme(legend.position="none"),mylegend2, nrow=4, heights=c(2.5,2.5,2.5,.8))


    ggsave("Comparison_Histograms.pdf",device="pdf",Reads_Length_Quality_GC,height=10, width=15)

    message("Done!")
  }
  
  if (GCC == FALSE) {

  	message("Creating comparison histograms without GC content...")
    
    for (nnn in 1:length(split.vars)) {
      Info_Table<-read.table(file.path(split.vars[nnn], paste0(Name_Of_Runs_To_Compare[nnn],"_Information_Table.txt")), header=T, sep="\t")
      Read_Length<-as.numeric(log10(Info_Table[,5]))
      value<-rep(Name_Of_Runs_To_Compare[nnn], length(Read_Length))
      group<-rep("Reads Length",length(Read_Length))
      Read_Length_DF<-data.frame(x=Read_Length, value, group)
      List_Reads_Length[[nnn]]<-Read_Length_DF
      if (length(List_Reads_Length) == length(split.vars)) {
        Table_Reads<-do.call(rbind,List_Reads_Length)
        Reads_Plot<-ggplot(Table_Reads, aes(x=x, fill=value)) + theme_bw() +geom_histogram(col="black")+scale_fill_brewer(palette="Set3") + labs(x="Length (bp)",y="Count")+theme(legend.position="none",axis.title.y=element_text(face="bold.italic"),axis.title.x=element_text(face="bold.italic"),strip.text.x = element_blank())+facet_wrap(~ value, scale="free_y",nrow=1)
      }
      Read_Quality<-as.numeric(Info_Table[,6])
      value<-rep(Name_Of_Runs_To_Compare[nnn], length(Read_Quality))
      group<-rep("Reads Quality",length(Read_Quality))
      Read_Quality_DF<-data.frame(x=Read_Quality, value, group)
      List_Reads_Quality[[nnn]]<-Read_Quality_DF
      if (length(List_Reads_Quality) == length(split.vars)) {
        Table_Quality<-do.call(rbind,List_Reads_Quality)
        Quality_Plot<-ggplot(Table_Quality, aes(x=x, fill=value)) + theme_bw()+ xlim(0,20)+geom_histogram(col="black")+scale_fill_brewer(palette="Set3") + labs(x="Quality (Phred)",y="Count")+theme(legend.text=element_text(size=8),legend.position="bottom",legend.title=element_blank(),axis.title.y=element_text(face="bold.italic"),axis.title.x=element_text(face="bold.italic"),strip.text.x = element_blank())+facet_wrap(~value, scale="free_y",nrow=1)+guides(fill = guide_legend(nrow = 1,label.position = "bottom"))
      }
    }
    
    mylegend2<-g_legend(Quality_Plot)
    
    Reads_Length_Quality<-grid.arrange(Reads_Plot,Quality_Plot+theme(legend.position="none"),mylegend2, nrow=3,heights=c(3.2,3.2,0.6))

    ggsave("Comparison_Histograms.pdf",device="pdf",Reads_Length_Quality,height=10, width=15)

    message("Done!")
  }
}
