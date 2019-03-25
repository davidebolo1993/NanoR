############################################ NanoCompare ############################################


#' @title Compare ONT experiments 
#' @description NanoCompare fast build violin plots, comparing MinION and GridION X5 experiments analyzed with the other functions from this package. Bins of 30 minutes are taken.
#' @param DataIn  Character vector with paths to "DataForComparison" for each experiment
#' @param DataOut Where results will be saved.
#' @param Labels Character vector containing ordered labels to name experiments in "DataIn"
#' @return Plot: \cr 
#' - Violins.pdf; \cr
#' @examples
#' #do not run
#' DataIn<-c("Path/To/AnalyzedFolder1/DataForComparison","Path/To/AnalyzedFolder2/DataForComparison",...)
#' Labels<-c("Label1","Label2") #labels
#' NanoCompare(DataIn=DataIn,DataOut="Path/To/DataOut",Labels=Labels) #compare


NanoCompare<-function(DataIn,DataOut,Labels) { #Now is a lot faster and can deal with theoretically any number of experiments
  
  library(ggplot2)
  library(RColorBrewer)
  library(grid)


  define_region <- function(row, col)
  {
    viewport(layout.pos.row = row, layout.pos.col = col)
  }

  if (length(DataIn) != length(Labels)) { #check
    
    stop("Different number of input data and labels")
  
  }


  Directory<-file.path(DataOut,"ComparisonPlots")
  dir.create(Directory, showWarnings = FALSE)


  #get data inside the DataForComparison folder.

  Reads_<-list()
  BasePairs_<-list()
  Length_<-list()
  Quality_<-list()

  for (i in 1:(length(DataIn))) {
    
    IsCorrect<-list.files(DataIn[i])

    if (length(IsCorrect) != 4) { #check

      stop("Expected different number of files in ", DataIn[i])

    }

    Reads<-read.table(file.path(DataIn[i], 'Reads.txt'), sep='\t', header=TRUE)
    Reads$sample <- as.character(Labels[i])
    
    Times<-c()

    for (l in 1:length(Reads$x)) {

      if (as.numeric(Reads$x[l]) >=0 && as.numeric(Reads$x[l]) < 10) {

        Times[l]<-'0-10 hrs'

      }

      if (as.numeric(Reads$x[l]) >= 10 && as.numeric(Reads$x[l]) < 20) {

        Times[l]<-'10-20 hrs'

      }

      if (as.numeric(Reads$x[l]) >= 20 && as.numeric(Reads$x[l]) < 30) {

        Times[l]<-'20-30 hrs'

      }

      if (as.numeric(Reads$x[l]) >= 30 && as.numeric(Reads$x[l]) < 40) {

        Times[l]<-'30-40 hrs'

      }

      if (as.numeric(Reads$x[l]) >= 40 && as.numeric(Reads$x[l]) < 50) {

        Times[l]<-'40-50 hrs'

      }

      if (as.numeric(Reads$x[l]) >= 50  && as.numeric(Reads$x[l]) < 60) {

        Times[l]<-'50-60 hrs'

      }

      if (as.numeric(Reads$x[l]) >= 60  && as.numeric(Reads$x[l]) < 70) {

        Times[l]<-'60-70 hrs'

      }

      if (as.numeric(Reads$x[l]) >= 70  && as.numeric(Reads$x[l]) <= 80) {

        Times[l]<-'70-80 hrs'

      }



    }

    Reads$times<-Times
    Reads_[[i]]<-Reads
    BasePairs<-read.table(file.path(DataIn[i], 'Bases.txt'), sep='\t', header=TRUE)
    BasePairs$sample <- as.character(Labels[i])
    BasePairs$times<-Times

    BasePairs_[[i]]<-BasePairs    
    Length<-read.table(file.path(DataIn[i], 'Length.txt'), sep='\t', header=TRUE)
    Length$sample <- as.character(Labels[i])
    Length$times<-Times
    Length_[[i]]<-Length
    Quality<-read.table(file.path(DataIn[i], 'Quality.txt'), sep='\t', header=TRUE)
    Quality$sample <- as.character(Labels[i])
    Quality$times<-Times
    Quality_[[i]]<-Quality

  
  }

  TotReads<-do.call(rbind,Reads_)
  TotBases<-do.call(rbind,BasePairs_)
  TotLength<-do.call(rbind,Length_)
  TotQuality<-do.call(rbind,Quality_)

  p_reads <- ggplot(TotReads, aes(x=sample, y=log10(y), fill=sample)) + #use log instead of number: so that there is not too much difference between experiments
  geom_violin(trim=FALSE,draw_quantiles = c(0.25, 0.5, 0.75), show.legend=FALSE)+
  #geom_boxplot(width=0.1, fill="white", show.legend=FALSE)+
  stat_summary(fun.data=mean_sdl,fun.args = list(mult = 1), geom="pointrange", color="black",show.legend=FALSE)+
  geom_jitter(shape=16, position=position_jitter(.05), size=.5, fill="grey30", show.legend=FALSE)+ ## check
  scale_fill_brewer(palette="RdBu")+
  labs(x="",y=expression("# reads "["(log10)"]))+
  theme_minimal()+
  theme(strip.background = element_blank(),panel.border = element_rect(colour = "grey40",fill=NA), axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~ times, scale="free_y", nrow=1)


  p_bases <- ggplot(TotBases, aes(x=sample, y=log10(y), fill=sample)) + 
  geom_violin(trim=FALSE,draw_quantiles = c(0.25, 0.5, 0.75), show.legend=FALSE)+
  #geom_boxplot(width=0.1, fill="white", show.legend=FALSE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange", color="black",show.legend=FALSE)+
  geom_jitter(shape=16, position=position_jitter(.05), size=.5, color="black", show.legend=FALSE)+ ## check
  scale_fill_brewer(palette="RdBu")+
  labs(x="",y=expression("# bps "["(log10)"]))+  
  theme_minimal()+
  theme(strip.text=element_blank(),panel.border = element_rect(colour = "grey40",fill=NA), axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~ times, scale="free_y", nrow=1)


  p_bases_single <- ggplot(TotBases, aes(x=sample, y=log10(y), fill=sample)) + 
  geom_violin(trim=FALSE,draw_quantiles = c(0.25, 0.5, 0.75), show.legend=FALSE)+
  #geom_boxplot(width=0.1, fill="white", show.legend=FALSE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange", color="black",show.legend=FALSE)+
  geom_jitter(shape=16, position=position_jitter(.05), size=.5, color="black", show.legend=FALSE)+ ## check
  scale_fill_brewer(palette="RdBu")+
  labs(x="",y=expression("# bps "["(log10)"]))+  
  theme_minimal()+
  theme(strip.background = element_blank(),panel.border = element_rect(colour = "grey40",fill=NA), axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~ times, scale="free_y", nrow=1)


  p_length <- ggplot(TotLength, aes(x=sample, y=y, fill=sample)) + 
  geom_violin(trim=FALSE,draw_quantiles = c(0.25, 0.5, 0.75), show.legend=FALSE)+
  #geom_boxplot(width=0.1, fill="white", show.legend=FALSE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange", color="black",show.legend=FALSE)+
  geom_jitter(shape=16, position=position_jitter(.05), size=.5, fill="grey30", show.legend=FALSE)+ ## check
  scale_fill_brewer(palette="RdBu")+
  labs(x="",y=expression("length "["(bps)"]))+  
  theme_minimal()+
  theme(strip.text=element_blank(),panel.border = element_rect(colour = "grey40",fill=NA), axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~ times, scale="free_y", nrow=1)

  p_length_single <- ggplot(TotLength, aes(x=sample, y=y, fill=sample)) + 
  geom_violin(trim=FALSE,draw_quantiles = c(0.25, 0.5, 0.75), show.legend=FALSE)+
  #geom_boxplot(width=0.1, fill="white", show.legend=FALSE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange", color="black",show.legend=FALSE)+
  geom_jitter(shape=16, position=position_jitter(.05), size=.5, fill="grey30", show.legend=FALSE)+ ## check
  scale_fill_brewer(palette="RdBu")+
  labs(x="",y=expression("length "["(bps)"]))+  
  theme_minimal()+
  theme(strip.background = element_blank(),panel.border = element_rect(colour = "grey40",fill=NA), axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~ times, scale="free_y", nrow=1)
   

   
  p_quality <- ggplot(TotQuality, aes(x=sample, y=y, fill=sample)) + 
  geom_violin(trim=FALSE,draw_quantiles = c(0.25, 0.5, 0.75), show.legend=FALSE)+
  #geom_boxplot(width=0.1, fill="white", show.legend=FALSE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange", color="black",show.legend=FALSE)+
  geom_jitter(shape=16, position=position_jitter(.05), size=.5, fill="grey30", show.legend=FALSE)+ ## check
  scale_fill_brewer(palette="RdBu")+
  labs(x="",y=expression("quality "["(phred)"]))+
  theme_minimal()+
  theme(strip.text=element_blank(),panel.border = element_rect(colour = "grey40",fill=NA), axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~ times, scale="free_y", nrow=1)

  p_quality_single <- ggplot(TotQuality, aes(x=sample, y=y, fill=sample)) + 
  geom_violin(trim=FALSE,draw_quantiles = c(0.25, 0.5, 0.75), show.legend=FALSE)+
  #geom_boxplot(width=0.1, fill="white", show.legend=FALSE)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange", color="black",show.legend=FALSE)+
  geom_jitter(shape=16, position=position_jitter(.05), size=.5, fill="grey30", show.legend=FALSE)+ ## check
  scale_fill_brewer(palette="RdBu")+
  labs(x="",y=expression("quality "["(phred)"]))+
  theme_minimal()+
  theme(strip.background = element_blank(),panel.border = element_rect(colour = "grey40",fill=NA), axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~ times, scale="free_y", nrow=1)

  pdf(file.path(Directory,'Violins_All.pdf'), height=5*length(DataIn), width=7*length(DataIn), onefile=TRUE) #teorethically, up to infinite number of experiments can be compared
    
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(4, 1)))
  print(p_reads, vp=define_region(1, 1))
  print(p_bases, vp = define_region(2, 1))
  print(p_length, vp = define_region(3, 1))
  print(p_quality,vp = define_region(4, 1))

  dev.off()

  pdf(file.path(Directory,'Violins_Reads.pdf'), height=5, width=7*length(DataIn), onefile=TRUE) #teorethically, up to infinite number of experiments can be compared

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 1)))
  print(p_reads, vp=define_region(1, 1))
  dev.off()

  pdf(file.path(Directory,'Violins_Bases.pdf'), height=5, width=7*length(DataIn), onefile=TRUE) #teorethically, up to infinite number of experiments can be compared

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 1)))
  print(p_bases_single, vp=define_region(1, 1))
  dev.off()

  pdf(file.path(Directory,'Violins_Length.pdf'), height=5, width=7*length(DataIn), onefile=TRUE) #teorethically, up to infinite number of experiments can be compared

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 1)))
  print(p_length_single, vp=define_region(1, 1))
  dev.off()


  pdf(file.path(Directory,'Violins_Quality.pdf'), height=5, width=7*length(DataIn), onefile=TRUE) #teorethically, up to infinite number of experiments can be compared

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 1)))
  print(p_quality_single, vp=define_region(1, 1))
  dev.off()




  ###SPEEDTEST ### 

  ## it would be nice to have a plot for speed test in terms of bps productivity but for the moment just rely on the others for a complete overview



  message("Done")
}
