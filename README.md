# NanoR

![alt text](NanoR.png)

NanoR is an up-to-date package for the statistical language and environment R, tested on Unix, MacOSX and Windows, that allows user-friendly analysis and comparison of 1D MinION and GridION X5 sequencing data within acceptable time frames.


If you are using NanoR, please consider citing "NanoR: A user-friendly R package to analyze and compare nanopore sequencing data" published in PLOS ONE.

#### Citation:

Bolognini D, Bartalucci N, Mingrino A, Vannucchi AM, Magi A (2019) NanoR: A user-friendly R package to analyze and compare nanopore sequencing data. PLoS ONE 14(5): e0216471. https://doi.org/10.1371/journal.pone.0216471


##### NanoR v 2.0 is out !

- Added support for multi-read .fast5 files
- Changed graphics design: NanoR goes minimal !
- Plotting is now faster
- NanoCompare is now faster (and histograms are not plotted anymore)
- The "M" version of NanoR now supports analysis of basecalled .fast5 files from MinION and GridION while the "G" version now supports analysis of sequencing summaries and .fastq files from the 2 platforms
- Added the possibility to extract/filter .fastq files on a minimum quality treshold
- Added the possibility to store tables behind ggplot2-plots (easier to hack colors ;))
- Removed seqinr dependency


## Background

NanoR was developed under R 3.1.3 and tested on Unix, MacOSX and Windows. It supports all the releases of MinION and GridION X5 instruments (and their outputs, of course !)


NanoR depends on the following R packages:

- ggplot2 (v 2.2.1)
- reshape2 (v 1.4.3)
- RColorBrewer (v 1.1.2)
- scales (v 0.5.0)
- gridExtra (v 2.3)
- rhdf5 (v 2.14)
- ShortRead (v 1.24.0)
- parallel (v 3.1.3)
- grid (v 3.1.3)

If you don't have these packages in your R library, you have to download and install them before intalling NanoR, as specified below.


In order to install the needed packages, you can copy, paste and run the following lines in your R console:


```R

install.packages(c("ggplot2","reshape2","ShortRead","RColorBrewer","scales","gridExtra"), repos= "http://cran.cnr.berkeley.edu/")
 
source("http://bioconductor.org/biocLite.R")

biocLite("rhdf5",suppressUpdates=TRUE, suppressAutoUpdate=TRUE)

biocLite("ShortRead",suppressUpdates=TRUE, suppressAutoUpdate=TRUE)

```
As rhdf5 version 2.14 (or higher) is needed, if biocLite does not automatically download this version, you can download it manually from bioconductor (http://bioconductor.org/packages/3.2/bioc/src/contrib/rhdf5_2.14.0.tar.gz) and install in R running from your R console:

```R

install.packages("/Path/To/rhdf5_2.14.0.tar.gz", repos=NULL)

```


Packages "parallel" and "grid" are "base packages" (they are usually present in your package list) but, if you don't find them in your available packages, quit R and add this line to your .Rprofile (or Rprofile.site) file:

```sh

options(defaultPackages=c(getOption("defaultPackages"),"parallel","grid"))

```

After that, save your new .Rprofile (or Rprofile.site) file and re-run R.
In order to install NanoR, you have to download the .tar.gz file from this repository using git:

```sh

git clone https://github.com/davidebolo1993/NanoR.git

```

If you have any issues with this command, try:

```sh

git clone git://github.com/davidebolo1993/NanoR.git

```

The previous command will clone in your folder the last release of NanoR (current, NanoR v2.0): check the releases to download R .tar.gz for previous versions. 


Then, from your R console:

```R

install.packages("/Path/To/NanoR.tar.gz", repos=NULL)

```

If you're using NanoR on Windows, be sure to specify paths using the following format, in order to avoid problems with "escape characters":

_"C:\\\Path\\\To\\\Data"_

Use NanoR on the complete set of .fast5 files you obtain from Nanopore MinION/GridION X5 in order to obtain significant statistics. A MinION and GridION X5 sample dataset on which NanoR can be tested is available to be downloaded at _https://faspex.embl.de/aspera/faspex/external_deliveries/823?passcode=c4aa66d379850be32163f8517fa04750aec7a3bf&expiration=MjAxOS0wNS0wMlQwNjo1Nzo0MFo=_. If the first link does not work, a more stable one is:
_https://oc.embl.de/index.php/s/JraNyCrcJhStcNB_. Get in touch with me if these links do not work anymore.



## Workflow examples

You can access all the informations on how to run functions from NanoR within R, using the following commands:

MinION data analysis

```R

?NanoPrepareM

?NanoTableM

?NanoStatsM

?NanoFastqM

```

GridION X5 data analysis

```R

?NanoPrepareG

?NanoTableG

?NanoStatsG

?NanoFastqG

```

Comparison between experiments

```R

?NanoCompare

```


Here is an example of how to run the aforementioned commands:



### MinION and GridION X5 basecalled .fast5 files, single/multi-read (M version)

```R


List<-NanoPrepareM(DataPass="/Path/To/PassedFast5Files",DataFail=NA,DataSkip=NA,Label="Exp", MultiRead=FALSE) # prepare data. To allow multi-read .fast5 files support simply switch MultiRead to TRUE

Table<-NanoTableM(NanoMList=List,DataOut="/Path/To/DataOut",Cores=6,GCC=FALSE) # extract metadata. To allow GC content computation, switch GCC to TRUE

NanoStatsM(NanoMList=List,NanoMTable=Table,DataOut="/Path/To/DataOut", KeepGGObj = FALSE) #plot statistics. To store table behind ggplot2-plots, switch KeepGGObj to TRUE

NanoFastqM(DataPass="/Path/To/PassedFast5Files",DataOut="/Path/To/DataOut",Label="Exp",Cores=6,FASTA=FALSE, Minquality=7, MultiRead=FALSE) # extract .fastq. To convert .fastq to .fasta as well, switch FASTA to TRUE; to extract .fastq only from .fast5 files with quality greater or equal than Minquality, increase the Minquality parameter; to allow support for multi-read .fast5 files, switch MultiRead to TRUE.


```

If working with folders containing passed, failed and skipped .fast5 files together, give this folder to the "DataPass" parameter of NanoPrepareM and NanoFastqM: NanoR will automatically filter out the low-quality sequences




###  MinION and GridION X5 sequencing summary and .fastq files (G version)

```R

List<-NanoPrepareG(DataSummary="/Path/To/DataSummary", DataFastq="/Path/To/DataFastq", Cores = 1, Label="Exp") #prepare data. Using multiple cores is only useful when dealing with multiple sequencing summary files (old behaviour of GridION X5)

Table<-NanoTableG(NanoGList=List,DataOut="/Path/To/DataOut",GCC=FALSE) #arrange metadata. To extract GC content from .fastq files, switch GCC to TRUE

NanoStatsG(NanoGList=List,NanoGTable=Table,DataOut="/Path/To/DataOut", KeepGGObj = FALSE) # plot statistics.  To store table behind ggplot2-plots, switch KeepGGObj to TRUE

NanoFastqG(DataSummary="/Path/To/DataSummary", DataFastq="/Path/To/DataFastq", DataOut="/Path/To/DataOut", Cores = 1, Label="Exp", FASTA=FALSE, Minquality = 7) #filter .fastq file on a minimum quality defined in Minquality. To filter .fastq files on higher quality, increase Minquality treshold; to convert .fastq to .fasta as well, switch FASTA to TRUE.Using multiple cores is only useful when dealing with multiple sequencing summary files (old behaviour of GridION X5)

```

Some of the plots generated by NanoR for MinION and GridION X5 data analysis are shown below:

**Reads number, basepairs number, reads length and reads quality for every 30 minutes of experimental run**

![alt text](Plots/RBLQ.png)

**Reads length vs reads quality**

![alt text](Plots/LvsQ.png)

**Channels and muxes activity**

![alt text](Plots/Activity.png)


In addition, NanoR plots the yield (cumulative reads and cumulative base pairs), the number and percentage of failed/passed reads and the GC content count histogram (if GCC has been computed).

### Compare MinION/GridION X5 data

```R

DataIn<-c("Path/To/AnalyzedFolder1/DataForComparison","Path/To/AnalyzedFolder2/DataForComparison","Path/To/AnalyzedFolder3/DataForComparison" #path to the NanoR-analyzed data

Labels<-c("Label1","Label2","Label3") #labels to identify the experiments

NanoCompare(DataIn=DataIn,DataOut="Path/To/DataOut",Labels=Labels) #compare

```

NanoCompare returns a Violins.pdf plot that compare reads number, base pairs number, reads length and reas quality every 10 hours of experimental run (maximum 80 hours of experimental run)
