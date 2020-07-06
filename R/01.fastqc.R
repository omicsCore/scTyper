
#' @title fastqc
#' @description A wrapper function to run fastQC
#' @usage fastqc(fastqc.path, fastq.dir, sample.name, fq1.idx="_R1_001.fastq", fq2.idx="_R2_001.fastq", output.dir, run.cmd=TRUE, mc.cores=1)
#' @param fastqc.path FastQC program path
#' @param fastq.dir FastQC output directory
#' @param sample.name sample name
#' @param fq1.idx Index of the FASTQ file (Read 1)
#' @param fq2.idx Index of the FASTQ file (Read 2)
#' @param output.dir Output directory
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details FastQC aims to provide a QC report that detects problems originating from either the sequencer or the starting library material.
#' @return Quality check report for sequence data. (e.g., .html)
#' @import parallel
#' @references FastQC: a quality control tool for high throughput sequence data. Andrews S. (2010).
#' @seealso \url{http://www.bioinformatics.babraham.ac.uk/projects/fastqc}
#' @export
fastqc=function(fastqc.path, fastq.dir, sample.name, fq1.idx="_R1_001.fastq", fq2.idx="_R2_001.fastq", output.dir, run.cmd=TRUE, mc.cores=1){
  fastq.sample.dir=file.path(fastq.dir,sample.name)
  fq1=list.files(fastq.sample.dir, fq1.idx, full.names = TRUE)
  fq2=list.files(fastq.sample.dir, fq2.idx, full.names = TRUE)

  fq2=ifelse(grepl(fq2.idx, fq2), fq2, "")
  # command
  cmd=paste(fastqc.path, "-o", output.dir, "--extract", fq1, fq2)

  #make directory

  fastqc.out=file.path(output.dir)
  dir.create(fastqc.out,showWarnings=FALSE)
  #  setwd(output.dir)
  # run
  message("[[",Sys.time(),"]] Run fastQC --------")

  message(cmd)

  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "run.fastqc.log"), sep="\n", append = FALSE)
  list.files(output.dir, full.names=TRUE)
  message("[[",Sys.time(),"]] Finished fastQC --------")
}

#' @title get.qc.report
#' @description QC Reports
#' @usage get.qc.report(qc.dir)
#' @param qc.dir qc directory
#' @details Provides FASTQC report that summarizes the QC processing steps
#' @import rmarkdown, fastqcr
#' @export
get.qc.report=function(qc.dir){
  suppressPackageStartupMessages(library(fastqcr))
  fns= dir(qc.dir, "fastqc.zip$", full.names = TRUE)

  res.summary=data.frame(qc_stats(qc_aggregate(qc.dir, progressbar = FALSE)))[,-2]

  qual.scores=lapply(fns, function(a) data.frame(qc_read(a, modules = "per sequence quality scores", verbose = FALSE)))

  res.summary$Q30.perc=sapply(qual.scores, function(b) sum(b[which(b[,1]>30),2])/sum(b[,2])*100)

  colnames(res.summary)=c("Filename","GC%","Total Sequence","Sequence length","Phred Score(>30)(%)")

  list(qc=res.summary)
}



#' @title fastqc.table
#' @description A wrapper function to make fastqc dataframe
#' @usage fastqc.table(qc.dir)
#' @param qc.dir fastqc output directory
#' @details make table of fastqc using qc outputs.
#' @return fastqc dataframe
#' @export
fastqc.table <- function(qc.dir){
  qc.res=get.qc.report(qc.dir)
  fastqc.df=qc.res$qc
  fastqc.df[,3]=format(as.numeric(fastqc.df[,3]),big.mark = ",")
  colnames(fastqc.df)[2]="GC(%)"
  colnames(fastqc.df)[3]='Total\nReads'
  colnames(fastqc.df)[4]='Read\nlength(bp)'
  colnames(fastqc.df)[5]='Phred Score\n(>30)(%)'
  colnames(fastqc.df)[1]='Filename'
  return(fastqc.df)
}


#' @title fastqc.summary
#' @description A wrapper function to make fastqc dataframe
#' @usage fastqc.summary(fastqc.df, fq1.idx="_R1_001.fastq", fq2.idx="_R2_001.fastq")
#' @param fastqc.df fastqc output directory
#' @param fq1.idx Index of the FASTQ file (Read 1)
#' @param fq2.idx Index of the FASTQ file (Read 2)
#' @details make table of fastqc using qc outputs.
#' @return fastqc dataframe
#' @export
fastqc.summary <- function(fastqc.df, fq1.idx="_R1_001.fastq", fq2.idx="_R2_001.fastq"){
  ### 7.3.2 FastQC summary

  colnames(fastqc.df)[5]=c("PhredScore")

  fastqc.df$Reads=NA
  fastqc.df$Reads[grep(".1$|_1.fq$", fastqc.df[,1])]=sub(".1$|_1.fq$","",fastqc.df$Filename[grep(".1$|_1.fq$", fastqc.df[,1])])
  fastqc.df$Reads[grep(".2$|_2.fq$", fastqc.df[,1])]=sub(".2$|_2.fq$","",fastqc.df$Filename[grep(".2$|_2.fq$", fastqc.df[,1])])

  if(length(unique(fastqc.df[,1]))>10) coord_flip=coord_flip() else coord_flip=NULL
  if(length(unique(fastqc.df[,1]))<10) opts=theme(axis.text.x=element_text(angle=90 , hjust = 1)) else opts=NULL

  plot.fastqc.df=fastqc.df[rev(rownames(fastqc.df)),]
  plot.fastqc.df$Filename=as.factor(plot.fastqc.df$Filename)


  g=ggplot(data=plot.fastqc.df, aes(x=Filename, y=PhredScore, fill=Reads)) + geom_bar(stat="identity", position=position_dodge() ,width=0.5)

  print(g+theme_minimal()+coord_flip+scale_fill_manual(values=rep(c("skyblue", "steelblue"),length(unique(fastqc.df[,6]))))+ylab("Proportion(%)") + opts + theme(legend.position="none",axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold",size=7)) + xlab("File Name") + labs(fill = "") +ggtitle("Proportion of reads with a quality score of 30 or higher") + theme(plot.title = element_text(hjust = 0.5,size = 18))+scale_x_discrete(limits = rev(levels(plot.fastqc.df$Filename))))

}

