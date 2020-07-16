
#' @title CellrangerCount
#' @description A wrapper function to run Cellranger Count process
#' @param cellranger.path Cell Ranger program path
#' @param fastq.dir FastQC output directory
#' @param cellranger.ref.dir Directory of Cell Ranger reference file
#' @param output.dir Output directory
#' @param sample.name sample name
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details CellrangerCount takes FASTQ files from fastQC and performs alignment, filtering, barcode counting, and UMI counting.
#' @return feature-barcode matrices and Secondary analysis (e.g., dimensionality reduction, cell clustering, and differential expression)
#' @import parallel
#' @references Massively parallel digital transcriptional profiling of single cells. GXY Zheng. (2017).
#' @seealso \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count}
#' @export
CellrangerCount <- function(cellranger.path, fastq.dir, cellranger.ref.dir, output.dir, sample.name, run.cmd=TRUE, mc.cores=1){
  fqs=file.path(fastq.dir, sample.name)
  # command
  cmd=paste0(cellranger.path, " count ",  "--id=", sample.name, " --fastqs=", fqs, " --transcriptome=", cellranger.ref.dir)

  #make directory
  dir.create(output.dir, showWarnings=FALSE)

  setwd(output.dir)
  # run
  message("[[",Sys.time(),"]] Run CellrangerCount --------")
  message(cmd)
  if(run.cmd) lapply(cmd, system)
  cat(cmd, file=file.path(output.dir, "run.CellrangerCount.log"), sep="\n", append=TRUE)
  message("[[",Sys.time(),"]] CellrangerCount finished --------")

  out.dirs=file.path(output.dir, sample.name)
  return(out.dirs)
}


#' @title make.stat_summary
#' @description A wrapper function to run make.stat_summary
#' @usage make.stat_summary(count.dir, sample.name, output.dir, pheno.df)
#' @param count.dir cellragner count ouput directory
#' @param sample.name cell sample name
#' @param output.dir Output directory
#' @param pheno.df phenotype dataframe(reference an instruction manual)
#' @details make data summary file
#' @return data summary file
#' @export
make.stat_summary <- function(count.dir, sample.name, output.dir, pheno.df){
  message("[[",Sys.time(),"]] Run make stat summary --------")

  fns=file.path(count.dir, sample.name, "outs", "metrics_summary.csv")
  summary=lapply(fns, read.csv)

  summary.df=data.frame(t(sapply(summary, function(a) as.numeric(sapply(1:ncol(a), function(i) gsub(",|%", "", a[1,i]))))))
  colnames(summary.df)=colnames(summary[[1]])
  rownames(summary.df)=sample.name

  summary.df=data.frame(sample.name=rownames(summary.df), summary.df)
  cols=c("Number.of.Reads", "Q30.Bases.in.Barcode", "Q30.Bases.in.RNA.Read", "Q30.Bases.in.Sample.Index", "Q30.Bases.in.UMI")


  write.csv(summary.df, file = file.path(output.dir,"metrics_summary_all.csv"), row.names = F)
  message("[[",Sys.time(),"]] Finished stat summary --------")
  return(summary.df)
}


