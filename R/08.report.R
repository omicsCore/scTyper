
#' @title report
#' @description Reports the result of using scTyper()
#' @usage report(envList, qc.dir, output.dir)
#' @param envList R environment list
#' @param qc.dir qc directory
#' @param output.dir output directory
#' @details Provides a report that summarizes the processing steps and visualized tables and plots for the processed results. The report file is automatically generated recording the workflows of the data processing steps, the options used in the processing, and the outcome results.
#' @return pdf file include data processing result information
#' @import rmarkdown, fastqcr
#' @export
report=function(envList, qc.dir, output.dir){
  message("[[",Sys.time(),"]] Run to make report --------")

  #param
  global.conf.param=c("proj.name", "wd", "output.name", "pheno.fn", "qc","run.cellranger","norm.seurat", "cell.typing.method", "level", "run.inferCNV",  "mc.cores", "report.mode")
  qc.param=c("fastq.dir", "fastqc.path", "fq1.idx", "fq2.idx")
  cellranger.param=c("cellranger.path", "cellranger.ref.dir")
  normalization.param=c("percent.min.cells", "min.features", "percent.mt", "vars.to.regress", "dims", "resolution")
  celltyping.param=c("marker", "slot", "assay")
  celltyping.ntp.param=c("NTP.g.filter.method", "NTP.gene.filter.cutoff", "NTP.distance", "NTP.norm.method")
  malig.celltyping.param=c("gene.ref.gtf", "feature.to.test", "cells.test_excluded", "cells.test_reference", "fc.cutoff", "malignant.cell.type", "cutoff.gene.cluster")


  #param df
  global.conf.param.df=data.frame(Parameters=global.conf.param, Values=as.character(envList[global.conf.param]), row.names = NULL)
  qc.param.df=data.frame(Parameters=qc.param, Values=as.character(envList[qc.param]), row.names = NULL)
  cellranger.param.df=data.frame(Parameters=normalization.param, Values=as.character(envList[normalization.param]), row.names = NULL)
  normalization.param.df=data.frame(Parameters=normalization.param, Values=as.character(envList[normalization.param]), row.names = NULL)
  scTyping.param=c(celltyping.param, celltyping.ntp.param, malig.celltyping.param)
  scTyping.param.df=data.frame(Parameters=scTyping.param, Values=as.character(envList[scTyping.param]), row.names = NULL)

  if(class(envList$seurat)=="Seurat") input.data="Seurat object" else input.data="raw data (.fastq)" #input data type


  input.param.names=c("seurat", "runtime", "qc.dir", "metrics_summary", global.conf.param, qc.param, cellranger.param, normalization.param, celltyping.param, celltyping.ntp.param, malig.celltyping.param)

  message("===Input Parameters===")
  if(envList$qc==T & envList$run.cellranger==T){
    qc.res=get.qc.report(qc.dir)
    input.params=envList[input.param.names]
    input.params=c(input.params,qc.res)
  }else{
    input.params=envList[input.param.names]
  }

  all.pipeline=c("QC", "CellRanger", "Normalization", "Cell Typing", "inferCNV")
  all.param=c("qc", "run.cellranger", "norm.seurat", "cell.typing.method", "run.inferCNV")
  all.param.val=as.character(unlist(envList[all.param]))
  pipeline.des=sapply(1:length(all.pipeline), function(i) paste(all.pipeline[i], " (",all.param.val[i], ")", sep = ''))
  pipeline.des=paste(pipeline.des,collapse = ", ")

  report.envList=list(global.conf.param.df=global.conf.param.df, qc.param.df=qc.param.df, cellranger.param.df=cellranger.param.df, normalization.param.df=normalization.param.df, scTyping.param.df=scTyping.param.df, pipeline.des=pipeline.des, input.data=input.data)
  input.params=c(input.params, report.envList)
  print(input.params)

  message("===Output directory===")
  output.dir = file.path(input.params$wd, input.params$output.name)

  message("===Result===")

  # summary
  report.summary=input.params
  save(report.summary, file = file.path(report.summary$wd, report.summary$output.name, "data/report.summary.rda"))
  report.rmd.fn=system.file("report/scTyper_report.Rmd", package = "scTyper")

  render(report.rmd.fn, output_format="html_document", output_dir=output.dir)

  message("[[",Sys.time(),"]] Finish making report file --------")
}

