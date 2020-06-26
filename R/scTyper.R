#' @title scTyper
#' @description Run scTyper
#' @usage scTyper(seurat.object, marker, wd, output.name, pheno.fn, qc = FALSE, run.cellranger=FALSE, norm.seurat=FALSE, cell.typing.method, level, run.inferCNV=TRUE, proj.name, fastqc.path, fastq.dir, fq1.idx, fq2.idx, cellranger.path, cellranger.ref.dir=NULL, percent.min.cells=0.1, min.features=200, percent.mt=10, vars.to.regress=c("nCount_RNA", "percent.mt"), dims=1:10, resolution=2, slot=c("scale.data", "count.data", "data"), assay="RNA", NTP.g.filter.method=c("sd", "mad", "none"), NTP.gene.filter.cutoff=0.3, NTP.distance=c("cosine","correlation"), NTP.norm.method, gene.ref.gtf=NULL, feature.to.test = c("tissue.type", "cell.type"), cells.test_excluded="Epithelial", cells.test_reference = "immune", fc.cutoff=0.05, cutoff.gene.cluster=0.05, malignant.cell.type="Epithelial", report.mode=TRUE, mc.cores = 1)
#' @param seurat.object Seurat object, if users have pre-processed seurat object, user have to insert seurat object as input
#' @param marker Cell markers to use in cell typing, character or List (identifier or StudyName or User defined gene list)
#' @param wd Working directory
#' @param output.name Output directory name
#' @param pheno.fn Phenotype file path
#' @param qc  Whether to execute FASTQC (default=FALSE)
#' @param run.cellranger whether to excute cellranger count (default=FALSE)
#' @param norm.seurat whether to normalize seurat object (default=FALSE)
#' @param cell.typing.method cell typing method, c("NTP", "ES", "Average"), (default = "NTP")
#' @param level Indicate the cell assignment level (cell or cluster)
#' @param run.inferCNV Indicate whether ‘malignant cell typing by inferCNV process run
#' @param proj.name Project name
#' @param fastqc.path FastQC program path
#' @param fastq.dir FastQC output directory
#' @param fq1.idx Index of the FASTQ file (Read 1)
#' @param fq2.idx Index of the FASTQ file (Read 2)
#' @param cellranger.path Cell Ranger program path
#' @param cellranger.ref.dir Directory of Cell Ranger reference file
#' @param percent.min.cells Cutoff to filter features containing minimum percent of cells
#' @param min.features Cutoff to filter cells containing minimum number of features
#' @param percent.mt Cutoff for filtering cells that have >n percent mitochondrial counts
#' @param vars.to.regress Variables to regress out
#' @param dims A vector of the dimensions to use in construction of the SNN graph.
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param slot Data type of Seurat object, c("scale.data", "count.data", "data")
#' @param assay Assay of Seurat object
#' @param NTP.g.filter.method Method to filter genes in NTP
#' @param NTP.gene.filter.cutoff Cutoff to filter genes of in NTP
#' @param NTP.distance NTP distance method, a character, either c("correlation" or "cosine").
#' @param NTP.norm.method NTP normalization method, either c("none", "row.std")
#' @param gene.ref.gtf Path of GTF file including genomic location for genes
#' @param feature.to.test Column header name of the meta data in Seurat object (select the cell groups for T.test) either "tissue.type" or "cell.type"
#' @param cells.test_excluded A value indicates the cells to be excluded in T.test
#' @param cells.test_reference A value indicates the cells to use as be excluded in T.test
#' @param fc.cutoff Cutoff of fold change
#' @param cutoff.gene.cluster A cutoff P-value for filtering out the gene clusters (calculated from GO analysis)
#' @param malignant.cell.type Cell type to assign malignant cell
#' @param report.mode Generate report file
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @return Seurat object
#' @export
scTyper<- function(seurat.object=NULL,
                   marker="Puram.2017.HNSCC.TME",

                   wd=getwd(),
                   output.name = "test.result",
                   pheno.fn,
                   qc = FALSE,
                   run.cellranger=FALSE,
                   norm.seurat=FALSE,
                   cell.typing.method=c("NTP", "ES", "Average"),
                   level=c("cell","cluster"),
                   run.inferCNV=TRUE,
                   proj.name = "scTyper",

                   fastqc.path=NULL,
                   fastq.dir=NULL,
                   fq1.idx="_R1_001.fastq",
                   fq2.idx="_R2_001.fastq",

                   cellranger.path=NULL,
                   cellranger.ref.dir=NULL,

                   percent.min.cells=0.1,
                   min.features=200,
                   percent.mt=10,
                   vars.to.regress=c("nCount_RNA", "percent.mt"),
                   dims=1:10,
                   resolution=2,

                   slot=c("scale.data", "count.data", "data"),
                   assay="RNA",
                   NTP.g.filter.method=c("sd", "mad", "none"),
                   NTP.gene.filter.cutoff=0.3,
                   NTP.distance=c("cosine","correlation"),
                   NTP.norm.method=c("none", "row.std"),

                   gene.ref.gtf=NULL,
                   feature.to.test = c("tissue.type", "cell.type"),
                   cells.test_excluded="Epithelial",
                   cells.test_reference = "immune",
                   fc.cutoff=0.05,
                   cutoff.gene.cluster=0.05,
                   malignant.cell.type="Epithelial",
                   report.mode=TRUE,
                   mc.cores = 1){

  start_time <- Sys.time()
  cell.typing.method=match.arg(cell.typing.method)
  level=match.arg(level)
  slot=match.arg(slot)
  NTP.gene.filter.method=match.arg(NTP.g.filter.method)
  NTP.distance=match.arg(NTP.distance)
  NTP.norm.method=match.arg(NTP.norm.method)
  feature.to.test=match.arg(feature.to.test)

  library(fastqcr)
  library(Seurat)
  library(infercnv)
  library(gProfileR)
  library(rmarkdown)
  library(parallel)
  library(perm)
  library(Biobase)
  library(reshape2)
  library(limma)
  library(GenomicRanges)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(grDevices)
  library(ComplexHeatmap)
  library(circlize)
  library(png)
  library(knitr)
  library(kableExtra)
  library(pander)
  library(colorspace)

  # default outdirs
  qc.dir <- file.path(wd, output.name, "00_qc")
  count.dir <- file.path(wd, output.name,　"01_count")
  ntp.dir <- file.path(wd, output.name,　"02_NTP")
  inferCNV.dir  <- file.path(wd, output.name,　"03_inferCNV")
  rda.dir <- file.path(wd, output.name,　"data")



  #######################
  ### start scTyper
  dir.create(file.path(wd, output.name), showWarnings = F)
  dir.create(file.path(rda.dir), showWarnings = F)

  pheno.df = read.csv(pheno.fn, header = T, sep = ",")
  sample.name = as.character(pheno.df$Sample_ID)

  if(qc==TRUE){
    fastqc.outs = fastqc(fastqc.path, fastq.dir=fastq.dir, sample.name, fq1.idx=fq1.idx, fq2.idx=fq2.idx, output.dir=qc.dir, run.cmd=TRUE, mc.cores=mc.cores)
  }

  if(run.cellranger==TRUE){
    count.out.dirs = CellrangerCount(cellranger.path, fastq.dir=fastq.dir, output.dir=count.dir, cellranger.ref.dir=cellranger.ref.dir, sample.name=sample.name, run.cmd=TRUE, mc.cores=mc.cores)
    setwd(wd)
    metrics_summary = make.stat_summary(count.dir = count.dir, sample.name = sample.name, pheno.df=pheno.df, output.dir = count.dir)
  }

  if(norm.seurat==TRUE){
    seurat = run.seurat.process(count.dir = count.dir, rda.dir = rda.dir, project = proj.name, sample.name = sample.name, metrics_summary,
                                percent.min.cells=percent.min.cells, min.features=min.features, percent.mt=percent.mt, vars.to.regress=vars.to.regress, dims=dims,
                                resolution=resolution)
  }

  if(class(seurat.object)=="Seurat") seurat=seurat.object

  ##cell typing
  seurat=cell.typing.seurat(seurat=seurat, marker=marker, cell.typing.method=cell.typing.method, level=level, wd=wd, slot=slot, assay=assay, ntp.dir=ntp.dir, rda.dir=rda.dir, NTP.g.filter.method=NTP.g.filter.method, NTP.gene.filter.cutoff=NTP.gene.filter.cutoff, NTP.distance=NTP.distance, NTP.norm.method=NTP.norm.method, mc.cores=mc.cores)

  if(run.inferCNV==TRUE){
    seurat@meta.data$tissue.type = pheno.df[match(seurat@meta.data$sample.name, pheno.df[,"Sample_ID"]),"TissueType"]
    fdata = make.seurat.fdata(seurat, rda.dir, gene.ref.gtf = gene.ref.gtf)

    seurat = run.inferCNV(seurat = seurat, assay=assay, fdata, pheno_info = pheno.df, output.dir = inferCNV.dir, rda.dir = rda.dir, feature.to.test=feature.to.test, cells.test_reference=cells.test_reference, cells.test_excluded=cells.test_excluded, fc.cutoff=fc.cutoff, cutoff.gene.cluster=cutoff.gene.cluster, mc.cores=mc.cores)

    seurat = malignant.cellTyper(seurat = seurat, rda.dir = rda.dir, cells.test_reference=cells.test_reference, malignant.cell.type=malignant.cell.type)

  }
  end_time <- Sys.time()
  runtime = format(end_time-start_time, format="%H")
  #report
  envList=as.list(environment())
  if(report.mode) {report(envList, qc.dir, output.dir=file.path(wd, output.name))}

  return(seurat)
}



