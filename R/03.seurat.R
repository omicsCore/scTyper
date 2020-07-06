
#' @title run.seurat.process
#' @description A wrapper function to run run.seurat.process
#' @usage run.seurat.process(count.dir = count.dir, rda.dir = rda.dir, project = proj.name, metrics_summary, sample.name = sample.name, percent.min.cells = 0.01, min.features = 200, scale.factor = 10000, vars.to.regress=c("nCount_RNA", "percent.mt"), selection.method = "vst", more_nFeature_RNA = 200, Less_nFeature_RNA = 8000, percent.mt = 10, normalize=TRUE, assay = 'RNA', dims=1:10, resolution=seq(0.6, 2, 0.2), random.seed=1234)ss_nFeature_RNA = 8000, percent.mt = 10, normalize=TRUE, assay = 'RNA', clust.dims=1:10, resolution=seq(0.6, 2, 0.2), random.seed=1234)
#' @param count.dir Ouput directory of cellragner count
#' @param rda.dir Path of the RData saving directory
#' @param project project name
#' @param metrics_summary cellranger summary metrics
#' @param sample.name Sample name
#' @param percent.min.cells Include genes with detected expression in at least this many cells. Will subset the raw.data matrix as well. To reintroduce excluded genes, create a new object with a lower cutoff.
#' @param min.features Include cells where at least this many genes are detected.
#' @param scale.factor Sets the scale factor for cell-level normalization(10,000 by default)
#' @param vars.to.regress Variables to regress out (previously latent.vars in RegressOut). For example, nUMI, or percent.mito.
#' @param selection.method How to choose top variable features. Choose one of 'vst', 'mean.var.plot', 'dispersion'
#' @param more_nFeature_RNA High cutoffs for filtering cells that have unique feature counts (default is 200)
#' @param Less_nFeature_RNA low cutoffs for filtering cells that have unique feature counts (default is 8000)
#' @param percent.mt low cutoffs for filtering cells that have >n percent mitochondrial counts (default is 10)
#' @param normalize use log normalization
#' @param assay Assay to use
#' @param dims A vector of the dimensions to use in construction of the SNN grouph.
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param random.seed Seed of the random number generator.
#' @details CellrangerCount takes FASTQ files from fastQC and  performs alignment, filtering, barcode counting, and UMI counting.
#' @return feature-barcode matrices and Secondary analysis (e.g., dimensionality reduction, cell clustering, and differential expression)
#' @references Massively parallel digital transcriptional profiling of single cells. GXY Zheng. (2017).
#' @import seurat, limma
#' @export
run.seurat.process <- function(count.dir = count.dir,
                               rda.dir = rda.dir,
                               project = proj.name,
                               metrics_summary,
                               sample.name = sample.name,
                               percent.min.cells = 0.01,
                               min.features = 200,
                               scale.factor = 10000,
                               vars.to.regress=c("nCount_RNA", "percent.mt"),
                               selection.method = "vst",
                               more_nFeature_RNA = 200,
                               Less_nFeature_RNA = 8000,
                               percent.mt = 10,
                               normalize=TRUE,
                               assay = 'RNA',
                               dims=1:10,
                               resolution=seq(0.6, 2, 0.2),
                               random.seed=1234){

  # run
  message("[[",Sys.time(),"]] make.seurat --------")
  seurat = make.seurat(count.dir,
                       sample.name = sample.name,
                       project = project,
                       min.cells = percent.min.cells,
                       min.features = min.features)

  message("[[",Sys.time(),"]] Run Seurat filtering and normalization --------")
  seurat = cell.filter.seurat(seurat,
                              sample.name,
                              metrics_summary,
                              more_nFeature_RNA ,
                              Less_nFeature_RNA,
                              percent.mt=percent.mt )
  if(normalize==TRUE){
    seurat = NormalizeData(object =seurat, normalization.method = "LogNormalize", scale.factor = scale.factor)
  }
  seurat = ScaleData(seurat, vars.to.regress=vars.to.regress)


  message("[[",Sys.time(),"]] Identify highly variable features --------")
  seurat = FindVariableFeatures(seurat, selection.method = selection.method, assay = assay)

  message("[[",Sys.time(),"]] Perform linear dimensional reduction --------")
  seurat = RunPCA(seurat, features = VariableFeatures(object = seurat), verbose = FALSE)

  message("[[",Sys.time(),"]] Make cluster of the cells --------")
  seurat = FindNeighbors(seurat, dims = dims, random.seed = random.seed)
  seurat = FindClusters(seurat, random.seed = random.seed)


  message("[[",Sys.time(),"]] Run non-linear dimensional reduction(t-SNE) --------")
  seurat = RunTSNE(object = seurat, seed.use=random.seed)

  #make directory
  dir.create(rda.dir, showWarnings=FALSE)
  save(seurat, file=paste0(rda.dir,"/seurat.rda"))
  save(seurat, file = file.path(rda.dir, "seurat.norm.rda"))
  message("[[",Sys.time(),"]] Seurat process finished --------")
  print(seurat)
  return(seurat)
}


#' @title make.seurat
#' @description A wrapper function to make.seurat
#' @usage make.seurat(count.dir, sample.name = sample.name, project = "SeuratProject", min.cells=0, min.features=0)
#' @param count.dir Path of the cellranger count directory
#' @param sample.name single cell RNA sequening sample name
#' @param project project name(string)
#' @param min.cells Include genes with detected expression in at least this many cells. Will subset the raw.data matrix as well. To reintroduce excluded genes, create a new object with a lower cutoff.
#' @param min.features Include cells where at least this many genes are detected.
#' @details Initializes the Seurat object and some optional filtering
#' @return  Seurat object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset.
#' @export
make.seurat <- function(count.dir,
                        sample.name = sample.name,
                        project = "SeuratProject",
                        min.cells=0,
                        min.features=0){

  fns=file.path(count.dir, sample.name, "outs", "filtered_gene_bc_matrices/GRCh38")


  expression_matrixList=lapply(fns, Read10X)
  expression_matrixList=lapply(1:length(expression_matrixList), function(i) {
    colnames(expression_matrixList[[i]])=paste0(colnames(expression_matrixList[[i]]), "-", i)
    expression_matrixList[[i]]
  })
  colnames(expression_matrixList[[1]])
  sum(sapply(expression_matrixList, ncol))
  expression_matrix=do.call(cbind, expression_matrixList)
  dim(expression_matrix)

  # Initialize the Seurat object with the raw (non-normalized data)
  seurat <- CreateSeuratObject(counts = expression_matrix, project =project, min.cells = min.cells, min.features = min.features)

  return(seurat)
}



#' @title cell.filter.seurat
#' @description A wrapper function to cell.filter.seurat
#' @usage cell.filter.seurat(seurat, sample.name, metrics_summary, more_nFeature_RNA = 200, Less_nFeature_RNA = 8000, percent.mt=10)
#' @param seurat Seurat object
#' @param sample.name single cell RNA sequening sample name
#' @param metrics_summary summary metrics
#' @param more_nFeature_RNA High cutoffs for filtering cells that have unique feature counts (default is 200)
#' @param Less_nFeature_RNA low cutoffs for filtering cells that have unique feature counts (default is 2500)
#' @param percent.mt low cutoffs for filtering cells that have >n percent mitochondrial counts (default is 5)
#' @details Creates a Seurat object containing only a subset of the cells in the original object.
#' @return a Seurat object containing only the relevant subset of cells
#' @export
cell.filter.seurat<- function(seurat,
                              sample.name,
                              metrics_summary,
                              more_nFeature_RNA = 200,
                              Less_nFeature_RNA = 8000,
                              percent.mt = 10){
  # sample name
  tbl=table(sub(".*-", "", rownames(seurat@meta.data)))
  tbl=tbl[match(1:length(tbl), names(tbl))]
  match(1:length(tbl), sub(".*-", "", rownames(seurat@meta.data)))
  seurat@meta.data$sample.name=unlist(lapply(1:length(sample.name), function(i) rep(sample.name[i], tbl[i])))
  table(seurat@meta.data$sample.name)

  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  seurat[["percent.mt"]] <- PercentageFeatureSet(object = seurat, pattern = "^MT-")
  tapply(seurat@meta.data$percent.mt, seurat@meta.data$sample.name, quantile)
  tapply(seurat@meta.data$nFeature_RNA, seurat@meta.data$sample.name, quantile)

  # summary qc metrix
  seurat@meta.data=data.frame(seurat@meta.data, metrics_summary[match(seurat@meta.data$sample.name, metrics_summary[,1]),setdiff(colnames(metrics_summary), colnames(seurat@meta.data))])

  ### filter cells
  seurat=seurat[,seurat$nFeature_RNA > more_nFeature_RNA]
  seurat=seurat[,seurat$nFeature_RNA < Less_nFeature_RNA]
  seurat=seurat[,seurat$percent.mt < percent.mt]

  return(seurat)
}




