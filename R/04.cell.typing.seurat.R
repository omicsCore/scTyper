

#' @title cell.typing.seurat
#' @description Wrapper function for cell typing of seurat object
#' @param seurat Seurat object
#' @param marker Cell markers to use in cell typing, character or List (identifier or StudyName or User defined gene list)
#' @param cell.typing.method cell typing method, c("NTP", "ES", "Average"), (default = "NTP")
#' @param wd working directory
#' @param slot assay data type of seurat object, c("scale.data", "count.data", "data")
#' @param assay Assay of seurat object
#' @param ntp.dir Output directory of NTP
#' @param rda.dir Path of the RData saving directory
#' @param NTP.g.filter.method Method of gene filtering in NTP c(sd (Default), mad, none)
#' @param NTP.gene.filter.cutoff Cut-off score of standard deviation in NTP
#' @param NTP.distance Method of calculating distance in NTP, either c("correlation" or "cosine").
#' @param NTP.norm.method Method of normalization in NTP, either c("none", "row.std")
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @export
cell.typing.seurat=function(seurat,
                            marker="Puram.2017.HNSCC",
                            cell.typing.method=c("NTP", "ES", "Average"),
                            level=c("cell", "cluster"),
                            wd,
                            slot=c("scale.data", "count.data", "data"),
                            assay='RNA',
                            ntp.dir,
                            rda.dir,
                            NTP.g.filter.method=c("sd", "mad", "none"),
                            NTP.gene.filter.cutoff=0.3,
                            NTP.distance=c("cosine","correlation"),
                            NTP.norm.method=c("none", "row.std"),
                            mc.cores=1){
  level=match.arg(level)
  expr=GetAssayData(seurat, slot = slot)
  markerList=get.markerList(marker)
  markerList=lapply(markerList, function(a) intersect(a, rownames(seurat)))

  if ((cell.typing.method!="Average") & ("1" %in% as.character(lapply(markerList, function(l) length(l)))) ==TRUE ){
    stop("### Error :: If there is only one marker gene for a particular cell type, please use the average method ###")
  }
  #select cell typing method
  if(cell.typing.method=="NTP"){
    seurat=cell_type_NTP(seurat, wd, markerList, slot=slot, assay=assay, output.dir=ntp.dir, rda.dir, NTP.g.filter.method=NTP.g.filter.method, NTP.gene.filter.cutoff=NTP.gene.filter.cutoff, NTP.distance=NTP.distance, NTP.norm.method=NTP.norm.method, mc.cores=mc.cores)
    cell.type=as.character(seurat$ntp.fdr)

  }else if(cell.typing.method=="ES"){
    message("[[",Sys.time(),"]] Run Cell typing using Enrichment Score --------")
    marker.ES=preRanked.GSEA(expr, SIGDB = markerList, mc.cores=mc.cores)

    cell.type=as.character(apply(t(marker.ES), 1, function(a) names(sort(-a)[1])) )
    seurat$es=cell.type
  }else if(cell.typing.method=="Average"){
    message("[[",Sys.time(),"]] Run Cell typing using average expression value--------")
    marker.mean=sapply(markerList, function(a) as.matrix(if(length(unlist(a))==1){expr[a,]}else{colMeans(expr[a,], na.rm=T)}))
    cell.type=apply(marker.mean, 1, function(a) names(sort(-a)[1]))
    seurat$average=cell.type
  }

  #select cell typing level
  if(level=="cluster"){
    df=table(cell.type, seurat$seurat_clusters)
    rownames(df)
    c.name=apply(df,2, function(a) rownames(df)[order(-a)[1]])
    c.cell.type=sapply(seurat$seurat_clusters, function(a) c.name[a])
    seurat$cell.type=as.factor(c.cell.type)
  }else if(level=="cell"){ #cell level
    seurat$cell.type=as.factor(cell.type)
  }else{
    stop("### You have to select parameter as either cell or cluster ###")
  }

  save(seurat, file = file.path(rda.dir, "seurat.rda"))
  return(seurat)
}


