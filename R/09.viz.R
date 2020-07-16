
#' @title draw.heatmap
#' @description Visualize heatmap of scTyper
#' @param seurat seurat object
#' @param wd working directory
#' @param run.inferCNV whether run inferCNV (default=TRUE)
#' @param slot assay data type of seurat object, c("scale.data", "count.data", "data")
#' @param marker cell type marker
#' @export
draw.heatmap=function(seurat,
                      wd,
                      run.inferCNV=TRUE,
                      slot=c("scale.data", "count.data", "data"),
                      marker="Puram.2017.HNSCC.TME"){

  slot=match.arg(slot)
  markerList=get.markerList(marker)
  markerList=lapply(markerList, function(a) intersect(a, rownames(seurat)))
  g = unlist(markerList)
  expr=GetAssayData(seurat, slot = slot)

  if(run.inferCNV==TRUE){
    o.st=order(seurat$cell.type, seurat$cnv.score, colMeans(expr[g,colnames(seurat)]))
  }else{ o.st=order(seurat$cell.type, colMeans(expr[g,colnames(seurat)]))}

  cell_type = seurat@meta.data$cell.type[o.st]



  annotation.bar=data.frame(cell_type)

  if(nrow(expr)>20){show_rownames=F}else{show_rownames=T}
  if(ncol(expr)>50){show_colnames=F}else{show_colnames=T}
  cluster_rows=F
  cluster_cols=F
  if(nrow(expr)==1){expr=rbind(expr,expr)}
  if(!is.null(annotation.bar)){rownames(annotation.bar)=colnames(expr)}
  ha = HeatmapAnnotation(df = annotation.bar)


  Heatmap(expr[g,o.st],
          name = "expression",
          col = colorRamp2(c(-0.5, 0, 0.5), c("SteelBlue", "white", "red")),
          top_annotation = ha,
          show_row_names = show_rownames,
          show_column_names = show_colnames,
          cluster_rows = cluster_rows,
          cluster_columns = cluster_cols)

}
#' @title cnv.distribution
#' @description Visualize cnv distribution of scTyper
#' @param seurat seurat object
#' @param wd working directory
#' @param marker cell type marker
#' @param slot assay data type of seurat object, c("scale.data", "count.data", "data")
#' @export
cnv.distribution <- function(seurat,
                             wd,
                             marker="Puram.2017.HNSCC.TME",
                             slot=c("scale.data", "count.data", "data")){
  scale.dat=GetAssayData(seurat, slot = slot)

  markerList=get.markerList(marker=marker)
  g = unlist(markerList)

  o.st=order(seurat$cell.type, seurat$cnv.score, colMeans(scale.dat[g,colnames(seurat)]))

  mal.st=seurat$malignant.st
  barplot(seurat$cnv.score[o.st], border=NA, ylim=c(0,500),col = ifelse(mal.st[o.st]==TRUE,"red","grey"), names.arg = NA)
}


