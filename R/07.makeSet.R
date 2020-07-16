#' @title make.eset
#' @description For the expression data are transformed to a file with extension .eSet
#' @usage make.eset(expr, pdata=NULL, fdata=NULL, verbose=TRUE)
#' @param expr expression data
#' @param pdata phenotype data
#' @param fdata feature data
#' @param verbose default TRUE, Class to writing verbose messages to a connection or file.
#' @details store expression data in ExpressionSet format for convenient analysis
#' @return expression set
#' @export
make.eset <- function(expr, pdata=NULL, fdata=NULL, verbose=TRUE){
  Sys.setlocale(category = "LC_ALL", locale = "us")

  if (inherits(expr, "data.frame")) expr=as.matrix(expr)
  if (inherits(expr, "ExpressionSet")) expr=as.matrix(exprs(expr))

  #for null data
  if(is.null(pdata))    pdata=data.frame(phenoID=colnames(expr), row.names=colnames(expr))
  if(is.null(rownames(pdata))) rownames(pdata)=colnames(expr)

  if(is.null(fdata))  fdata=data.frame(featureID=rownames(expr), row.names=rownames(expr))
  if(is.null(rownames(fdata))) rownames(fdata)=rownames(expr)

  #match
  if(ncol(pdata)>1){
    p.st=match(colnames(expr), rownames(pdata))
    p.col=colnames(pdata)
    class(pdata)
    pdata=data.frame(pdata[p.st,p.col],row.names=rownames(pdata)[p.st])
    colnames(pdata)=p.col
  }


  f.st=match(rownames(expr), rownames(fdata))
  f.col=colnames(fdata)
  fdata=data.frame(fdata[f.st,], row.names=rownames(fdata)[f.st])
  colnames(fdata)=f.col

  #create eset
  metadata <- data.frame(labelDescription = colnames(pdata), row.names=colnames(pdata))
  phenoData<-new("AnnotatedDataFrame", data=as.data.frame(pdata), varMetadata=metadata)
  fmetadata <- data.frame(labelDescription = colnames(fdata), row.names=colnames(fdata))
  featureData<-new("AnnotatedDataFrame", data=as.data.frame(fdata), varMetadata=fmetadata)

  eset<-new("ExpressionSet", exprs=expr, phenoData=phenoData, featureData=featureData)
  if(verbose) print(eset)
  return(eset)
}


#' @title make.seurat.fdata
#' @description A function to make fdata using seurat object
#' @usage make.seurat.fdata(seurat, gene.ref.gtf, rda.dir)
#' @param seurat Seurat object
#' @param gene.ref.gtf gene reference gtf file
#' @param rda.dir rData directory
#' @details make feature data using seurat
#' @return feature data
#' @export
make.seurat.fdata <- function(seurat, gene.ref.gtf, rda.dir){

  message("[[",Sys.time(),"]] Run to make seurat.fdata --------")
  gtf=read.delim(gene.ref.gtf, comment.char = "#", header = F)
  gene.gtf=gtf[gtf$V3=="gene",]

  gene.annot=strsplit2(gene.gtf$V9, split = ";| ")
  gene.annot=data.frame("gene_id"=gene.annot[,2], "gene_name"=gene.annot[,8], "gene_source"=gene.annot[,11], "gene_biotype"=gene.annot[,14])
  head(gene.annot)


  mat.st=match(rownames(seurat), gene.annot$gene_name)

  rownames(seurat)[is.na(mat.st)]
  mat.st[is.na(mat.st)]=match(sub("\\.1", "", rownames(seurat)[is.na(mat.st)]), gene.annot$gene_name)
  sum(is.na(mat.st))

  fdata=data.frame(gene.gtf[mat.st,c(1,4,5,7)], gene.annot[mat.st,])

  colnames(fdata)[1:4]=c("chr", "str", "end", "strand")
  rownames(fdata)=rownames(seurat)

  save(fdata, file=file.path(rda.dir,"fdata.rda"))
  message("[[",Sys.time(),"]] Finish making seurat.fdata --------")

  return(fdata)
}

#' @title make.seurat.eset
#' @description A wrapper function to make.seurat.eset
#' @usage make.seurat.eset(seurat, slot=c("scale.data", "count.data", "count.data"), fdata, output.dir=NULL, save=T)
#' @param seurat Seurat object
#' @param slot assay data type of seurat object, c("scale.data", "count.data", "data")
#' @param fdata feature data
#' @param output.dir output directory
#' @param save whether save
#' @details make expression set using seurat
#' @return expression set
#' @export
make.seurat.eset <- function(seurat, slot=c("scale.data", "count.data", "count.data"), fdata, output.dir=NULL, save=T){

  message("[[",Sys.time(),"]] Run to make seurat.eset --------")
  expr=GetAssayData(seurat, slot = slot)
  seurat.eset=make.eset(expr = expr, pdata = seurat@meta.data, fdata = seurat[['RNA']]@meta.features, verbose = FALSE)
  sum(!rownames(fdata)==featureNames(seurat.eset))
  fData(seurat.eset)=fdata

  if(save) save(seurat.eset, file = file.path(output.dir, "seurat.eset.rda"))
  message("[[",Sys.time(),"]] Finish making seurat.eset --------")

  return(seurat.eset)
}

