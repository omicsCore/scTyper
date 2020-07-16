
#' @title run.inferCNV
#' @description A wrapper function to run inferCNV
#' @param seurat Seurat object
#' @param assay Name of assay to pull data from seurat object
#' @param output.dir output directory
#' @param rda.dir rData directory
#' @param fdata feature information
#' @param pheno_info phenotype information
#' @param feature.to.test featuest to test either "tissue.type" or "cell.type"
#' @param cells.test_reference a vector containing the classifications of the reference (normal) cells to use for infering cnv
#' @param cells.test_excluded cell type to exclude functional enrichment analysis
#' @param fc.cutoff fold change cutoff
#' @param cutoff.gene.cluster A cutoff P-value for filtering out the gene clusters (calculated from GO analysis)
#' @param min_mean_expr_cutoff the minimum mean value allowed for a gene to be retained in the expression matrix.
#' @param window_length length of window (number of genes) for the moving average
#' @param smooth_ends perform smoothing at the ends of the chromosomes (default:TRUE)
#' @param recenter_method method to select the center of the cell expression value. (default:'mean', option:'mean', 'median')
#' @param ordered order bool, default FALSE. Sort descending vs. ascending
#' @param inv_log mean values will be determined based on (2^x -1)
#' @param sd_amplifier multiplicative factor applied to the standard deviation to alter the noise range (default: 1.5)
#' @param bp base pair
#' @param sd.cut standard deviation cutoff
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details CNV inference
#' @return Seurat object
#' @export
run.inferCNV <- function(seurat,
                         assay='RNA',
                         output.dir = "./",
                         rda.dir = "./data",
                         fdata,
                         pheno_info = pheno.df,
                         feature.to.test = c("tissue.type", "cell.type"),
                         cells.test_reference = "Normal",
                         cells.test_excluded = c("Epithelial"),
                         fc.cutoff = 0.05,
                         cutoff.gene.cluster=0.05,

                         min_mean_expr_cutoff = 0.1,
                         window_length = 101,
                         smooth_ends = TRUE,
                         recenter_method = "median",
                         ordered = FALSE,
                         inv_log = TRUE,
                         sd_amplifier = 1.5,
                         bp = 1000000,
                         sd.cut = 0.3,
                         mc.cores=1){

  dir.create(output.dir, showWarnings = F)
  if(feature.to.test=="tissue.type"){
    seurat@meta.data$cell.group = pheno_info[match(seurat@meta.data$sample.name, pheno_info[,"Sample_ID"]),"TissueType"]
  }else if(feature.to.test=="cell.type"){
    if(is.null(seurat@meta.data$cell.type.ref)){
      stop("### You have to assign separately cell type reference to seurat@meta.data$cell.type.ref. ###")
    }
    seurat@meta.data$cell.group = seurat@meta.data$cell.type.ref
  }


  message("[[",Sys.time(),"]] Prepare inferCNV input data --------")
  ## make input data
  # 1. raw.count matrix
  raw.mat=as.matrix(seurat[[assay]]@counts)
  dim(raw.mat)
  # 2. cell annotation
  annot=data.frame(cell.names=colnames(seurat), group.st=seurat@meta.data$cell.group)
  dim(annot)
  # 3. gene feature data
  head(fdata)

  seurat[[assay]]@meta.features <- data.frame(seurat[[assay]]@meta.features, fdata)
  names(seurat[[assay]]@meta.features)

  gene.od=data.frame(gene=rownames(seurat[[assay]]@meta.features), seurat[[assay]]@meta.features[,c("chr", "str", "end")])


  table(gene.od$chr)
  fil.st=is.element(gene.od$chr, c(1:22, "X"))

  sum(is.na(match(gene.od$gene, rownames(seurat))))

  raw.mat=raw.mat[fil.st,]
  dim(raw.mat)
  gene.od=gene.od[fil.st,]
  dim(gene.od)
  table(as.character(gene.od$chr))
  gene.od$chr=factor(as.character(gene.od$chr), levels=c(1:22,"X"))
  od.st=order(gene.od$chr, gene.od$str)
  match(c(1:22, "X"), gene.od[od.st,]$chr)
  sum(is.na(match(gene.od$gene, rownames(seurat))))

  #make directory
  dir.create(output.dir, showWarnings=FALSE)

  # write files
  raw.mat.fn=file.path(output.dir, "raw_counts_matrix.txt")
  annot.fn=file.path(output.dir, "annotations.txt")
  gene.od.fn=file.path(output.dir, "gene_order.txt")

  write.table(raw.mat[od.st,], raw.mat.fn, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(annot, annot.fn, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(gene.od[od.st,], gene.od.fn, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  ##run inferCNV
  message("[[",Sys.time(),"]] Create infercnv object --------")
  # 1. create infercnv obj
  infercnv_obj = infercnv::CreateInfercnvObject(raw_counts_matrix=raw.mat[od.st,],
                                                annotations_file=annot.fn,
                                                gene_order_file=gene.od.fn,
                                                ref_group_names=cells.test_reference)

  message("[[",Sys.time(),"]] Filter and normalize infercnv object --------")
  # 2. filter out low expressed genes
  infercnv_obj <- infercnv:::require_above_min_mean_expr_cutoff(infercnv_obj, min_mean_expr_cutoff=min_mean_expr_cutoff)

  # 3. Normalize each cell’s counts for sequencing depth
  infercnv_obj <- infercnv:::normalize_counts_by_seq_depth(infercnv_obj)

  # 4. perform Anscombe normalization
  infercnv_obj <- infercnv:::anscombe_transform(infercnv_obj)

  # 5. log transform the normalized counts:
  infercnv_obj <- infercnv:::log2xplus1(infercnv_obj)

  # 6. Apply maximum bounds to the expression data to reduce outlier effects
  infercnv_obj <- infercnv:::apply_max_threshold_bounds(infercnv_obj, threshold=mean(abs(infercnv:::get_average_bounds(infercnv_obj))) )

  message("[[",Sys.time(),"]] Perform smoothing across chromosomes --------")

  # 7. perform smoothing across chromosomes
  infercnv_obj <- infercnv:::smooth_by_chromosome(infercnv_obj, window_length=window_length, smooth_ends=smooth_ends)

  message("[[",Sys.time(),"]] Finish smoothing across chromosomes --------")

  # 8. re-center each cell
  infercnv_obj <- infercnv:::center_cell_expr_across_chromosome(infercnv_obj, method = recenter_method)

  ### 9. cell type specific CNV gene filter
  cset=infercnv2cset(infercnv_obj = infercnv_obj, pdata = seurat@meta.data)
  save(cset, file=file.path(rda.dir, 'cset.rda'))

  message("[[",Sys.time(),"]] Run sample permutation t-test  --------")

  cell.type.set=setdiff(levels(cset$cell.type), cells.test_excluded)

  perm.t.resList = perm.subcset.t(cset, cell.type.set, rda.dir, ordered= ordered, mc.cores=mc.cores)

  message("[[",Sys.time(),"]] Run cell type specific GO analysis  --------")

  cts.geneClustList = cts.geneSetCluster(cset, rda.dir = rda.dir, perm.t.resList, fc.cutoff = fc.cutoff, bp = bp)
  gprofiler.resList = cts.GO(cell.type.set, rda.dir = rda.dir, cts.geneClustList)

  infercnv_obj =fil.infercnv_obj(infercnv_obj, cset, rda.dir = rda.dir, gprofiler.resList, cell.type.set, cts.geneClustList, cutoff.gene.cluster=cutoff.gene.cluster)

  # 10. subtract the reference values from observations, now have log(fold change) values
  infercnv_obj <- infercnv:::subtract_ref_expr_from_obs(infercnv_obj, inv_log=inv_log)

  # 11. invert log values
  infercnv_obj <- infercnv:::invert_log2(infercnv_obj)

  # 12. Removing noise
  infercnv_obj <- infercnv:::clear_noise_via_ref_mean_sd(infercnv_obj, sd_amplifier = sd_amplifier)

  # 13. Remove outlier data points
  infercnv_obj = infercnv:::remove_outliers_norm(infercnv_obj)

  message("[[",Sys.time(),"]] Filter out by standard deviations in reference groups  --------")
  ### 14. Filter out by standard deviations in reference groups
  ref.sd=apply(infercnv_obj@expr.data[,infercnv_obj@reference_grouped_cell_indices[[1]]], 1, sd)

  save(infercnv_obj, file=file.path(rda.dir, 'infercnv_obj.o.rda'))

  sd.fil.st <- ref.sd < sd.cut
  infercnv_obj@expr.data=infercnv_obj@expr.data[sd.fil.st,]
  infercnv_obj@count.data=infercnv_obj@count.data[sd.fil.st,]
  infercnv_obj@gene_order=infercnv_obj@gene_order[sd.fil.st,]

  save(infercnv_obj, file=file.path(rda.dir, 'infercnv_obj.sd_filtered.rda'))
  save(infercnv_obj, file=file.path(rda.dir, 'infercnv_obj.rda'))

  message("[[",Sys.time(),"]] Make filterd cnv set  --------")
  fil.cset=infercnv2cset(infercnv_obj = infercnv_obj, pdata = seurat@meta.data)
  fil.cset

  # cnv score
  pData(fil.cset)$cnv.score=apply(exprs(fil.cset)-1, 2, function(a) sum(abs(a)))

  seurat$cnv.score=pData(fil.cset)$cnv.score
  save(fil.cset, file=file.path(rda.dir, 'fil.cset.rda'))

  return(seurat)
}


#' @title infercnv2cset
#' @description A wrapper function to infercnv2cset
#' @usage infercnv2cset(infercnv_obj, pdata)
#' @param infercnv_obj an infercnv object
#' @param pdata phenotype data
#' @details Creation of cset using infercnv object.
#' @return cset
#' @export
infercnv2cset=function(infercnv_obj, pdata){
  ### 9. cell type specific CNV gene filter
  cset=make.eset(expr = infercnv_obj@expr.data, pdata = pdata[colnames(infercnv_obj@expr.data),], fdata = infercnv_obj@gene_order)
  cset
}


#' @title perm.subcset.t
#' @description A wrapper function to perm.subcset.t
#' @usage perm.subcset.t(cset, cell.type.set, rda.dir, ordered= FALSE, levels= c(1,0), mc.cores=5)
#' @param cset cnv Set
#' @param cell.type.set cell type set
#' @param rda.dir rData directory
#' @param ordered order bool, default FALSE. Sort descending vs. ascending
#' @param levels confidence level of the interval.
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details Performs t-tests on cnv subset
#' @return list of t-test result
#' @import parallel
#' @export
perm.subcset.t <- function(cset, cell.type.set, rda.dir ,ordered = FALSE, levels= c(1,0), mc.cores=1){
  #make subcset
  sub.cset=cset[,is.element(cset$cell.type, cell.type.set)]
  sub.cset
  # permutation t-test for cell-type vs. others
  perm.t.resList=lapply(cell.type.set, function(a){
    perm.ttest(sub.cset, g.st = factor(as.numeric(is.element(sub.cset$cell.type, a)), levels = levels), mc.cores = mc.cores, ordered = ordered)
  })

  names(perm.t.resList)=cell.type.set
  save(perm.t.resList, file = file.path(rda.dir, "perm.t.resList.rda"))

  return(perm.t.resList)

}


#' @title get.geneClust
#' @description A wrapper function to get.geneClust
#' @usage get.geneClust(gr, bp=1000000)
#' @param gr GRanges object
#' @param bp base pair
#' @details get gene cluster
#' @return gene clust List
#' @export
get.geneClust=function(gr, bp=1000000){
  seqnames=names(which(table(as.character(seqnames(gr)))>1))

  geneClustList=lapply(seqnames, function(chr){
    sub.gr=gr[as.character(seqnames(gr))==chr]
    dist=0; for(i in 2:length(sub.gr)){ dist=c(dist, start(sub.gr)[i]-end(sub.gr)[i-1]) }
    sub.gr$dist=dist

    st=which(sub.gr$dist>bp)
    if(length(st)==0) g=list(names(sub.gr))
    else {
      g=NULL
      for(i in 1:length(st)){
        if(i==1) g=list(names(sub.gr)[1:(st[i]-1)])
        else if(i<length(st)) g=c(g, list(names(sub.gr)[st[i]:(st[i+1]-1)]))

        if(i==length(st)) g=c(g, list(names(sub.gr)[st[i]:length(sub.gr)]))
      }
    }
    return(g)
  })

  names(geneClustList)=seqnames
  return(geneClustList)
}

#' @title cts.geneSetCluster
#' @description A wrapper function to cts.geneSetCluster
#' @usage cts.geneSetCluster(cset, rda.dir, perm.t.resList, fc=0.05, bp = 1000000)
#' @param cset cnv Set
#' @param rda.dir rData directory
#' @param perm.t.resList list of t-test result
#' @param fc fold change
#' @param bp base pair
#' @details get gene cluster
#' @return cell type specific geneClustList
#' @export
cts.geneSetCluster <- function(cset, rda.dir, perm.t.resList, fc.cutoff=0.05, bp = 1000000){
  # identify cell type specific geneSetCluster
  cts.geneList=lapply(perm.t.resList, function(perm.t.res){
    list(up=rownames(perm.t.res)[which(perm.t.res$fc > fc.cutoff & perm.t.res$perm.p < 10^-10)],
         dn=rownames(perm.t.res)[which(perm.t.res$fc < -fc.cutoff & perm.t.res$perm.p < 10^-10)])
  })

  fData(cset)$strand="*"
  cset.gr=df2gr(df = fData(cset), seqnames = "chr", start = "start", end = "stop", strand = "strand")
  names(cset.gr)=cset.gr$symbol=featureNames(cset)
  range(width(cset.gr))

  cts.geneClustList=lapply(cts.geneList, function(cts.gene) c(get.geneClust(gr = cset.gr[cts.gene$up], bp = bp),
                                                              get.geneClust(gr = cset.gr[cts.gene$dn], bp = bp)))

  save(cts.geneClustList, file = file.path(rda.dir, "cts.geneClustList.rda"))

  return(cts.geneClustList)
}


#' @title cts.GO
#' @description A wrapper function to cts.GO
#' @usage cts.GO(cell.type.set, rda.dir,cts.geneClustList)
#' @param cell.type.set cell type set
#' @param rda.dir rData directory
#' @param cts.geneClustList cell type specific geneClustList
#' @details Interface to the g:Profiler tool for finding enrichments in gene lists.
#' @return gprofiler result List
#' @export
cts.GO <- function( cell.type.set, rda.dir, cts.geneClustList){
  gprofiler.resList=mclapply(cell.type.set, function(cell.type) lapply(cts.geneClustList[[cell.type]], function(a) lapply(a, function(b) gprofiler(b, src_filter = c("GO:BP", "KEGG")))))

  names(gprofiler.resList)=cell.type.set
  save(gprofiler.resList, file = file.path(rda.dir, "gprofiler.resList.rda"))

  return(gprofiler.resList)
}


#' @title fil.infercnv_obj
#' @description A wrapper function to fil.infercnv_obj
#' @usage fil.infercnv_obj(infercnv_obj, cset, rda.dir, gprofiler.resList, cell.type.set, cts.geneClustList)
#' @param infercnv_obj infercnv object
#' @param cset cnv Set
#' @param rda.dir rData directory
#' @param gprofiler.resList gprofiler result List
#' @param cell.type.set cell type set
#' @param cts.geneClustList cell type specific geneClustList
#' @details get gene cluster
#' @return infercnv object
#' @export
fil.infercnv_obj<- function(infercnv_obj,
                            cset,
                            rda.dir,
                            gprofiler.resList,
                            cell.type.set,
                            cts.geneClustList,
                            cutoff.gene.cluster=0.05){


  go.df=lapply(cell.type.set, function(cell.type){
    chrs=names(gprofiler.resList[[cell.type]])
    df=lapply(chrs, function(chr){
      len=length(gprofiler.resList[[cell.type]][[chr]])
      dfList=lapply(1:len, function(i) {
        if(nrow(gprofiler.resList[[cell.type]][[chr]][[i]])>0) df=data.frame(gprofiler.resList[[cell.type]][[chr]][[i]], chr=chr, idx=i) else df=NULL
      })
      do.call(rbind, dfList)
    })

    if(!is.null(do.call(rbind, df))) return(data.frame(do.call(rbind, df), cell.type=cell.type))
    else return(NULL)
  })
  go.df=do.call(rbind, go.df)

  fil.df=unique(data.frame(chr=go.df$chr, idx=go.df$idx, cell.type=go.df$cell.type, p.value=go.df$p.value))
  fil.df=fil.df[(fil.df$p.value<cutoff.gene.cluster),]
  dim(fil.df)
  fil.df=do.call(rbind, lapply(1:nrow(fil.df), function(i) data.frame(fil.df[i,], gene=cts.geneClustList[[as.character(fil.df$cell.type[i])]][[as.character(fil.df$chr)[i]]][[fil.df$idx[i]]] )))

  cts.cnv.df=fil.df
  save(cts.cnv.df, file = file.path(rda.dir, "cts.cnv.df.rda"))
  head(cts.cnv.df)

  head(cts.cnv.df)
  fil.st=!is.element(featureNames(cset), unique(fil.df$gene))
  fil.cset=cset[fil.st,]
  fil.cset

  sum(!featureNames(cset)==rownames(infercnv_obj@expr.data))
  infercnv_obj@expr.data=infercnv_obj@expr.data[fil.st,]
  infercnv_obj@count.data=infercnv_obj@count.data[fil.st,]
  infercnv_obj@gene_order=infercnv_obj@gene_order[fil.st,]
  save(infercnv_obj, file=file.path(rda.dir, 'infercnv_obj.cell_type_specific_csv_filtered.rda'))

  return(infercnv_obj)
}

