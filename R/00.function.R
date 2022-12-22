#' @title update.sig.db
#' @description Update sig.db
#' @usage update.sigTyper.db(sig.db.path, db.name=c("CellMarker", "sigTyper.db"), output.dir=system.file("/data",package = "scTyper") )
#' @param sig.db.path Path of sig.db.txt
#' @param db.name database name to update. either c("sigTyper.db", "CellMarker")
#' @param output.dir storage path of sig.db marker
#' @export
update.sig.db <- function(sig.db.path, db.name=c("sigTyper.db", "CellMarker"), output.dir=system.file("/data",package = "scTyper") ){
  sig.db <- read.delim(sig.db.path)
  names(sig.db)
  lst=rep(list(NULL), length(unique(sig.db$StudyName)))
  names(lst)=unique(sig.db$StudyName)
  for(i in 1:nrow(sig.db)){
    varname=sig.db$StudyName[i]
    ind=which(names(lst)==varname)
    cell.name=as.character(sig.db[i,]$CellName)
    symbol=as.character(sig.db[i,]$GeneSymbol)
    sigList=lapply(symbol, function(a) unlist(strsplit(a, split = ", ")))
    names(sigList)=cell.name
    lst[[ind]]=c(lst[[ind]], (sigList))
  }

  if(db.name=="CellMarker"){
    CellMarker=lst
    save(CellMarker, file = file.path(output.dir,'CellMarker.rda') )
  }else if(db.name=="sigTyper.db"){
    sigTyper.db=lst
    save(sigTyper.db, file = file.path(output.dir,'sigTyper.db.rda') )
  }
}


#' @title get.markerList
#' @usage get.markerList(marker)
#' @description get markerList from scTyper database
#' @param marker Signature_list or Signature name of scTyper db or User-defined list of marker genes
#' @return marker list
#' @export
get.markerList=function(marker="Puram.2017.HNSCC.TME"){
  if(class(marker)=="character"){
    CellMarker=scTyper::CellMarker
    sigTyper.db=scTyper::sigTyper.db
    markerDB=c(CellMarker, sigTyper.db)

    markerList=list()
    for(i in 1:length(marker)){
      varname=strsplit(marker, split=":")[[i]][1]
      one_markerlist=markerDB[[varname]]
      if(!is.na(strsplit(marker, split=":")[[i]][3])) {one_markerlist=one_markerlist[strsplit(marker, split=":")[[i]][3]]}
      markerList=c(markerList, one_markerlist)
    }

    ##remove overlapped gene
    if(any(as.vector(table(as.character(unlist(markerList))))>1)==TRUE){
      df=data.frame(table(as.character(unlist(markerList))))
      df_freq1gene=df[df$Freq==1,]
      lstdiff_gene=df_freq1gene$Var1
      markerList.o=markerList
      markerList=sapply(1:length(markerList), function(i) markerList[[i]]=intersect(markerList.o[[i]], lstdiff_gene) )

      names(markerList)=names(markerList.o)
    }

    markerList_name=unique(names(markerList))
    for(i in 1:length(unique(names(markerList)))){
      markerList_markers=unique(as.character(unlist(markerList[names(markerList) == markerList_name[i]])))
      markerList[names(markerList) == markerList_name[i]]=NULL
      markerList[[markerList_name[i]]]=markerList_markers
    }

  }else if(class(marker)=="list"){
    markerList=marker
  }else{
    stop("### Error :: Please write upper case letter and Class of marker is not list or character###")
  }
  return(markerList)
}

#' @title list2matrix
#' @usage list2matrix(List)
#' @description convert list to matrix
#' @param List list
#' @export
list2matrix<-function(List){
  maxn  =max(unlist(lapply(List, length)))
  for (i in 1:length(List))  length(List[[i]])=maxn
  res=NULL
  for (i in 1:length(List)){ res=cbind(res, as.character(List[[i]]))}
  colnames(res)=names(List)
  res[which(is.na(res))]=""
  return(res)
}


#' @title preRanked.GSEA
#' @usage preRanked.GSEA(expr, SIGDB, weighted.score.type = 0, correl.vector = NULL, n.cutoff=1,mc.cores=1)
#' @description a universal gene set enrichment analysis tools
#' @param expr expression metrix
#' @param SIGDB gene signature list
#' @param weighted.score.type Type of weight score
#' @param correl.vector correlation vector
#' @param n.cutoff number of cutoff
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @export
preRanked.GSEA<-function (expr, SIGDB, weighted.score.type = 0, correl.vector = NULL, n.cutoff=1, mc.cores=1){
  #preprocess(expr, method="c", cen.method="g")
  geneList = lapply(SIGDB, function(a) intersect(a, rownames(expr)))
  gene.set.order = lapply(geneList, function(a) match(a, rownames(expr)))
  gene.set.order=gene.set.order[lapply(gene.set.order, length)>=n.cutoff]
  ES = mclapply(as.list(1:ncol(expr)), function(a) as.numeric(lapply(gene.set.order,  function(b) as.numeric(GSEA.EnrichmentScore2(order(-expr[, a]), b, weighted.score.type = weighted.score.type, correl.vector = correl.vector)))), mc.cores=mc.cores)
  ES = apply(list2matrix(ES), 2, as.numeric)
  rownames(ES) = names(gene.set.order)
  colnames(ES) = colnames(expr)
  return(ES)
}


#' @title GSEA.EnrichmentScore2
#' @usage GSEA.EnrichmentScore2(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL)
#' @description Run GSEA(Gene Set Enrichment Analysis)
#' @param gene.list gene signature list
#' @param gene.set gene set
#' @param weighted.score.type Type of weighted score
#' @param correl.vector correlation vector
#' @export
GSEA.EnrichmentScore2<-function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {
  N <- length(gene.list)
  Nh <- length(gene.set)
  Nm <- N - Nh
  loc.vector <- vector(length = N, mode = "numeric")
  peak.res.vector <- vector(length = Nh, mode = "numeric")
  valley.res.vector <- vector(length = Nh, mode = "numeric")
  tag.correl.vector <- vector(length = Nh, mode = "numeric")
  tag.diff.vector <- vector(length = Nh, mode = "numeric")
  tag.loc.vector <- vector(length = Nh, mode = "numeric")
  loc.vector[gene.list] <- seq(1, N)
  tag.loc.vector <- loc.vector[gene.set]
  tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
  if (weighted.score.type == 0) {
    tag.correl.vector <- rep(1, Nh)
  }
  else if (weighted.score.type == 1) {
    tag.correl.vector <- correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  }
  else if (weighted.score.type == 2) {
    tag.correl.vector <- correl.vector[tag.loc.vector] *
      correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  }
  else {
    tag.correl.vector <- correl.vector[tag.loc.vector]^weighted.score.type
    tag.correl.vector <- abs(tag.correl.vector)
  }
  norm.tag <- 1/sum(tag.correl.vector)
  tag.correl.vector <- tag.correl.vector * norm.tag
  norm.no.tag <- 1/Nm
  tag.diff.vector[1] <- (tag.loc.vector[1] - 1)
  tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh -
                                                                      1)] - 1
  tag.diff.vector <- tag.diff.vector * norm.no.tag
  peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
  valley.res.vector <- peak.res.vector - tag.correl.vector
  max.ES <- max(peak.res.vector)
  min.ES <- min(valley.res.vector)
  ES <- signif(ifelse(max.ES > -min.ES, max.ES, min.ES), digits = 5)
  return(list(ES = ES))
}

#' @title perm.ttest
#' @description Run permutation T-test
#' @usage perm.ttest(eset, g.st,  level=NULL, t.test=F, permp=T, permp.exact=NULL, ordered=T, mc.cores=1,...)
#' @param eset expression set
#' @param g.st group subset
#' @param levels confidence level of the interval.
#' @param permp bool, defalut TRUE. Calculating permuted T test p-values or not
#' @param ordered order bool, default FALSE. Sort descending vs. ascending
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details Performs sample t-tests on vectors of data.
#' @return t test result
#' @import parallel
#' @export
perm.ttest =function(eset, g.st,  level=NULL, t.test=FALSE, permp=TRUE, ordered=T, mc.cores=1,...){
  if(inherits(eset, "ExpressionSet"))   expr=exprs(eset) else    expr=eset

  if(ncol(expr)!=length(g.st))  cat("class labels has a different length")
  if(!inherits(g.st, "factor")) g.st=factor(g.st)
  if(!is.null(level)) g.st=factor(g.st, level=level)
  res=NULL
  if(t.test){
    message("Calculating T test p-values")
    if(mc.cores>1){
      if(Sys.info()[['sysname']]=="Windows"){
        res=mclapply(1:nrow(expr), function(a) try(t.test(as.numeric(expr[a,])~g.st), silent=T),mc.cores=mc.cores,expr=expr,g.st=g.st, packageToLoad=c("stat","perm"))
      }else{
        res=mclapply(1:nrow(expr), function(a) try(t.test(as.numeric(expr[a,])~g.st), silent=T),mc.cores=mc.cores,...=...)
      }
    }
    if(mc.cores==1) res=lapply(1:nrow(expr), function(a) try(t.test(as.numeric(expr[a,])~g.st), silent=T))

    tval=as.numeric(sapply(res, function(a) try(a$stat, silent=T)))
    test.p=as.numeric(sapply(res, function(a) try(a$p.val, silent=T)))
    res=data.frame(t.stat=(tval), ttest.p=test.p)
    rownames(res)=rownames(expr)
  }
  if(permp)  {

    message("Calculating permuted T test p-values")
    if(mc.cores>1){
      if(Sys.info()[['sysname']]=="Windows") {
        res$perm.p=as.numeric(mclapply(1:nrow(expr), function(a) try(permTS(as.numeric(expr[a,]) ~ g.st)$p.value, silent=T),mc.cores=mc.cores,cluster.export=F, expr=expr,g.st=g.st, packageToLoad=c("stat","perm")))
      }else{
        res$perm.p=as.numeric(parallel::mclapply(1:nrow(expr), function(a) try(permTS(as.numeric(expr[a,]) ~ g.st)$p.value, silent=T),mc.cores=mc.cores))
      }
    }
    if(mc.cores==1) res$perm.p=as.numeric(lapply(1:nrow(expr), function(a) try(permTS(as.numeric(expr[a,]) ~ g.st)$p.value, silent=T)))
    res$FDR = p.adjust(res$perm.p, "BH")
  }

  class.mean=sapply(levels(g.st), function(a) rowMeans(expr[,which(g.st==a)], na.rm=T))
  colnames(class.mean) = paste(colnames(class.mean), "(mean)")
  fc=as.matrix(class.mean[,1]-class.mean[,2])
  res=cbind(as.data.frame(res), class.mean, fc)

  if(ordered)  res=res[order(-fc),]
  return(res)
}

#' @title df2gr
#' @description A wrapper function to make dataframe to GRange
#' @usage df2gr(df, seqnames, start, end, strand)
#' @param df dataframe
#' @param seqnames A character vector of recognized names for the column in df that contains the chromosome name (sequence name) associated with each genomic range.
#' @param start A character vector of recognized names for the column in df that contains the start positions of the genomic ranges.
#' @param end A character vector of recognized names for the column in df that contains the end positions of the genomic ranges.
#' @param strand A character vector of recognized names for the column in df that contains the strand associated with each genomic range.
#' @details the workhorse behind the coercion method from data.frame to GRanges.
#' @return GRanges object
#' @export
df2gr=function(df, seqnames, start, end, strand){
  meta.st=setdiff(names(df),c(seqnames, start, end,strand))
  start=as.numeric(df[,start])
  end=as.numeric(df[,end])
  strand=as.character(df[,strand])
  seqnames=update.seqnames(df[,seqnames])
  if(invalid(strand)) strand=NULL
  GRanges(seqname=seqnames,strand, ranges=IRanges(start,end), df[,meta.st])
}


#' @title update.seqnames
#' @description updata the column in df that contains the chromosome name (sequence name)
#' @usage update.seqnames(seqnames)
#' @param seqnames A character vector of recognized names for the column in df that contains the chromosome name (sequence name) associated with each genomic range.
#' @return sequence name
#' @export
update.seqnames=function(seqnames){
  seqnames=as.character(seqnames)
  seqnames=sub("chr", "",seqnames, ignore.case=T)
  seqnames=sub("23", "X", seqnames, ignore.case=T)
  seqnames=sub("24", "Y", seqnames, ignore.case=T)
  seqnames=paste0("chr", seqnames)
  seqnames
}


#' @title invalids
#' @description Test if a value is missing, empty, or contains only NA or NULL values.
#' @usage invalid(x)
#' @param x value to be tested
#' @return Bool
#' @export
invalid= function (x)
{
  if (missing(x) || is.null(x) || length(x) == 0)
    return(TRUE)
  if (is.list(x))
    return(all(sapply(x, invalid)))
  else if (is.vector(x))
    return(all(is.na(x)))
  else return(FALSE)
}

#' @title loadTestData
#' @description load test data in 'extdata' directory
#' @usage loadTestData()
#' @return test.seurat
#' @export
loadTestData = function(){
  file_name = system.file('extdata/GSE103322.seurat.reduced.rda', package = 'scTyper')
  test.seurat=get(load(file_name))

  return(test.seurat)
}


#' @title make.color.set
#' @description Get random colors
#' @usage make.color.sample(n, col.label)
#' @param n the number  of color sample, integer
#' @param col.label color label, character
#' @export
make.color.set <- function(n, col.label){
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  cl=sample(col_vector, n)
  names(cl)=col.label
  return(cl)
}

