
#' @title cell_type_NTP
#' @description A wrapper function to cell_type_NTP
#' @usage cell_type_NTP(seurat, wd, markerList, assay="RNA", slot=c("scale.data", "count.data", "data"), output.dir = "./", rda.dir = "./data", NTP.g.filter.method=c("sd", "mad", "none"), NTP.gene.filter.cutoff=0.3, NTP.distance=c("cosine","correlation"), NTP.norm.method=c("none", "row.std"), mc.cores=1)
#' @param seurat Seurat object
#' @param wd working directory
#' @param markerList List of cell type marker
#' @param assay Assay to use
#' @param slot seurat object expression data, c("scale.data", "count.data", "data")()
#' @param output.dir output directory
#' @param rda.dir Path of the RData saving directory
#' @param NTP.g.filter.method Method of gene filtering in NTP c(sd (Default), mad, none)
#' @param NTP.gene.filter.cutoff Cut-off score of standard deviation in NTP
#' @param NTP.distance Method of calculating distance in NTP, either c("correlation" or "cosine").
#' @param NTP.norm.method Method of normalization in NTP, either c("none", "row.std")
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details cell type annotation using NTP
#' @return  Seurat object
#' @export
cell_type_NTP <- function(seurat,
                          wd,
                          markerList,
                          assay="RNA",
                          slot=c("scale.data", "count.data", "data"),
                          output.dir = "./",
                          rda.dir = "./data",
                          NTP.g.filter.method=c("sd", "mad", "none"),
                          NTP.gene.filter.cutoff=0.3,
                          NTP.distance=c("cosine","correlation"),
                          NTP.norm.method=c("none", "row.std"),
                          mc.cores=1){
  message("[[",Sys.time(),"]] Run NTP --------")
  dir.create(output.dir, showWarnings = F)
  expr=GetAssayData(seurat, slot = slot)

  if(NTP.g.filter.method=="sd"){
    seurat[[assay]]@meta.features$sd = apply(expr, 1, sd)
    fil.st = seurat[[assay]]@meta.features$sd > NTP.gene.filter.cutoff
  }else if(NTP.g.filter.method=="mad"){
    seurat[[assay]]@meta.features$mad = apply(expr, 1, mad)
    fil.st = seurat[[assay]]@meta.features$mad > NTP.gene.filter.cutoff
  }else if(NTP.g.filter.method=="none"){
    fil.st = TRUE
  }

  ntp.dist.dir=file.path(output.dir, assay, NTP.distance)
  dir.create(ntp.dist.dir, showWarnings = F, recursive = T)

  res.ntp=NTP(eset = expr[fil.st,], sigList=markerList, out.dir = ntp.dist.dir, norm.method = NTP.norm.method, dist.selection = NTP.distance, mc.cores=mc.cores)
  if(NTP.distance=="cosine"){
    seurat@meta.data$ntp.cos=names(markerList)[res.ntp$pred.summary$predict.label]
    seurat@meta.data$ntp.fdr.cos=seurat$ntp.cos
    seurat@meta.data$ntp.fdr.cos[res.ntp$pred.summary$BH.FDR>0.05]="Unresolved_cell"
    seurat@meta.data$ntp.fdr.cos=factor(as.character(seurat$ntp.fdr.cos), levels=c(names(markerList), "Unresolved_cell"))
    seurat$ntp=seurat$ntp.cos
    seurat$ntp.fdr=seurat$ntp.fdr.cos
  }else if(NTP.distance=="correlation"){
    seurat@meta.data$ntp.cor=names(markerList)[res.ntp$pred.summary$predict.label]
    seurat@meta.data$ntp.fdr.cor=seurat$ntp.cor
    seurat@meta.data$ntp.fdr.cor[res.ntp$pred.summary$BH.FDR>0.05]="Unresolved_cell"
    seurat@meta.data$ntp.fdr.cor=factor(as.character(seurat$ntp.fdr.cor), levels=c(names(markerList), "Unresolved_cell"))
    seurat$ntp=seurat$ntp.cor
    seurat$ntp.fdr=seurat$ntp.fdr.cor
  }
  return(seurat)
}




#' @title NTP
#' @description A wrapper function to NTP
#' @usage NTP(eset, sigList, out.dir, output.name, dist.selection, norm.method, nresmpl, rnd.seed, mc.cores)
#' @param eset expressio
#' @param sigList gene signiture list
#' @param out.dir output directory
#' @param dist.selection calculating distance, a character, either c("correlation" or "cosine").
#' @param norm.method normalization method, either c("none", "row.std")
#' @param nresmpl an integer, number of permutations for \eqn{p}-value
#' @param rnd.seed Seed of the random number generator.
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details Nearest Template Prediction (NTP) based on predefined class templates.
#' @references Hoshida, Y. (2010). Nearest Template Prediction: A Single-Sample-Based Flexible Class Prediction with Confidence Assessment. PLoS ONE 5, e15543.
#' @return Nearest Template Prediction (NTP) result
#' @import parallel
#' @export
NTP<-function(
  eset,
  sigList, # sig n>=2
  out.dir,
  # temp.nn distance & row normalize
  dist.selection="cosine", # "correlation" or "cosine"
  #temp.nn.wt="T",   # only for 2 cls
  norm.method=c("none", "row.std"), # row.std

  # resampling to generate null dist
  nresmpl=1000,

  # Seed
  rnd.seed=1234,

  # mc.cores
  mc.cores=min(detectCores()-1, ncol(eset))
)
{

  # Advanced setting
  within.sig="F"

  # set random seed
  set.seed(rnd.seed)


  ### input ###

  # selected features used for prediction
  features=melt(sigList)

  features=data.frame(ProbeID=as.character(features$value), GeneName=as.character(features$value), Class=as.numeric(factor(features$L1, levels = names(sigList))))

  num.features<-nrow(features)
  num.cls<-length(sigList)

  print(paste("num.features:", num.features))
  print(paste("num.cls:"))
  print(table(features$Class))

  features<-data.frame(ord=seq(1:num.features), features)
  head(features)


  # expression data

  ## Other dataset's mean & SD for row normalization (optional)

  ProbeID<-rownames(eset)
  num.samples<-ncol(eset)
  exp.dataset<-eset
  sample.names<-colnames(eset)


  # row normalize
  normed.exp.dataset<-exp.dataset

  if (norm.method=="row.std"){
    exp.mean <- apply(exp.dataset,1,mean,na.rm=T)
    exp.sd <- apply(exp.dataset,1,sd,na.rm=T)
    normed.exp.dataset<-(exp.dataset-exp.mean)/exp.sd
  } else if(norm.method=="none"){}

  normed.exp.dataset <- data.frame(ProbeID, normed.exp.dataset)

  # extract features from normed.exp.dataset

  exp.dataset.extract<-merge(features,normed.exp.dataset,sort=F)

  if (length(exp.dataset.extract[,1])<1){
    stop("### No matched probes! ###")
  }

  exp.dataset.extract<-exp.dataset.extract[order(exp.dataset.extract$ord),]
  order.extract.after<-exp.dataset.extract$ord
  exp.dataset.extract<-exp.dataset.extract[-2]


  features.extract<-exp.dataset.extract[,1:3]
  features.extract<-cbind(order.extract.after, features.extract) # order:ProbeID:gene name:cls:wt(if any)
  num.features.extract<-nrow(features.extract)

  ProbeID.extract<-as.vector(exp.dataset.extract$ProbeID)
  exp.dataset.extract<-exp.dataset.extract[,-c(1:3)]

  if(sum(duplicated(ProbeID.extract))>0){
    stop("### Marker genes overlap. Use the ES or MEAN method! ###")
  }else{rownames(exp.dataset.extract)<-ProbeID.extract}



  # make template

  for (i in 1:num.cls){
    temp.temp<-as.vector(features.extract$Class)
    temp.temp[temp.temp!=i]<-0
    temp.temp[temp.temp==i]<-1
    eval(parse(text=paste("temp.",i,"<-temp.temp",sep="")))
    #    eval(parse(text=paste("temp\.",i,"<-temp\.temp",sep="")))  ### for < R-2.4.0
  }
  table(temp.1)
  table(temp.2)

  ### compute distance and p-value ###

  resList=mclapply(1:num.samples, function(i){
    print(paste("sample # ", i, sep=""))

    current.sample <- as.numeric(as.vector(exp.dataset.extract[,i]))

    # compute original distance

    orig.dist.to.all.temp <- vector(length=num.cls,mode="numeric")

    if (dist.selection=="cosine"){
      for (o in 1:num.cls){      # compute distance to all templates
        eval(parse(text=paste("current.temp <- temp.",o,sep="")))
        orig.dist.to.all.temp[o]<-sum(current.temp*current.sample)/(sqrt(sum(current.temp^2))*sqrt(sum(current.sample^2)))
      }

    }
    if (dist.selection=="correlation"){
      for (o in 1:num.cls){      # compute distance to all templates
        eval(parse(text=paste("current.temp <- temp.",o,sep="")))
        orig.dist.to.all.temp[o] <- cor(current.temp,current.sample,method="pearson",use="complete.obs")
      }
    }

    if (num.cls==2){           # find nearest neighbor (2 classes)
      if (orig.dist.to.all.temp[1]>=orig.dist.to.all.temp[2]){
        predict.label<-1
        dist.to.template<-1-orig.dist.to.all.temp[1]
        dist.to.cls1<--(orig.dist.to.all.temp[1]+1)
      }
      if (orig.dist.to.all.temp[1]<orig.dist.to.all.temp[2]){
        predict.label<-2
        dist.to.template<-1-orig.dist.to.all.temp[2]
        dist.to.cls1<-orig.dist.to.all.temp[2]+1
      }
    }

    if (num.cls>2){
      for (o in 1:num.cls){       # find nearest neighbor (>2 classes)
        if (is.na(orig.dist.to.all.temp[o])!=T){
          if (orig.dist.to.all.temp[o]==max(orig.dist.to.all.temp,na.rm=T)){
            predict.label<-o
            dist.to.template<-1-orig.dist.to.all.temp[o]
            dist.to.cls1<-(1-orig.dist.to.all.temp[o])+o
          }
        }
      }
    }

    # permutation test
    rnd.feature.matrix<-matrix(0,nrow=num.features.extract,ncol=nresmpl)
    perm.dist.vector<-vector(length=nresmpl*num.cls,mode="numeric")

    if (within.sig=="F"){     # generate resampled features from all probes
      for (p in 1:nresmpl){
        rnd.feature.matrix[,p]<-sample(as.numeric(normed.exp.dataset[,(i+1)]),num.features.extract,replace=F)
      }
    }
    if (within.sig=="T"){     # generate resampled features from only signature genes
      for (p in 1:nresmpl){
        rnd.feature.matrix[,p]<-sample(as.numeric(exp.dataset.extract[,i]),num.features.extract,replace=F)
      }
    }

    # compute distance to all templates
    if (dist.selection=="cosine"){          # cosine
      for (res in 1:num.cls){
        eval(parse(text=paste("temp.resmpl<-temp.",res,sep="")))
        prod.sum<-apply(t(t(rnd.feature.matrix)*temp.resmpl),2,sum)
        data.sq.sum<-apply(rnd.feature.matrix^2,2,sum)
        temp.sq.sum<-sum(temp.resmpl^2)

        perm.dist.vector[(1+(nresmpl*(res-1))):(nresmpl*res)]<-(1-(prod.sum/(sqrt(data.sq.sum)*sqrt(temp.sq.sum))))
      }
    }

    if (dist.selection=="correlation"){          # correlation
      for (res in 1:num.cls){
        eval(parse(text=paste("temp.resmpl<-temp.",res,sep="")))
        perm.dist.vector[(1+(nresmpl*(res-1))):(nresmpl*res)]<-(1-as.vector(cor(rnd.feature.matrix,temp.resmpl,method="pearson",use="complete.obs")))
      }
    }

    # compute nominal p-value

    combined.stats.rank<-rank(c(dist.to.template,perm.dist.vector))
    nominal.p<-combined.stats.rank[1]/length(combined.stats.rank)

    c(predict.label=predict.label, dist.to.template=dist.to.template, dist.to.cls1=dist.to.cls1, nominal.p=nominal.p)
  }, mc.cores=mc.cores) # main sample loop END

  predict.label=sapply(resList, function(a) a["predict.label"])
  dist.to.template=sapply(resList, function(a) a["dist.to.template"])
  dist.to.cls1=sapply(resList, function(a) a["dist.to.cls1"])
  nominal.p=sapply(resList, function(a) a["nominal.p"])

  # MCT correction

  BH.FDR<-nominal.p*num.samples/rank(nominal.p)
  Bonferroni.p<-nominal.p*num.samples

  BH.FDR[BH.FDR>1]<-1
  Bonferroni.p[Bonferroni.p>1]<-1

  ### output ###

  # prediction results

  dist.to.cls1.rank <- rank(dist.to.cls1)
  pred.summary <- cbind(sample.names,predict.label,dist.to.template,dist.to.cls1.rank,
                        nominal.p,BH.FDR,Bonferroni.p)

  write.table(pred.summary, file.path(out.dir,"NTP_prediction_result.txt"), quote=FALSE,sep="\t",row.names=FALSE)

  # extracted features
  write.table(features.extract[,2:4], file.path(out.dir,"NTP_features.txt"), quote=FALSE,sep="\t",row.names=FALSE)
  features=features.extract[,2:4]

  pred.summary=data.frame(sample.names=pred.summary[,1], apply(pred.summary[,-1], 2, as.numeric))
  list(pred.summary=pred.summary, features=features)

}  # END main
