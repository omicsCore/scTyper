#' @title update.sig.db
#' @description Update sig.db
#' @usage update.sigTyper.db(sig.db.path, db.name=c("CellMarker", "sigTyper.db"), output.dir=system.file("/data",package = "scTyper") )
#' @param sig.db.path Path of sig.db.txt
#' @param output.dir storage path of sig.db marker
#' @export
update.sig.db <- function(sig.db.path, db.name=c("CellMarker", "sigTyper.db"), output.dir=system.file("/data",package = "scTyper") ){
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
