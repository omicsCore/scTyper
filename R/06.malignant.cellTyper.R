#' @title malignant.cellTyper
#' @description A function to malignant cell typing
#' @param seurat Seurat object
#' @param rda.dir rData directory
#' @param malignant.cell.type Cell type to assign malignant cell
#' @param feature.to.test features to test as reference
#' @param cells.test_reference cells to test as reference
#' @details classification of malignant and non malignant seurat object.
#' @return Seurat object
#' @export
malignant.cellTyper <- function(seurat,
                                rda.dir = "./data",
                                malignant.cell.type="Epithelial",
                                feature.to.test = c("cell.type","tissue.type"),
                                cells.test_reference="immune"){


  message("[[",Sys.time(),"]] Run malignant.cellTyper --------")
  cell.type=as.character(seurat$cell.type)
  mal.fil.st = cell.type==malignant.cell.type

  if(feature.to.test=="tissue.type"){
    cnv.cut=quantile(seurat$cnv.score[seurat$cell.group %in% cells.test_reference], probs=0.90, na.rm = TRUE)
  }else if(feature.to.test=="cell.type"){
    cnv.cut=quantile(seurat$cnv.score[seurat$cell.group %in% c(cells.test_reference, "Unresolved_cell")], probs=0.90, na.rm = TRUE)
  }

  seurat$cnv.st=seurat$cnv.score>cnv.cut
  cnv.fil.st=seurat$cnv.st

  fil.st= (mal.fil.st | cnv.fil.st)
  table(fil.st)

  seurat$malignant.st=fil.st
  seurat$nonmalignant.st = !fil.st

  cell.type[fil.st]="Malignant_cell"
  seurat$cell.type=as.factor(cell.type)

  save(seurat, file = file.path(rda.dir, "seurat.rda"))
  message("[[",Sys.time(),"]] Finish malignant.cellTyper --------")
  return(seurat)
}
