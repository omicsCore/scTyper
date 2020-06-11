
# 1. Installation
`scTyper` is implemented as an R package, scTyper, which can be installed from GitHub by:
R version 3.5 or higher is required.

```r
# install BiocManager if necessary
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# install devtools if necessary
BiocManager::install("devtools")
library(devtools)

# install the scTyper package
devtools::install_github("omicsCore/scTyper")

# load
library("scTyper") 
```

# 2. External dependency

## 2.1 Required Program

`scTyper` uses `FastQC` and `Cell Ranger` to process scRNA-seq dataset, based on reference data.

The Program and reference data can be found in the links below: </br>
**FastQC** : http://www.bioinformatics.babraham.ac.uk/projects/fastqc  </br>
**Cell Ranger** : https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest#loupe 


## 2.2 Required libraries

`scTyper` relies on the following dependencies which should be downloaded/updated.The current version of scTyper supports Linux operating system, due to the compatibility of the required software packages. Parallel computations on a multi-core machines can be used by calling the ‘parallel’ flag in R package.


| Processing                  | library                   | 
|----------------------------|--------------------|
|  **Global** | S4Vectors, e1071, IRanges, parallel |
|  **Quality Check** | fastqcr, grid |
|  **Cell Ranger** | parallel |
|  **Normalization** | Seurat |
|  **NTP** | Seurat |
|  **ES** | Seurat |
|  **Average** | Seurat |
|  **malignant.celltyping** | infercnv, utils, perm, gProfileR, GenomicRanges  |
|  **Make Set** | limma, Biobase |
|  **Report** | seurat, rmarkdown, ComplexHeatmap, reshape2, pander, ggplot2, grid, gridExtra, circlize, magrittr, knitr, kableExtra, RColorBrewer, png, grDevices, colorspace |
