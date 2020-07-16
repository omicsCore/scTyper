
# scTyper: a comprehensive pipeline for the cell typing analysis of single-cell RNA-seq data


# 1. Overview 
&nbsp;&nbsp;&nbsp;&nbsp; scTyper is a comprehensive pipeline for the cell typing and scRNA-Seq data analysis. It has been equipped with both database of cell type markers such as scTyper.db, CellMarker. Of note, markers for malignant cells, cancer-associated fibroblasts, and tumor-infiltrating T cells were collected in this database, which will be helpful in analyzing data from cancer tissues. In addition, scTyper provides three customized methods for estimating cell type marker expression, including the nearest template prediction (NTP), gene set enrichment analysis (GSEA), and average expression values. DNA copy number inference method (inferCNV), with improved modification, was also implemented in scTyper, which can be used for the typing of malignant cells. The package also supports the data preprocessing pipelines of Cell Ranger from 10X Genomics. Reporting system for analysis summary is also implemented which may facilitate users to perform reproducible analyses.  

# 2. Workflow

![](https://user-images.githubusercontent.com/36435306/84363831-3cec7000-ac0a-11ea-802d-41de1b953835.png)
</br>
&nbsp;&nbsp;&nbsp;&nbsp;scTyper is comprised of the modularized processes of ‘QC’, ‘Cell Ranger’, ‘Seurat processing’, ‘cell typing’, and ‘malignant cell typing’. These processes can be customized by manipulating the parameters for each process. If users want to perform only the cell typing process and a preprocessed input file with Seurat object is already prepared, the processing steps of ‘QC’, ‘Cell Ranger’ and ‘Seurat processing’ can be skipped by setting the parameters ‘qc’, ‘run.cellranger’ and ‘norm.seurat’ to ‘FALSE’. The processes and their parameters implemented in scTyper are summarized.

# 3. Getting Started 

## 3.1. External dependency

### 3.1.1 Program

`scTyper` uses `FastQC` and `Cell Ranger` to process scRNA-seq dataset, based on reference data.
If raw data processing is required, it should be installed (no need to install if not needed).

The Program and reference data can be found in the links below: </br>
**FastQC** : http://www.bioinformatics.babraham.ac.uk/projects/fastqc  </br>
**Cell Ranger** : https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest#loupe 


### 3.1.2 Required libraries

`scTyper` relies on the following dependencies which should be downloaded/updated.The current version of scTyper supports Linux operating system, due to the compatibility of the required software packages. Parallel computations on a multi-core machines can be used by calling the ‘parallel’ flag in R package.


| Processing                  | library                   | 
|----------------------------|--------------------|
|  **Quality Check** | fastqcr,  parallel |
|  **Cell Ranger** | parallel  |
|  **Seurat processing** | Seurat |
|  **Cell typing** | Seurat,  parallel |
|  **malignant.celltyping** | infercnv, perm, reshape2, gProfileR, GenomicRanges, grDevices, parallel  |
|  **Report** | rmarkdown, ComplexHeatmap, reshape2, pander, png, ggplot2, grid, gridExtra, circlize, knitr, kableExtra, colorspace |


## 3.2 Installation and loading
Source codes for scTyper are available at : https://github.com/omicsCore/scTyper

### 3.2.1 Installation package
scTyper runs in the R statistical computing environment. R version 3.5 or higher is required.

```r
# install BiocManager if necessary
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# install devtools if necessary
BiocManager::install("devtools")
library(devtools)

# install the scTyper package
devtools::install_github("omicsCore/scTyper")
```

### 3.2.2 Loading package and documentation


```r
# load
library("scTyper") 
library(help="scTyper")
```
 
# 4. More information
 
- [Sample analysis](http://htmlpreview.github.io/?https://github.com/omicsCore/scTyper/blob/master/vignettes/Sample_analysis.html)
- [Documentation](https://github.com/omicsCore/scTyper/files/4928754/Reference_manual.pdf)


