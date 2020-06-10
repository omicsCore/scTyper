
# scTyper: an analysis pipeline for cell typing of single-cell RNA-seq data
____________________________________________________________________________________________________________________________</br>


# 1. Overview
&nbsp;&nbsp;&nbsp;&nbsp;Advancements in single-cell RNA sequencing (scRNA-Seq) technology has enabled us identifying individual cell types from a complex cell population in tissues. Cell typing is one of the key challenges in scRNA-Seq data analysis, which is usually performed by estimating cell marker gene expression. However, there is no standard practice for cell typing analysis and using different cell makers and cell typing algorithms results in variable and inaccurate results. scTyper is a comprehensive pipeline for the cell typing and scRNA-Seq data analysis. It has been equipped with a database of cell type markers i.e., scTyper.db, containing 213 cell markers, collected from previous studies. Of note, markers for malignant cells, cancer-associated fibroblasts, and tumor-infiltrating T cells were collected in this database, which will be helpful in analyzing data from cancer tissues. In addition, scTyper provides three customized methods for estimating cell type marker expression, including the nearest template prediction (NTP), gene set enrichment analysis (GSEA), and average expression values. DNA copy number inference method (inferCNV), with improved modification, was also implemented in scTyper, which can be used for the typing of malignant cells. The package also supports the data preprocessing pipelines of Cell Ranger from 10X Genomics. Reporting system for analysis summary is also implemented which may facilitate users to perform reproducible analyses.  





# 2. Implementation
&nbsp;&nbsp;&nbsp;&nbsp;scTyper is an R package that supports cell typing analysis from scRNA-seq data. The package is designed to be run by a single command, otherwise the pipeline can be customized at each step by choosing user-defined options. scTyper also includes the pipeline for quality control and sequence alignment, supported by FASTQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc) and Cell Ranger [14], respectively. Data normalization, clustering, and visualization steps are also included in the pipelines, which are linked with- and performed by Seurat R package [15].




# 3. Workflow


![](https://user-images.githubusercontent.com/36435306/83315120-db66f180-a258-11ea-83d6-716428dd9e87.png)
</br>
&nbsp;&nbsp;&nbsp;&nbsp;First, we constructed a cell marker database scTyper.db, which included the manually curated 213 gene sets and the 122 cell types from 22 studies. For user convenience, the package also implemented the data from CellMarker database including 2,867 gene sets for 467 cell types [10]. Users can readily select the cell markers from the DB and apply them to the cell typing pipeline in the package. For cell typing analysis, cell marker gene expression can be estimated by three different methods of NTP [16], pre-ranked GSEA [17], and average gene expression. For malignant cell typing, inferred DNA copy numbers can be used in the pipeline, which were estimated by an inferCNV R package [18] with modification (for details see Results). 
 The package contained analyses pipeline steps for quality control, alignment and quantification of raw sequencing data, and cell typing. These processes can be executed by one-step command. Data processing steps for log transformation, normalization, and clustering were performed by the functions in ‘Seurat’ R package, which were wrapped in the package. Analysis results are automatically documented as a report with the summaries of the processing steps, visualization tables and clusters visualized by t-SNE plots. 




# 4. Getting Started 

## 4.1 Required libraries
scTyper is an R package and requires various Bioconductor libraries. Before running scTyper, the user must install the appropriate reference data including `FastQC`, `Cell Ranger` and several libraries. The current version of scTyper supports Linux operating system, due to the compatibility of the required software packages. Parallel computations on a multi-core machines can be used by calling the ‘parallel’ flag in R package.

The reference data packages required for scTyper can be found in the links below: </br>
**FastQC** : http://www.bioinformatics.babraham.ac.uk/projects/fastqc  </br>
**Cell Ranger** : https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest#loupe </br>

scTyper relies on the following dependencies which should be downloaded/updated before installing scTyper.

| Processing                  | library                   | 
|----------------------------|--------------------|
|  **Global** | S4Vectors, e1071, IRanges, parallel |
|  **Quality Check** | fastqcr, grid |
|  **Cell Ranger** | parallel |
|  **Normalization** | Seurat |
|  **NTP** | Seurat |
|  **ES** | Seurat, perm |
|  **Average** | Seurat |
|  **malignant.celltyping** | infercnv, utils, perm, gProfileR, GenomicRanges  |
|  **Make Set** | limma, Biobase |
|  **Report** | seurat, rmarkdown, ComplexHeatmap, reshape2, pander, ggplot2, grid, gridExtra, circlize, magrittr, knitr, kableExtra, RColorBrewer, png, grDevices, colorspace |


## 4.2 Installation and loading
Source codes for scTyper are available at : https://github.com/omicsCore/scTyper

### 4.2.1 Required software
scTyper runs in the R statistical computing environment. R version 3.5 or higher is required.

```r
# Enter commands in R (or R studio, if installed)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("devtools")
```

### 4.2.2 Installing the latest version


```r
library(devtools)
devtools::install_github("omicsCore/scTyper")
```

### 4.2.3 Loading package and documentation


```r
library("scTyper") 
library(help="scTyper")
```


## 4.3 Prepare input data
### 4.3.1 Test dataset

The raw data for the test data set used here can be downloaded from  “https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322” and further organized to run it as a test. The pre-processed (QC and other initial steps have been already processed), consists of 5,902 cells (4,541 cancer cells and 1,361 lymph node cells) and 23,686 genes extracted from 17 patients is available at GitHub **(scTyper/data/GSE103322.seurat.rda)** .



### 4.3.2 Phenotype data
Users need to make the sample phenotype information table as input for pipeline in a tabular form  **‘.csv’** format. The order of the columns must follow the example template and the names of the columns be formatted as <span style=" font-weight: bold;    color: #130CA3 !important;" >Sample_ID</span> and <span style=" font-weight: bold;    color: #130CA3 !important;" >TissueType</span>.
We use [Pheotype information (click link)](https://github.com/omicsCore/scTyper/raw/master/data/pheno_info_public.csv) to testdata. </br>
</br>
&nbsp;&nbsp;**Sample_ID** : Sample name of single cell RNA sequencing dataset.  </br>
&nbsp;&nbsp;**TissueType** : Type of tissue that is composed of cells with similar structure and act together to perform a specific function. </br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ex) Normal, Cancer etc. </br>

 - **Example** </br>
 
| Sample_ID                  | TissueType       | 
|-------------------|-------------|
|  sample_1 | Cancer |
|  sample_2 | Cancer |
|  sample_3 | Normal |
|  sample_4 | Normal |
|  sample_5 | Normal |

### 4.3.4 Cell type marker set
&nbsp;&nbsp;&nbsp;&nbsp;After the data processing, the cell typing process could be performed using cell markers, collected in **scTyper.db** and **CellMarker** databases. User can easily search the markers of interest from these databases, which are uniquely labeled using unified nomenclature. The input of marker parameter is elements of `#Identifier`, `#StudyName` columns or user can define marker list as input parameter. While selecting the combined or a single markers, the duplicated genes, if there is any, are excluded in cell typing process. If there is one gene in one cell type, the “Average” cell typing method should be used.

Also, user can customize markers by grouping variable names with c(). </br>

 - **Example** </br>

```r
    marker="Puram.2017.HNSCC.TME"
    marker=c("Kawai.2018.Liver", "Li.2017.CRC.TME") #StudyName
    marker=c("Costea.2013.OSCC.CAF:Normal cell:NF", "Costea.2013.OSCC.CAF:Cancer cell:CAF_D", #Identifier
              "Elyada.2019.PDAC.CAF:Cancer cell:iCAF", "Elyada.2019.PDAC.CAF:Cancer cell:myCAF")
    marker=list("T_cell"=c("CD3D", "TRAC", "TRBC1", "TRBC2"),  #user defined marker list 
                "B_cell"=c("CD79A", "IGKC", "IGLC3", "IGHG3"), 
                "Endothelial_cell"=c("CLDN5", "FLT1", "CDH5", "RAMP2"),
                "Epithelial_cell"=c("KRT14", "KRT17", "KRT6A", "KRT5", "KRT19", "KRT8", "KRT16",
                                  "KRT18", "KRT6B", "KRT15", "KRT6C", "KRTCAP3", "EPCAM", "SFN")
```

[<span style=" font-weight: bold;    color: #0060CC !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #C1CDCD !important;" >'cell marker database.xlsx (click link)'</span>](https://github.com/omicsCore/scTyper/files/4723417/cell.marker.database.xlsx)
The link is marker databases comprised of "scTyper.db" and "CellMarker DB".
The marker databases are comprised of total 3,080 cell marker gene set collected from 1,764 CellMarker studies and 22 scTyper.db studies.

### 4.4.5 scTyper() arguments 
To execute **scTyper**, some basic input arguments are required for proper execution. 
User should get some knowledge and set the corresponding argument for the required process.

<table class="table table-striped" style="width: auto !important; ">
 <thead>
  <tr>
   <th style="text-align:left;"> Process </th>
   <th style="text-align:left;"> Parameters </th>
   <th style="text-align:left;"> Description </th>
   <th style="text-align:left;"> Values </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> global configuration </td>
   <td style="text-align:left;"> wd </td>
   <td style="text-align:left;"> Working directory </td>
   <td style="text-align:left;"> Character </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> output.name </td>
   <td style="text-align:left;"> Output directory name </td>
   <td style="text-align:left;"> Character </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> pheno.fn </td>
   <td style="text-align:left;"> Phenotype file path </td>
   <td style="text-align:left;"> File path </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> qc, run.cellranger , norm.seurat </td>
   <td style="text-align:left;"> Indicate whether the process run </td>
   <td style="text-align:left;"> Logical (Default = ‘FALSE’) </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> cell.typing.method </td>
   <td style="text-align:left;"> Cell typing method </td>
   <td style="text-align:left;"> ‘NTP’ (Default), ‘ES’, ‘Average’ </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> level </td>
   <td style="text-align:left;"> Indicate the cell assignment level (cell or cluster) </td>
   <td style="text-align:left;"> ‘cell’ (Default), ‘cluster’ </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> run.inferCNV </td>
   <td style="text-align:left;"> Indicate whether ‘malignant cell typing by inferCNV process run </td>
   <td style="text-align:left;"> Logical (Default = ‘TRUE’) </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> project.name </td>
   <td style="text-align:left;"> Project name </td>
   <td style="text-align:left;"> Character </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> mc.cores </td>
   <td style="text-align:left;"> Number of cores </td>
   <td style="text-align:left;"> Numeric (Default = ‘1’) </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> report.mode </td>
   <td style="text-align:left;"> Generate report file </td>
   <td style="text-align:left;"> Logical (Default = ‘TRUE’) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> QC </td>
   <td style="text-align:left;"> fastqc.path </td>
   <td style="text-align:left;"> FastQC program path </td>
   <td style="text-align:left;"> File path </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> fastq.dir </td>
   <td style="text-align:left;"> FastQC output directory </td>
   <td style="text-align:left;"> File path </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> fq1.idx </td>
   <td style="text-align:left;"> Index of the FASTQ file (Read 1) </td>
   <td style="text-align:left;"> Character (Default = ‘_R1_001.fastq’) </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> fq2.idx </td>
   <td style="text-align:left;"> Index of the FASTQ file (Read 2) </td>
   <td style="text-align:left;"> Character (Default = ‘_R2_001.fastq’) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CellRanger </td>
   <td style="text-align:left;"> cellranger.path </td>
   <td style="text-align:left;"> Cell Ranger program path </td>
   <td style="text-align:left;"> File path </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> cellranger.ref.dir </td>
   <td style="text-align:left;"> Directory of Cell Ranger reference file </td>
   <td style="text-align:left;"> File path </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Seurat processing </td>
   <td style="text-align:left;"> percent.min.cells </td>
   <td style="text-align:left;"> Cutoff to filter features containing minimum percent of cells </td>
   <td style="text-align:left;"> 0.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> min.features </td>
   <td style="text-align:left;"> Cutoff to filter cells containing minimum number of features </td>
   <td style="text-align:left;"> 200 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> percent.mt </td>
   <td style="text-align:left;"> Cutoff for filtering cells that have &gt;n percent mitochondrial counts </td>
   <td style="text-align:left;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> vars.to.regress </td>
   <td style="text-align:left;"> Variables to regress out </td>
   <td style="text-align:left;"> Default=c(‘nCount_RNA’, ‘percent.mt’) </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> dims </td>
   <td style="text-align:left;"> A vector of the dimensions to use in construction of the SNN graph. </td>
   <td style="text-align:left;"> 1:100 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> resolution </td>
   <td style="text-align:left;"> Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. </td>
   <td style="text-align:left;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cell typing </td>
   <td style="text-align:left;"> seurat.object </td>
   <td style="text-align:left;"> Seurat object </td>
   <td style="text-align:left;"> Seurat object </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> slot </td>
   <td style="text-align:left;"> Data type of Seurat object </td>
   <td style="text-align:left;"> ‘scale.data’ (Default), ‘count.data’, ‘data’ </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> marker </td>
   <td style="text-align:left;"> Cell markers to use cell typing </td>
   <td style="text-align:left;"> Character or List (Signature names or Study names or User defined gene set list) </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> assay </td>
   <td style="text-align:left;"> Assay of Seurat object </td>
   <td style="text-align:left;"> Character (Default=’RNA’) </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NTP.g.filter.method </td>
   <td style="text-align:left;"> Method to filter genes in NTP </td>
   <td style="text-align:left;"> ‘sd’ (Default),’mad’, ‘none’ </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NTP.gene.filter.cutoff </td>
   <td style="text-align:left;"> Cutoff to filter genes of in NTP </td>
   <td style="text-align:left;"> Numeric (Default = ‘0.3’) </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NTP.distance </td>
   <td style="text-align:left;"> NTP distance method </td>
   <td style="text-align:left;"> ‘cosine’ (Default), ‘correlation’ </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NTP.norm.method </td>
   <td style="text-align:left;"> NTP normalization method </td>
   <td style="text-align:left;"> ‘none’ (Default), ‘row.std’ </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Malignant cell typing (inferCNV </td>
   <td style="text-align:left;"> gene.ref.gtf </td>
   <td style="text-align:left;"> Path of GTF file including genomic location for genes </td>
   <td style="text-align:left;"> File path </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> feature.to.test </td>
   <td style="text-align:left;"> Column header name of the meta data in Seurat object (select the cell groups for T.test) </td>
   <td style="text-align:left;"> Character (Default = ‘cell.type’), ‘tissue.type’ </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> cells.test_excluded </td>
   <td style="text-align:left;"> A value indicates the cells to be excluded in T.test </td>
   <td style="text-align:left;"> Character (Default = ‘Epithelial’) </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> cells.test_reference </td>
   <td style="text-align:left;"> A value indicates the cells to use as be excluded in T.test </td>
   <td style="text-align:left;"> Character (Default = ‘immune’) </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> fc.cutoff </td>
   <td style="text-align:left;"> Cutoff of fold change </td>
   <td style="text-align:left;"> Numeric (Default = ‘0.05’) </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> cutoff.gene.cluster </td>
   <td style="text-align:left;"> A cutoff P-value for filtering out the gene clusters (calculated from GO analysis) </td>
   <td style="text-align:left;"> Numeric (Default = ‘0.05’) </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> malignant.cell.type </td>
   <td style="text-align:left;"> Cell type to assign malignant cell </td>
   <td style="text-align:left;"> Character (Default = ’Epithelial’) </td>
  </tr>
</tbody>
</table>
Note: When users have pre-processed seurat object (completed with qc and cellranger), **seurat.object** parameter must be assigned.
</br>
</br>
</br>
&nbsp;&nbsp;&nbsp;&nbsp;After setting, users need to set parameters such as fastqc.path, cellranger.path, cellranger.ref.dir, gene.ref.gtf. The executable file paths of softwares should be setted on this script file. If you don't know the executable file path, type 'which fastqc' on the shell or system (e.g.'which fastqc') on your R session. 

</br>
 - **Example** </br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;fastq.dir="/data/Rpackage/scTyper/test" </br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;fastqc.program.path= "/data/program/bin/fastqc" </br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cellranger.program.path = "/data/program/bin/cellranger" </br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cellranger.ref.dir = "/data/pubdata/ngs_ref/cellranger/refdata-cellranger-GRCh38-1.2.0" </br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;gene.ref.gtf="/data/pubdata/ngs_ref/cellranger/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf" </br>


## 4.5 Workflow process at a glance
Below, you can see the code that is designed to run by a single command, `scTyper()`.  </br>
This function is based on the arguments mentioned above and is run in Rscript.  </br>

### 4.5.1 Load test data


```r
library(scTyper)
test.seurat = get(load("path.test.data/GSE103322.seurat.rda"))
```


### 4.5.2 Run scTyper() 
scTyper can be run as a single command  like a following example.
 - ***Runing (pre-processing skip, NTP cell typing with inferCNV at cell level)***

```r
celltyped.seurat=scTyper(wd = "/data/Rpackage/scTyper",
                         output.name = "cell.typed.seurat.result",
                         pheno.fn = "/data/Rpackage/scTyper/data/pheno_info_public.csv",
                         qc = FALSE,
                         run.cellranger=FALSE,
                         norm.seurat=FALSE,
                         cell.typing.method="NTP",
                         level="cell",
                         run.inferCNV=TRUE,
                         proj.name = "scTyper",
                         seurat.object=test.seurat,
                         marker="Puram.2017.HNSCC.TME",
                         gene.ref.gtf="/data/pubdata/ngs_ref/cellranger/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf",
                         feature.to.test = "cell.type",
                         cells.test_excluded=c("Epithelial_cell", "Unresolved_cell"),
                         cells.test_reference = c("Fibroblasts", "T.cells", "B_Plasma.cells", "Macrophages", "Dendritic.cells", "Endothelial.cells", "Myocytes", "Mast"),
                         malignant.cell.type="Epithelial_cell",
                         report.mode=TRUE,
                         mc.cores = 10)
```

### 4.5.3 scTyper() output 

 - ***final output of scTyper() is a Seurat object*** : scTyper is implemented as a wrapper function in R package, which may facilitate the subsequent analysis under R environment. Furthermore, scTyper provides a function to transform the processed data into a **Seurat object**, which is a popular data type for the subsequent analyses and biological interpretation of scRNA-seq. The Seurat object is a class allowing for the storage and manipulation of single-cell data. Seurat object was designed to allow for greater flexibility to work with all these data types in a cohesive framework. 

 - ***Report file*** : Finally report.mode is a function used to combine the results into one unified file. scTyper generates a report summary automatically. The document summarizes each step of the processes, the parameters used, and the results of cell typing and clustering, and visualizes the results by heatmaps and t-SNE plots (Supplementary Data). This may help users reproduce their analysis workflows.
 
 - ***Automatically generated output directory*** : The result files are created automatically in output directory that users setted by parameter ‘wd’ and ‘output.name’ according to your processing step.</br>
 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- 00_qc : Directory of FastQC output </br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- 01_count : Directory of CellRanger output </br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- 02_NTP : Directory of NTP output </br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- 03_inferCNV : Directory of inferCNV output </br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- data : Directory of RData (rda). 'seurat.rda' (final output) saved. </br>


### 4.5.4 Result visualization
If ***“report.mode=TRUE”*** is set true in the scTyper(), the user will get a seurat object as an output.

- **Final seurat object to use for analysis**


```r
celltyped.seurat
```

```
An object of class Seurat 
23686 features across 5902 samples within 1 assay 
Active assay: RNA (23686 features, 2000 variable features)
 2 dimensional reductions calculated: pca, tsne
```

- **Cell type statistics across samples**


```r
t.tbl=table(celltyped.seurat$cell.type)
tbl=cbind(t.tbl, table(celltyped.seurat$cell.type, celltyped.seurat$sample.name))
colnames(tbl)[1]="Total"
panderOptions('table.split.table', Inf)
pander(tbl, split.cells = 10)
```

- ***Make table of inferred cell types***


```r
tbl=table(celltyped.seurat$cell.type)
s.tbl=table(celltyped.seurat$cell.type, celltyped.seurat$sample.name); r.tbl=apply(s.tbl, 2, function(a) a/sum(a))
```

- **Distribution of inferred cell types**


```r
par(mar=c(7,5,0,0));barplot(tbl, las=2, ylab="Number of cells", cex.names=0.8, cex.axis=0.8)
```
![](https://user-images.githubusercontent.com/36435306/83992767-58e2de00-a98c-11ea-9d32-3b63bb8dfad6.png)


- **Distibution of cell types across samples**


```r
par(mar=c(7,5,0,10))
barplot(r.tbl, las=2, ylab="Proportion (%)", cex.names=0.8, cex.axis=0.8, legend.text = rownames(r.tbl), args.legend = list(x = 'right', bty='n', inset=c(-0.3,0), xpd = TRUE, cex=0.7))
```
![](https://user-images.githubusercontent.com/36435306/83992781-64360980-a98c-11ea-9e1e-5d1bb6b6698e.png)


- **inferred cell types by cell typing method**


```r
cols=rainbow_hcl(length(levels(celltyped.seurat$cell.type))); names(cols)=levels(celltyped.seurat$cell.type)
p1=DimPlot(celltyped.seurat, reduction = 'tsne',group.by="cell.type",label=F,cols=cols, pt.size=0.2); p1=LabelClusters(plot = p1, id = "cell.type", size = 3)
p1
```
![](https://user-images.githubusercontent.com/36435306/83992830-99daf280-a98c-11ea-8081-8f334671a8e3.png)

- **Malignant cells by inferCNV**


```r
p4=DimPlot(celltyped.seurat, reduction = 'tsne',group.by="malignant.st",label=F, pt.size=0.2, cols=c("grey90", "red")) 
p4
```
![](https://user-images.githubusercontent.com/36435306/83992848-a2cbc400-a98c-11ea-9365-9ae5ac515e1d.png)

- **CNV score by inferCNV**


![](https://user-images.githubusercontent.com/36435306/83992887-b8d98480-a98c-11ea-99f7-9f94943ba558.png)

- **Seurat clusters**


```r
p2=DimPlot(celltyped.seurat, reduction = 'tsne',group.by="seurat_clusters",label=F, pt.size=0.2); p2=LabelClusters(plot = p2, id = "seurat_clusters", size = 3)
p2
```
![](https://user-images.githubusercontent.com/36435306/83992903-ce4eae80-a98c-11ea-948e-bffe52ceae39.png)
- **Samples**


```r
p3=DimPlot(celltyped.seurat, reduction = 'tsne',group.by="sample.name",label=F, pt.size=0.2)
p3
```
![](https://user-images.githubusercontent.com/36435306/83992908-d27acc00-a98c-11ea-96db-932ea1836619.png)

- **Cell markers heatmap**
A heatmap shows the cell typing result and the gene expression levels of cell marker gene sets from Puram.2017.HNSCC.TME. For each method, the assigned cell types are indicated by color bars.  


```r
draw.heatmap(seurat = celltyped.seurat,
             wd = "/data/Rpackage/scTyper",
             run.inferCNV = TRUE,
             slot = "scale.data",
             marker="Puram.2017.HNSCC.TME")
```
![](https://user-images.githubusercontent.com/36435306/83992910-d4448f80-a98c-11ea-8aba-0f403cd2b994.png)




## 4.6  Example script for public dataset with options.
 - ***Runing (pre-processing skip, NTP cell typing with inferCNV at cell level)***

```r
celltyped.seurat=scTyper(wd = "/data/Rpackage/scTyper",
                         output.name = "cell.typed.seurat.result",
                         pheno.fn = "/data/Rpackage/scTyper/data/pheno_info_public.csv",
                         qc = FALSE,
                         run.cellranger=FALSE,
                         norm.seurat=FALSE,
                         cell.typing.method="NTP",
                         level="cell",
                         run.inferCNV=TRUE,
                         proj.name = "scTyper",
                         seurat.object=test.seurat,
                         marker="Puram.2017.HNSCC.TME",
                         gene.ref.gtf="/data/pubdata/ngs_ref/cellranger/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf",
                         feature.to.test = "cell.type",
                         cells.test_excluded=c("Epithelial_cell", "Unresolved_cell"),
                         cells.test_reference = c("Fibroblasts", "T.cells", "B_Plasma.cells", "Macrophages", "Dendritic.cells", "Endothelial.cells", "Myocytes", "Mast"),
                         malignant.cell.type="Epithelial_cell",
                         report.mode=TRUE,
                         mc.cores = 10)
```

 - ***Runing (pre-processing skip, ES cell typing with inferCNV at cell level)***

```r
celltyped.seurat=scTyper(wd = "/data/Rpackage/scTyper",
                         output.name = "cell.typed.seurat.result",
                         pheno.fn = "/data/Rpackage/scTyper/data/pheno_info_public.csv",
                         qc = FALSE,
                         run.cellranger=FALSE,
                         norm.seurat=FALSE,
                         cell.typing.method="ES",
                         level="cell",
                         run.inferCNV=TRUE,
                         proj.name = "scTyper",
                         seurat.object=test.seurat,
                         marker="Puram.2017.HNSCC.TME",
                         gene.ref.gtf="/data/pubdata/ngs_ref/cellranger/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf",
                         feature.to.test = "cell.type",
                         cells.test_excluded=c("Epithelial_cell", "Unresolved_cell"),
                         cells.test_reference = c("Fibroblasts", "T.cells", "B_Plasma.cells", "Macrophages", "Dendritic.cells", "Endothelial.cells", "Myocytes", "Mast"),
                         malignant.cell.type="Epithelial_cell",
                         report.mode=TRUE,
                         mc.cores = 10)
```
 - ***Runing (pre-processing skip, Average cell typing with inferCNV at cell level)***

```r
celltyped.seurat=scTyper(wd = "/data/Rpackage/scTyper",
                         output.name = "cell.typed.seurat.result",
                         pheno.fn = "/data/Rpackage/scTyper/data/pheno_info_public.csv",
                         qc = FALSE,
                         run.cellranger=FALSE,
                         norm.seurat=FALSE,
                         cell.typing.method="Average",
                         level="cell",
                         run.inferCNV=TRUE,
                         proj.name = "scTyper",
                         seurat.object=test.seurat,
                         marker="Puram.2017.HNSCC.TME",
                         gene.ref.gtf="/data/pubdata/ngs_ref/cellranger/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf",
                         feature.to.test = "cell.type",
                         cells.test_excluded=c("Epithelial_cell", "Unresolved_cell"),
                         cells.test_reference = c("Fibroblasts", "T.cells", "B_Plasma.cells", "Macrophages", "Dendritic.cells", "Endothelial.cells", "Myocytes", "Mast"),
                         malignant.cell.type="Epithelial_cell",
                         report.mode=TRUE,
                         mc.cores = 10)
```
 - ***Runing (pre-processing(QC, CellRanger, Normalization) , NTP cell typing with inferCNV at cell level)***

```r
###### Running all step (User's data)
library(scTyper)
processed.celltyped.seurat=scTyper(wd="/data/Rpackage/scTyper",
                                   output.name = "test.result",
                                   pheno.fn="/data/Rpackage/scTyper/data/pheno_info_test.csv",
                                   qc = TRUE,
                                   run.cellranger=TRUE,
                                   norm.seurat=TRUE,
                                   cell.typing.method="NTP",
                                   level="cell",
                                   run.inferCNV=TRUE,
                                   proj.name = "scTyper",
                                   fastqc.path="/data/program/bin/fastqc",
                                   fastq.dir="/data/Rpackage/scTyper/test",
                                   fq1.idx="_R1_001.fastq",
                                   fq2.idx="_R2_001.fastq",
                                   cellranger.path="/data/program/bin/cellranger",
                                   cellranger.ref.dir="/data/pubdata/ngs_ref/cellranger/refdata-cellranger-GRCh38-1.2.0",
                                   marker="Puram.2017.HNSCC.TME",
                                   gene.ref.gtf="/data/pubdata/ngs_ref/cellranger/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf",
                                   feature.to.test = "tissue.type",
                                   cells.test_excluded=c("Epithelial_cell", "Unresolved_cell"),
                                   cells.test_reference = "Normal",
                                   malignant.cell.type="Epithelial_cell",
                                   report.mode=TRUE,
                                   mc.cores = 10)
```
Help information is available with `scTyper --help`.





# 9 References
1.	Hwang B, Lee JH, Bang D: Single-cell RNA sequencing technologies and bioinformatics pipelines. Experimental & molecular medicine 2018, 50(8):96.
2.	Abdelaal T, Michielsen L, Cats D, Hoogduin D, Mei H, Reinders MJT, Mahfouz A: A comparison of automatic cell identification methods for single-cell RNA sequencing data. Genome Biol 2019, 20(1):194-194.
3.	Pliner HA, Shendure J: Supervised classification enables rapid annotation of cell atlases. 2019, 16(10):983-986.
4.	Ma F, Pellegrini M: ACTINN: automated identification of cell types in single cell RNA sequencing. Bioinformatics (Oxford, England) 2020, 36(2):533-538.
5.	Alquicira-Hernandez J, Sathe A, Ji HP, Nguyen Q, Powell JE: scPred: accurate supervised method for cell-type classification from single-cell RNA-seq data. Genome Biol 2019, 20(1):264.
6.	Kim T, Lo K, Geddes TA, Kim HJ, Yang JYH, Yang P: scReClassify: post hoc cell type classification of single-cell rNA-seq data. 2019, 20(Suppl 9):913.
7.	Ceder JA, Jansson L, Helczynski L, Abrahamsson PA: Delta-like 1 (Dlk-1), a novel marker of prostate basal and candidate epithelial stem cells, is downregulated by notch signalling in intermediate/transit amplifying cells of the human prostate. European urology 2008, 54(6):1344-1353.
8.	Ma S, Chan KW, Hu L, Lee TK, Wo JY, Ng IO, Zheng BJ, Guan XY: Identification and characterization of tumorigenic liver cancer stem/progenitor cells. Gastroenterology 2007, 132(7):2542-2556.
9.	Zhang X, Lan Y, Xu J, Quan F, Zhao E, Deng C, Luo T, Xu L, Liao G, Yan M et al: CellMarker: a manually curated resource of cell markers in human and mouse. Nucleic acids research 2019, 47(D1):D721-d728.
10.	Costea DE, Hills A, Osman AH, Thurlow J, Kalna G, Huang X, Pena Murillo C, Parajuli H, Suliman S, Kulasekara KK et al: Identification of two distinct carcinoma-associated fibroblast subtypes with differential tumor-promoting abilities in oral squamous cell carcinoma. Cancer research 2013, 73(13):3888-3901.
11.	Navab R, Strumpf D, Bandarchi B, Zhu CQ, Pintilie M, Ramnarine VR, Ibrahimov E, Radulovich N, Leung L, Barczyk M et al: Prognostic gene-expression signature of carcinoma-associated fibroblasts in non-small cell lung cancer. Proceedings of the National Academy of Sciences of the United States of America 2011, 108(17):7160-7165.
12.	Zhang Q, He Y, Luo N, Patel SJ, Han Y, Gao R, Modak M, Carotta S, Haslinger C, Kind D et al: Landscape and Dynamics of Single Immune Cells in Hepatocellular Carcinoma. Cell 2019, 179(4):829-845.e820.
13.	Zheng GXY, Terry JM, Belgrader P, Ryvkin P, Bent ZW, Wilson R, Ziraldo SB, Wheeler TD, McDermott GP, Zhu J et al: Massively parallel digital transcriptional profiling of single cells. Nature Communications 2017, 8(1):14049.
14.	Satija R, Farrell JA, Gennert D, Schier AF, Regev A: Spatial reconstruction of single-cell gene expression data. Nature biotechnology 2015, 33(5):495-502.
15.	Hoshida Y: Nearest template prediction: a single-sample-based flexible class prediction with confidence assessment. PloS one 2010, 5(11):e15543.
16.	Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES et al: Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences of the United States of America 2005, 102(43):15545-15550.
17.	Patel AP, Tirosh I, Trombetta JJ, Shalek AK, Gillespie SM, Wakimoto H, Cahill DP, Nahed BV, Curry WT, Martuza RL et al: Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma. Science (New York, NY) 2014, 344(6190):1396-1401.
18.	Puram S, Tirosh I, Parikh A, Patel A, Yizhak K, Gillespie S, Rodman C, Luo C, Mroz E, Emerick K et al: Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor Ecosystems in Head and Neck Cancer. Cell 2017, 171.


