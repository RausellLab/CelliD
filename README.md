# CelliD v0.99
R package for gene signature extraction and cell identity recognition at individual cell level from single-cell RNA-seq.

<img src=blob/legacy/tools/sticker.png height="100">

----------------------------------------

Welcome to the official Github repository of the **CelliD** software presented at the BioRxiv preprint [CelliD: gene signature extraction and cell identity recognition at individual cell level. Cortal A, Martignetti L, Six E, Rausell A. BioRxiv 2020](https://www.biorxiv.org/content/10.1101/2020.07.23.215525v1)

## Overview

CelliD is a robust statistical method that performs gene signature extraction and functional annotation for each individual cell in a single-cell RNA-seq dataset. CelliD is based on Multiple Correspondence Analysis (MCA) and produces a simultaneous representation of cells and genes in a low dimension space. Genes are then ranked by their distance to each individual cell, providing unbiased per-cell gene signatures. Such signatures proved valuable to (i) correctly predict cell type labels at individual cell resolution, (ii) correctly match cells from the same cell type across independent datasets, overcoming batch effects arising from different technologies, tissues-of-origin and donors, and (iii) uncover functionally relevant cell heterogeneity that would have been missed by clustering-based approaches. CelliD enables the robust identification of rare or even unique cells whose gene signatures are reproducible across diverse single-cell omics datasets. 

----------------------------------------

## Installation

CelliD has recently moved to [Bioconductor](https://bioconductor.org/packages/devel/bioc/html/CelliD.html) but is still in the devel branch and therefore can be installed only with R 4.1. The master branch of this repository is the mirror of the bioconductor package.
In order to use CelliD with R version 3.6 ~ 4.0, please install the legacy branch of the repository. 

Within R, set first:
```r
install.packages("devtools")
setRepositories(ind = c(1,2,3))
```
To install CelliD then just type:
```r
devtools::install_github("RausellLab/CelliD", ref = "legacy")
```
## Known installation issues

MAC OS users might experience installation issues related to Gfortran library. To solve such issue download and install the appropriate gfortran dmg file from https://github.com/fxcoudert/gfortran-for-macOS

## Data input formats

CelliD use as input single cell data in the form of specific S4 objects. Currently supported files are SingleCellExperiment from Bioconductor and Seurat Version 3 or 4 from CRAN.

## Vignettes

A vignette illustrating CelliD step-by-step procedures is provided [here](https://rauselllab.github.io/CelliD//vignettes/vign.html). Applications include MCA dimensionality reduction, per-cell gene signatures extraction, automatic cell type prediction using marker gene lists, label-transferring across datasets and functional enrichment analysis.

## Authors

* **Akira Cortal** - [akira.cortal@institutimagine.org](akira.cortal@institutimagine.org)
* **Antonio Rausell** -  [antonio.rausell@institutimagine.org](antonio.rausell@institutimagine.org)


## License

This project is licensed under the GNU General Public License 3 - see the [LICENSE](LICENSE) file for details

## References

CelliD: gene signature extraction and cell identity recognition at individual cell level. Cortal A, Martignetti L, Six E, Rausell A. BioRxiv 2020 doi: [https://doi.org/10.1101/2020.07.23.215525](https://doi.org/10.1101/2020.07.23.215525)

## Companion Github repository CelliDPaperScript

Companion Github repository with R scripts and intermediate data representations required to reproduce all figures from the Cell-ID manuscript can be found here https://github.com/RausellLab/CellIDPaperScript.

## Updates 
You may follow us in Twitter for regular updates: https://twitter.com/AntonioRausell
