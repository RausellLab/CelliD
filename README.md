# CelliD
Gene signature extraction and cell-identity recognition at individual cell level from single-cell RNA-seq. 

Cell-ID v0.99


----------------------------------------

Welcome to the official Github repository of the Cell-ID software.

**Cell-ID enables the identification of rare individual cells and their reproducible gene signatures across independent single-cell RNA-seq datasets.**

Cell-ID has been recently presented at the America Society of Human Genetics annual meeting (ASHG 2019) as a platform talk with the following abstact: https://eventpilotadmin.com/web/page.php?page=Session&project=ASHG19&id=204004

Cell-ID software will be officially released in this repository as an open-source R package in the next weeks. 

If you would like to be informed about the release please send an email to antonio.rausell@institutimagine.org with subject: “Cell-ID release announcement”. 

You may also follow up in Twitter for regular updates: https://twitter.com/AntonioRausell

----------------------------------------

# Abstract

Akira Cortal 1, Antonio Rausell 1,2

1. Paris Descartes University-Sorbonne Paris Cité, Imagine Institute, Clinical Bioinformatics Lab, Paris, France; 
2. INSERM UMR 1163, Institut Imagine, Paris, France

Single-cell transcriptome profiling of patient’s biological samples may help identifying abnormal rare cell fractions potentially associated to disease. Nonetheless, the computational identification of bona-fide rare cells, representing <2%, is challenged by high levels of biological and technical noise. Current computational methods for single-cell data analysis often rely on a low-dimensional representation of cells. The characterization of cell types is then typically carried out through a clustering step followed by differential gene expression analysis among groups. However, such approach is sensitive to batch effects that may compromise the replication of findings across independent donors. Moreover, an exhaustive exploration of cellular heterogeneity requires a per-cell gene signature assessment rather than a group-based analysis. Such possibility was lacking in the scientific literature.

Here we present Cell-ID, a robust statistical method that performs gene signature extraction and functional annotation for each individual cell in a single-cell RNA-seq dataset. Cell-ID is based on Multiple Correspondence Analysis and produces a simultaneous representation of cells and genes in a low dimension space. Genes are then ranked by their distance to each individual cell providing unbiased per-cell gene signatures. Such signatures proved valuable to estimate cell similarities across independent datasets and overcomed batch effects arising from different technologies, tissues-of-origin and donors. We evaluated Cell-ID on a diverse collection of single-cell RNA-seq libraries including blood cells, airway epithelial cells, as well as pancreatic islet cells from both healthy donors and Type 2 diabetes patients. Cell-ID correctly predicted well-established rare cell types at individual cell resolution. We then showed that the unbiased per-cell gene signatures obtained by Cell-ID allowed the identification of the corresponding cells of the same type both within and across independent datasets. Finally, we applied Cell-ID to pancreatic islet single-cell RNA-seq data where we illustrated the ability of our per-cell signature analysis to uncover functionally relevant cell heterogeneity that would have been missed by a clustering-based approach. Overall, we demonstrate that Cell-ID enables the robust identification of rare or even unique cells with a potential role in human disease. Cell-ID is freely distributed as an R package with a user-friendly shiny interface.



## Contact
Please address comments and questions about Cell-ID to:
* **Akira Cortal** - *Package Maintainer* - [akira.cortal@institutimagine.org](akira.cortal@institutimagine.org)
* **Antonio Rasuell** - *Supervisor* - [antonio.rausell@institutimagine.org](antonio.rausell@institutimagine.org)

