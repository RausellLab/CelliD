#' Hallmark Pathways from MSigDB
#'
#' A dataset containing the Hallmark gene sets from MSigDB.
#' @format A named list of length 50 containing Hallmark gene sets.
#' @references Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst. 2015 Dec 23;1(6):417-425.
#' @source \url{http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.2/h.all.v6.2.symbols.gmt}
"Hallmark"

#' Seurat object of 400 PBMC cells
#'
#' A subset of the PBMC3k data from Seurat vignette. Normalisation, VariableFeatures, ScaleData and PCA has alreay been computed with default Seurat parameter.
#' @format A seurat object.
#' @references Butler et al., Nature Biotechnology 2018.
#' @source \url{https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz}
"seuratPbmc"

#' Mus Musculus Protein Coding Genes
#'
#' A gene list of mouse protein coding genes extracted from biomaRt.
#' @format A list of 3857 gene onthology terms with the corresponding genes.
#' @references The Gene Ontology project in 2008, The Gene Ontology Consortium Nucleic Acids Research, Volume 36, Issue suppl_1, January 2008, Pages D440–D444,
#' @source \url{http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C5}
"MgProteinCodingGenes"


#' Homo Sapiens Protein Coding Genes
#'
#' A gene list of human protein coding genes extracted from biomaRt.
#' @format A list of 19308 gene onthology terms with the corresponding genes.
#' @references The Gene Ontology project in 2008, The Gene Ontology Consortium Nucleic Acids Research, Volume 36, Issue suppl_1, January 2008, Pages D440–D444,
#' @source \url{http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C5}
"HgProteinCodingGenes"
