##  ............................................................................
##  MC Dimensionality Reduction                                             ####

#' Diffusion Map on MCA coordinates
#' 
#' @description (!EXPERIMENTAL) Run DiffusionMap on MCA cell and feature coordinates. 
#' This will allow to draw the trajectory of both cells and the genes at the same time.  
#'
#' @param X Seurat or SingleCellExperiment object
#' @param reduction Which dimensionality reduction to use, must be based on MCA.
#' @param features Character vector of feature names to subset feature coordinates. If not specified will take all features available from specified reduction Loadings.
#' @param dims A vector of integers indicating which dimensions to use with reduction embeddings and loadings for distance calculation.
#' @param reduction.name name of the created dimensionlaity reduction, default set to "mca" for Seurat and "MCA" for SCE.
#' @param ... other arguments passed to methods or DiffusionMap
#'
#' @return Seurat or SingleCellExperiment object with MCDMAP stored in the reduction slot
#' @export
#'
#' @examples
#' seuratPbmc <- RunMCA(seuratPbmc, nmcs = 5)
#' seuratPbmc <- RunMCDMAP(seuratPbmc, dims = seq(5), k = 5)
RunMCDMAP <-
    function(X, reduction, features, dims, reduction.name, ...) {
        UseMethod("RunMCDMAP", X)
    }

#' @rdname RunMCDMAP
#' @param assay Seurat Asssay slot name.
#' @export
RunMCDMAP.Seurat <-
    function(X, reduction = "mca", features = NULL, dims = seq(50), reduction.name = "mcdmap", assay = DefaultAssay(X), ...) {
        requireNamespace("destiny", quietly = TRUE)
        GeneCellCoordinates <-
            GetGeneCellCoordinates(
                X = X,
                reduction = reduction,
                dims = dims,
                features = features
            )
        if (any(duplicated(GeneCellCoordinates))) {
            GeneCellCoordinates[duplicated(GeneCellCoordinates), ncol(GeneCellCoordinates)] <- GeneCellCoordinates[duplicated(GeneCellCoordinates), ncol(GeneCellCoordinates)] + runif(min = 10^-6, max = 10^-5, n = sum(duplicated(GeneCellCoordinates)))
        }
        MCDMAP <-
            destiny::DiffusionMap(data = GeneCellCoordinates, suppress_dpt = TRUE, ...)
        Emb <- MCDMAP@eigenvectors
        rownames(Emb) <- rownames(GeneCellCoordinates)
        cellEmb <- Emb[rownames(Emb) %in% rownames(Embeddings(X, reduction)), ]
        geneEmb <- Emb[!rownames(Emb) %in% rownames(Embeddings(X, reduction)), ]
        X <-
            setDimMCSlot(
                X = X,
                cellEmb = cellEmb,
                geneEmb = geneEmb,
                assay = assay,
                reduction.name = reduction.name
            )
        return(X)
    }

#' @rdname RunMCDMAP
#' @export
RunMCDMAP.SingleCellExperiment <-
    function(X, reduction = "MCA", features = NULL, dims = seq(50), reduction.name = "MCDMAP", ...) {
        requireNamespace("destiny", quietly = TRUE)
        GeneCellCoordinates <-
            GetGeneCellCoordinates(
                X = X,
                reduction = reduction,
                dims = dims,
                features = features
            )
        MCDMAP <-
            destiny::DiffusionMap(data = GeneCellCoordinates, ...)
        Emb <- MCDMAP@eigenvectors
        rownames(Emb) <- rownames(GeneCellCoordinates)
        geneEmb <- Emb[!rownames(Emb) %in% rownames(reducedDim(X, reduction)), ]
        cellEmb <- Emb[rownames(Emb) %in% rownames(reducedDim(X, reduction)), ]
        X <-
            setDimMCSlot(
                X = X,
                cellEmb = cellEmb,
                geneEmb = geneEmb,
                reduction.name = reduction.name
            )
        return(X)
    }

#' tSNE on MCA coordinates
#' 
#' @description (!EXPERIMENTAL) Run TSNE on MCA fetures and cells coordinates
#' This will allow to embbed in 2D both cells and the genes at the same time. 
#'
#' @param X Seurat or SingleCellExperiment object
#' @param reduction Which dimensionality reduction to use, must be based on MCA.
#' @param features Character vector of feature names to subset feature coordinates. If not specified will take all features available from specified reduction Loadings.
#' @param dims A vector of integers indicating which dimensions to use with reduction embeddings and loadings for distance calculation.
#' @param reduction.name name of the created dimensionlaity reduction, default set to "mca" for Seurat and "MCA" for SCE.
#' @param ... other arguments passed to methods or Rtsne::Rtsne
#'
#' @return Seurat or SingleCellExperiment object with MCTSNE stored in the reduction slot
#' @importFrom Rtsne Rtsne
#' @export
#'
#' @examples
#' seuratPbmc <- RunMCA(seuratPbmc, nmcs = 5)
#' seuratPbmc <- RunMCTSNE(seuratPbmc, dims = seq(5))
RunMCTSNE <-
    function(X, reduction, dims, features, reduction.name, ...) {
        UseMethod("RunMCTSNE", X)
    }

#' @rdname RunMCTSNE
#' @param assay Seurat assay slot. When not specified set with DefaultAssay(X)
#' @export
RunMCTSNE.Seurat <-
    function(X, reduction = "mca", dims = seq(50), features = NULL, reduction.name = "mctsne", assay = DefaultAssay(X), ...) {
        GeneCellCoordinates <-
            GetGeneCellCoordinates(
                X = X,
                reduction = reduction,
                dims = dims,
                features = features
            )
        message("\nrunning TSNE\n")
        MCTSNE <-
            Rtsne::Rtsne(
                X = GeneCellCoordinates,
                pca = FALSE,
                check_duplicates = FALSE,
                ...
            )
        message("\nreturning seurat object\n")
        Emb <- MCTSNE$Y
        rownames(Emb) <- rownames(GeneCellCoordinates)
        geneEmb <- Emb[!rownames(Emb) %in% rownames(Embeddings(X, reduction)), ]
        cellEmb <- Emb[rownames(Emb) %in% rownames(Embeddings(X, reduction)), ]
        X <-
            setDimMCSlot(
                X = X,
                cellEmb = cellEmb,
                geneEmb = geneEmb,
                assay = assay,
                reduction.name = reduction.name
            )
        return(X)
    }

#' @rdname RunMCTSNE
#' @export
RunMCTSNE.SingleCellExperiment <-
    function(X, reduction = "MCA", dims = seq(50), features = NULL, reduction.name = "MCTSNE", ...) {
        GeneCellCoordinates <-
            GetGeneCellCoordinates(
                X = X,
                reduction = reduction,
                dims = dims,
                features = features
            )
        message("\nrunning TSNE\n")
        MCTSNE <-
            Rtsne::Rtsne(
                X = GeneCellCoordinates,
                pca = FALSE,
                check_duplicates = FALSE,
                ...
            )
        message("\nreturning Single Cell Experiment object\n")
        Emb <- MCTSNE$Y
        colnames(Emb) <- paste0(reduction.name, seq(ncol(Emb)))
        rownames(Emb) <- rownames(GeneCellCoordinates)
        geneEmb <- Emb[!rownames(Emb) %in% rownames(reducedDim(X, reduction)), ]
        cellEmb <- Emb[rownames(Emb) %in% rownames(reducedDim(X, reduction)), ]
        X <- setDimMCSlot(
            X = X,
            cellEmb = cellEmb,
            geneEmb = geneEmb,
            reduction.name = reduction.name
        )
        return(X)
    }



#' UMAP on MCA coordinates
#' 
#' @description (!EXPERIMENTAL) Run UMAP on MCA fetures and cells coordinates.
#' This will allow to embbed in 2D both cells and the genes at the same time. 
#'
#' @param X Seurat or SingleCellExperiment object
#' @param reduction Which dimensionality reduction to use, must be based on MCA.
#' @param features Character vector of feature names to subset feature coordinates. If not specified will take all features available from specified reduction Loadings.
#' @param dims A vector of integers indicating which dimensions to use with reduction embeddings and loadings for distance calculation.
#' @param reduction.name name of the created dimensionlaity reduction, default set to "mca" for Seurat and "MCA" for SCE.
#' @param ... other arguments passed to methods or Rtsne::Rtsne
#'
#' @return Seurat or SingleCellExperiment object with MCUMAP stored in the reduction slot
#' @importFrom reticulate py_module_available
#' @importFrom umap umap
#' @export
#'
#' @examples
#' seuratPbmc <- RunMCA(seuratPbmc, nmcs = 5)
#' seuratPbmc <- RunMCUMAP(seuratPbmc, dims = seq(5))
RunMCUMAP <-
    function(X, reduction, dims, features, reduction.name, ...) {
        UseMethod("RunMCUMAP", X)
    }

#' @rdname RunMCUMAP
#' @param assay Seurat assay slot to assign MCUMAP. When not specified set to DefaultAssay(X)
#' @export
RunMCUMAP.Seurat <-
    function(X, reduction = "mca", dims = seq(50), features = NULL, reduction.name = "mcumap", assay = DefaultAssay(X), ...) {
        GeneCellCoordinates <-
            GetGeneCellCoordinates(
                X = X,
                reduction = reduction,
                dims = dims,
                features = features
            )
        message("\nrunning UMAP\n")
        if (py_module_available("umap")) {
            method <- "umap-learn"
        }
        else {
            message("\numap-learn not detected\n")
            method <- "naive"
        }
        MCUMAP <- umap(d = GeneCellCoordinates, method = method, ...)
        message("\nreturning Seurat object\n")
        Emb <- MCUMAP$layout
        rownames(Emb) <- rownames(GeneCellCoordinates)
        cellEmb <- Emb[rownames(Emb) %in% rownames(Embeddings(X, reduction)), ]
        geneEmb <- Emb[!rownames(Emb) %in% rownames(Embeddings(X, reduction)), ]
        X <- setDimMCSlot(
            X = X,
            cellEmb = cellEmb,
            geneEmb = geneEmb,
            assay = assay,
            reduction.name = reduction.name
        )
        return(X)
    }

#' @rdname RunMCUMAP
#' @export
RunMCUMAP.SingleCellExperiment <-
    function(X, reduction = "MCA", dims = seq(50), features = NULL, reduction.name = "MCUMAP", ...) {
        GeneCellCoordinates <-
            GetGeneCellCoordinates(
                X = X,
                reduction = reduction,
                dims = dims,
                features = features
            )
        message("\nrunning UMAP\n")
        if (reticulate::py_module_available("umap")) {
            method <- "umap-learn"
        }
        else {
            message("\numap-learn not detected\n")
            method <- "naive"
        }
        MCUMAP <- umap(d = GeneCellCoordinates, method = method, ...)
        message("\nreturning Single Cell Experiment object\n")
        Emb <- MCUMAP$layout
        rownames(Emb) <- rownames(GeneCellCoordinates)
        cellEmb <- Emb[rownames(Emb) %in% rownames(reducedDim(X, reduction)), ]
        geneEmb <- Emb[!rownames(Emb) %in% rownames(reducedDim(X, reduction)), ]
        X <-
            setDimMCSlot(
                X = X,
                cellEmb = cellEmb,
                geneEmb = geneEmb,
                reduction.name = reduction.name
            )
        return(X)
    }


#' GeneCellCoordinates
#' 
#' @description Get coordinates of both cells and features in a matrix
#'
#' @param X Seurat or SingleCellExperiment Object
#' @param reduction Which dimensionality reduction to use, must be based on MCA.
#' @param dims A vector of integers indicating which dimensions to use with reduction embeddings and loadings for distance calculation.
#' @param features Character vector of feature names to subset feature coordinates. If not specified will take all features available from specified reduction Loadings.
#' @importFrom stats runif
#'
#' @return A matrix with gene and cell coordinates of MCA
GetGeneCellCoordinates <- function(X, reduction, dims, features) {
    UseMethod("GetGeneCellCoordinates", X)
}


GetGeneCellCoordinates.Seurat <-
    function(X, reduction, dims, features) {
        message("\ngetting feature and cell coordinates\n")
        check <-
            checkCelliDArg(X,
                reduction = reduction,
                dims = dims,
                features = features
            )
        features <- check$features
        cells <- check$cells
        dims <- check$dims
        GeneCoordinates <- Loadings(X, reduction)[features, dims]
        CellCoordinates <- Embeddings(X, reduction)[cells, dims]
        GeneCellCoordinates <-
            rbind(GeneCoordinates, CellCoordinates)
        return(GeneCellCoordinates)
    }

GetGeneCellCoordinates.SingleCellExperiment <-
    function(X, reduction, dims, features) {
        message("\ngetting feature and cell coordinates\n")
        check <-
            checkCelliDArg(X,
                reduction = reduction,
                dims = dims,
                features = features
            )
        features <- check$features
        cells <- check$cells
        dims <- check$dims
        GeneCoordinates <-
            attr(reducedDim(X, reduction), "genesCoordinates")[features, dims]
        CellCoordinates <- reducedDim(X, reduction)[, dims]
        GeneCellCoordinates <-
            rbind(GeneCoordinates, CellCoordinates)
        return(GeneCellCoordinates)
    }
