##  ............................................................................
##  Group Coordinates                                                       ####

#' Centroids calculation for a given group
#'
#' @param X  Seurat or SingleCellExperiment object, alternatively a matrix.
#' @param group.by  column name of meta.data (Seurat) or ColData (SingleCellExperiment). For Seurat object if NULL active.ident slot will be taken.
#' @param reduction Which dimensionality reduction to use, must be based on MCA.
#' @param dims A vector of integers indicating which dimensions to use with reduction embeddings and loadings for distance calculation.
#' @param ... Other arguments passed to methods
#'
#' @return A data.table with coordinates of the group centroids for the specidied dims.
GetGroupCoordinates <- function(X, group.by, reduction, dims, ...) {
    UseMethod("GetGroupCoordinates", X)
}

#' @rdname GetGroupCoordinates
#' @export
GetGroupCoordinates.matrix <- function(X, group.by, reduction = NULL, dims, ...) {
    group <- NULL
    group.by <- group.by[rownames(X)]
    DT_cells <- as.data.table(X[, dims], keep.rownames = "cells")
    DT_cells[, group := group.by]
    DT_cells <- DT_cells[order(group), ]
    setcolorder(DT_cells, c(c("cells", "group"), setdiff(
        colnames(DT_cells), c("cells", "group")
    )))
    axis <- tail(names(DT_cells), -2)
    group_centroids <-
        DT_cells[, lapply(.SD, mean), by = group, .SDcols = axis]
    group_coordinates <-
        as.matrix(group_centroids, rownames = "group")
    return(group_coordinates)
}


#' @rdname GetGroupCoordinates
#' @export
GetGroupCoordinates.Seurat <-
    function(X, group.by = NULL, reduction = "mca", dims = seq(50), ...) {
        check <-
            checkCellIDArg(
                X,
                reduction = reduction,
                dims = dims,
                features = NULL,
                group.by = group.by
            )
        dims <- check$dims
        group.by <- check$group.by.vec
        Emb <- Embeddings(X, reduction)
        group_coordinates <-
            GetGroupCoordinates(
                X = Emb,
                group.by = group.by,
                dims = dims
            )
        return(group_coordinates)
    }

#' @rdname GetGroupCoordinates
#' @export
GetGroupCoordinates.SingleCellExperiment <-
    function(X, group.by = NULL, reduction = "MCA", dims, ...) {
        check <-
            checkCellIDArg(
                X,
                reduction = reduction,
                dims = dims,
                features = NULL,
                cells = NULL,
                group.by = group.by
            )
        dims <- check$dims
        group.by <- check$group.by.vec
        Emb <- reducedDim(X, reduction)
        group_coordinates <-
            GetGroupCoordinates(
                X = Emb,
                group.by = group.by,
                dims = dims
            )
        return(group_coordinates)
    }


##  ............................................................................
##  Group Gene Euclidean Distances                                          ####

#' Distance calculation between genes and group
#'
#' @param X  Seurat or SingleCellExperiment object, alternatively a matrix.
#' @param group.by  column name of meta.data (Seurat) or ColData (SingleCellExperiment)
#' @param reduction Which dimensionality reduction to use, must be based on MCA.
#' @param dims A vector of integers indicating which dimensions to use with reduction embeddings and loadings for distance calculation.
#' @param features A character vector of features name to subset feature coordinates for distance calculation.
#'
#' @return Distance Matrix between groups (column) and genes (row)
GetGroupGeneDistance <-
    function(X, group.by, reduction, dims, features) {
        UseMethod("GetGroupGeneDistance", X)
    }

#' @rdname GetGroupGeneDistance
#' @export
GetGroupGeneDistance.Seurat <-
    function(X, group.by = NULL, reduction = "mca", dims = seq(50), features = NULL) {
        check <-
            checkCellIDArg(
                X,
                reduction = reduction,
                dims = dims,
                features = features,
                cells = NULL,
                group.by = group.by
            )
        dims <- check$dims
        features <- check$features
        group.by <- check$group.by
        group_coordinates <-
            GetGroupCoordinates(X,
                group.by = group.by,
                reduction = reduction,
                dims = dims
            )
        if (is.null(features)) {
            genes_coordinates <- Loadings(X, reduction)[, dims]
        }
        else {
            genes_coordinates <- Loadings(X, reduction)[features, dims]
        }
        GroupGeneDistance <- t(pairDist(genes_coordinates, group_coordinates))
        return(GroupGeneDistance)
    }

#' @rdname GetGroupGeneDistance
#' @export
GetGroupGeneDistance.SingleCellExperiment <-
    function(X, group.by, reduction = "MCA", dims = seq(50), features = NULL) {
        check <-
            checkCellIDArg(
                X,
                reduction = reduction,
                dims = dims,
                features = features,
                cells = NULL,
                group.by = group.by
            )
        dims <- check$dims
        features <- check$features
        group.by <- check$group.by
        group_coordinates <-
            GetGroupCoordinates(X, group.by, reduction = reduction, dims = dims)
        if (is.null(features)) {
            genes_coordinates <-
                attr(reducedDim(X, reduction), "genesCoordinates")[, dims]
        }
        else {
            genes_coordinates <-
                attr(reducedDim(X, reduction), "genesCoordinates")[features, dims]
        }
        GroupGeneDistance <-t(pairDist(genes_coordinates, group_coordinates))
        return(GroupGeneDistance)
    }

##  ............................................................................
##  Group Gene Ranking                                                      ####


#' Gene Specificity Ranking Calculation
#'
#' @param X  Seurat or SingleCellExperiment object, alternatively a matrix.
#' @param group.by  column name of meta.data (Seurat) or ColData (SingleCellExperiment)
#' @param reduction Which dimensionality reduction to use, must be based on MCA.
#' @param dims A vector of integers indicating which dimensions to use with reduction embeddings and loadings for distance calculation.
#' @param features A character vector of features name to subset feature coordinates for distance calculation.
#'
#' @return List of genes ranking for each groups
#' @export
#'
#' @examples
#' seuratPbmc <- RunMCA(seuratPbmc, nmcs = 5)
#' GroupGeneRanking <- GetGroupGeneRanking(seuratPbmc, group.by = "seurat_clusters", dims = 1:5)
GetGroupGeneRanking <-
    function(X, group.by, reduction, dims, features) {
        UseMethod("GetGroupGeneRanking", X)
    }


#' @rdname GetGroupGeneRanking
#' @export
GetGroupGeneRanking.Seurat <-
    function(X, group.by = NULL, reduction = "mca", dims = seq(50), features = NULL) {
        GroupGeneDistance <-
            GetGroupGeneDistance(
                X = X,
                group.by = group.by,
                dims = dims,
                reduction = reduction,
                features = features
            )
        GroupGeneRanking <- DistSort(GroupGeneDistance)
        return(GroupGeneRanking)
    }

#' @rdname GetGroupGeneRanking
#' @export
GetGroupGeneRanking.SingleCellExperiment <-
    function(X, group.by, reduction = "MCA", dims = seq(50), features = NULL) {
        GroupGeneDistance <-
            GetGroupGeneDistance(
                X = X,
                group.by = group.by,
                dims = dims,
                reduction = reduction,
                features = features
            )
        GroupGeneRanking <- DistSort(GroupGeneDistance)
        return(GroupGeneRanking)
    }


##  ............................................................................
##  Group Gene Set                                                          ####


#' Extract cluster/group gene sets from MCA
#'
#' @param X  Seurat or SingleCellExperiment object, alternatively a matrix.
#' @param group.by  column name of meta.data (Seurat) or ColData (SingleCellExperiment).
#' @param reduction Which dimensionality reduction to use, must be based on MCA.
#' @param dims A vector of integers indicating which dimensions to use with reduction for distance calculation.
#' @param features A character vector of features name to subset feature coordinates for distance calculation.
#' @param n.features A single integer specifying how many top features will be extracted from ranking.
#' @return Distance Matrix between groups (column) and genes (row)
#' @export
#' @examples
#' seuratPbmc <- RunMCA(seuratPbmc, nmcs = 5)
#' GroupGeneSet <- GetGroupGeneSet(seuratPbmc, dims = 1:5, group.by = "seurat_clusters")
GetGroupGeneSet <-
    function(X, group.by, reduction, dims, features, n.features) {
        UseMethod("GetGroupGeneSet", X)
    }

#' @rdname GetGroupGeneSet
#' @export
GetGroupGeneSet.Seurat <-
    function(X,group.by = NULL,reduction = "mca",dims = seq(50),features = NULL,n.features = 200) {
        GroupGeneRanking <-
            GetGroupGeneRanking(
                X = X,
                group.by = group.by,
                dims = dims,
                reduction = reduction,
                features = features
            )
        GroupGeneSet <-
            lapply(GroupGeneRanking, function(x)
                names(head(x, n.features)))
        return(GroupGeneSet)
    }

#' @rdname GetGroupGeneSet
#' @export
GetGroupGeneSet.SingleCellExperiment <-
    function(X, group.by = NULL, reduction = "MCA", dims = seq(50), features = NULL, n.features = 200) {
        GroupGeneRanking <-
            GetGroupGeneRanking(
                X = X,
                group.by = group.by,
                dims = dims,
                reduction = reduction,
                features = features
            )
        GroupGeneSet <-
            lapply(GroupGeneRanking, function(x)
                names(head(x, n.features)))
        return(GroupGeneSet)
    }

# misc --------------------------------------------------------------------

#' Sort Gene Cell Distance Matrix
#'
#' @param distance distance matrix with features at rows and cell at columns
#'
#' @return list of ranking of genes by cells
DistSort <- function(distance) {
    message("\ncreating ranking\n")
    GroupGeneDistance_DF <- as.data.table(distance)
    GroupGeneRanking <-
        pbapply::pblapply(
            GroupGeneDistance_DF,
            FUN = function(x) {
                names(x) <- rownames(distance)
                x <- sort(x, method = "quick")
                return(x)
            }
        )
    return(GroupGeneRanking)
}
