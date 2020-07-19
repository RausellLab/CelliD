#' Seurat DimPlot for MCA like Dimensionality Reduction
#'
#' Small modification of the regular Seurat DimPlot function to enable plotting features for mca like dimensionality reduction.
#'
#' @param X a seurat object
#' @param reduction Which dimensionality reduction to use. If not specified, searches for mca.
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param features character vector of features to plot, must be present in the specified dimension loadings
#' @param size.feature integer indicating size of geom_point for features
#' @param size.feature.text integer indicating size of geom_text for features
#' @param as.text logical indicating as to include text label for feature plotting, will produce warning if TRUE and length(features) > 50
#' @param ... Other arguments passed to DimPlot
#'
#' @importFrom Seurat DimPlot
#' @importFrom ggrepel geom_text_repel
#' @return A ggplot object
#' @export
#'
#' @examples
#' seuratPbmc <- RunMCA(seuratPbmc, nmcs = 5)
#' seuratPbmc <- DimPlotMC(seuratPbmc,  features = Seurat::VariableFeatures(seuratPbmc))
DimPlotMC <-
    function(X, reduction = "mca", dims = c(1, 2), features = NULL, size.feature = 2, size.feature.text = 5, as.text = FALSE, ...) {
        check <-checkCellIDArg(
            X = X,
            dims = dims,
            reduction = reduction,
            features = features
        )
        features <- check$features
        dims <- check$dims
        CellData <- as.data.table(Embeddings(X, reduction)[, dims])
        featureData <- as.data.table(Loadings(X, reduction)[features, dims])
        if (length(features) == 1) {
            featureData <- as.data.frame(t(featureData))
        }
        featureData$features <- features
        featureData$genes <- "black"
        DIM1 <- colnames(featureData)[1]
        DIM2 <- colnames(featureData)[2]
        MCPlot <-
            DimPlot(X, dims = dims, reduction = reduction, ...) + geom_point(
                data = featureData,
                mapping = aes_string(x = DIM1, y = DIM2, text = "features",fill = "genes"),
                size = size.feature,
                shape = 4
            ) + scale_fill_identity(name = 'genes', labels = c(""), guide = 'legend')
        if (as.text) {
            if (length(features) > 50) {
                if (menu(
                    choices = c("Yes", "No"),
                    title = "you are plotting as Text more than 50 genes, this can take a long time to plot, continue?"
                ) == 2) {
                    stop()
                }
            }
            MCPlot <-
                MCPlot + ggrepel::geom_text_repel(
                    data = featureData,
                    mapping = aes_string(
                        x = DIM1,
                        y = DIM2,
                        text = "features",
                        label = "features"
                    ),
                    size = size.feature.text,
                    seed = 1,
                    min.segment.length = 0.01,
                    point.padding = 0.2
                )
        }
        return(MCPlot)
    }


#' Scater plotReducedDim for MCA like dimensionality Reduction
#'
#' Small modification of the Scater plotReducedDim function to enable plotting features for mca like dimensionality reduction.
#'
#' @param X a Single Cell Experiment Object
#' @param reduction Which dimensionality reduction to use. If not specified, searches for mca.
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param features character vector of features to plot, must be present in the specified dimension loadings
#' @param size.feature integer indicating size of geom_point for features
#' @param size.feature.text integer indicating size of geom_text for features
#' @param as.text logical indicating as to include text label for feature plotting, will produce warning if TRUE and length(features) > 50.
#' @param ... Other arguments passed to plotReducedDim
#'
#' @return A ggplot object
#'
#' @export
#' @importFrom scater plotReducedDim
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' scePBMC <- as.SingleCellExperiment(seuratPbmc)
#' scePBMC <- RunMCA(scePBMC, nmcs = 5)
#' plotReducedDimMC(scePBMC)
plotReducedDimMC <-
    function(X, reduction = "MCA", dims = c(1,2), features = NULL, size.feature = 3, size.feature.text = 5, as.text = FALSE, ...) {
        check <-
            checkCellIDArg(
                X,
                reduction = reduction,
                dims = dims,
                features = features,
                group.by = NULL
            )
        features <- check$features
        dims <- check$dims
        featureData <-
            as.data.frame(attr(reducedDim(X, reduction), "genesCoordinates"))[features, dims]
        featureData <-
            as.data.table(featureData, keep.rownames = "features")
        DIM1 <- colnames(featureData)[2]
        DIM2 <- colnames(featureData)[3]
        MCPlot <-
            plotReducedDim(X,
                           ncomponents = dims,
                           dimred = reduction,
                           ...
            ) + geom_point(
                data = featureData,
                mapping = aes_string(x = DIM1, y = DIM2, text ="features"),
                size = size.feature,
                shape = 4
            )
        if (as.text) {
            if (length(features) > 50) {
                if (menu(
                    choices = c("Yes", "No"),
                    title = "you are plotting as Text more than 50 genes, this can take along time to plot, continue?"
                ) == 2) {
                    stop()
                }
            }
            MCPlot <-
                MCPlot + ggrepel::geom_text_repel(
                    data = featureData,
                    mapping = aes_string(
                        x = DIM1,
                        y = DIM2,
                        label = "features"
                    ),
                    size = size.feature.text,
                    point.padding = 0.2
                )
        }
        # Plotly Hover
        # MCPlot <- plotly_build(MCPlot)
        # MCPlot$x$data[[length(MCPlot$x$data)]]$text <- str_remove(MCPlot$x$data[[length(MCPlot$x$data)]]$text, "<br.*")
        return(MCPlot)
    }