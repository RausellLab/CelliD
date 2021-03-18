#' Check for CelliD arguments
#'
#' @param X Seurat or SingleCell Experiment Object
#' @param group.by Name of meta.data or ColData column.
#' @param reduction Which dimensionality reduction to use, must be based on MCA.
#' @param dims A vector of integers indicating which dimensions to use of specified reduction embeddings and loadings.
#' @param features Character vector of feature names to subset feature coordinates. If not specified will take all features available from specified reduction loadings.
#' @param cells Character vector of cell names to subset cell coordinates. If not specified will take all features available from specified reduction Embeddigns.
#'
#' @return  list of corrected arguments if no error is thrown.
checkCelliDArg <-
    function(X, group.by, reduction, dims, features, cells) {
        UseMethod("checkCelliDArg", X)
    }
#' @rdname checkCelliDArg
checkCelliDArg.Seurat <-
    function(X, group.by = NULL, reduction, dims, features = NULL, cells = NULL) {
        # dimension check ---------------------------------------------------------
        if (length(reduction) != 1) {
            stop("Length reduction > 1, please specify a single reduction name")
        }
        if (is.null(X@reductions[[reduction]])) {
            stop(glue("specified reduction {reduction} does not exist"))
        }
        if (is.null(X@reductions[[reduction]]@misc$mca.flag)) {
            stop("selected reduction is not an mca")
        }
        if (length(dims) == 1) {
            dims <- seq(dims)
        }
        if (!all(dims %in% seq(ncol(Embeddings(X, reduction))))) {
            stop(
                glue(
                    "specified dimension not covered in {reduction}, only {ncol(Embeddings(X, reduction))} dimension computed"
                )
            )
        }
        # feature check -----------------------------------------------------------
        if (is.null(features)) {
            features <- rownames(Loadings(X, reduction))
        }
        else {
            Tfeatures <- rownames(Loadings(X, reduction))
            featureCheck <- features %in% Tfeatures
            if (!all(featureCheck)) {
                missingFeatures <- features[!featureCheck]
                matchingFeatures <- features[featureCheck]
                message(
                    glue(
                        "\nfeatures {paste(missingFeatures, collapse = \",\")} not in chosen reduction {reduction}\n"
                    )
                )
                featuresMenu <-
                    menu(
                        choices = c("continue without", "cancel"),
                        graphics = FALSE
                    )
                if (featuresMenu == 1) {
                    features <- matchingFeatures
                }
                else {
                    stop("function aborted by user")
                }
            }
        }
        # cells check -------------------------------------------------------------
        if (is.null(cells)) {
            cells <- rownames(Embeddings(X, reduction))
        } else {
            Tcells <- rownames(Embeddings(X, reduction))
            cellCheck <- cells %in% Tcells
            if (!all(cellCheck)) {
                missingcells <- cells[!cellCheck]
                matchingcells <- cells[cellCheck]
                message(glue(
                    "\ncells {missingcells} not in reduction {reduction}"
                ))
                cellsMenu <-
                    menu(
                        choices = c("continue without", "cancel"),
                        graphics = FALSE
                    )
                if (cellsMenu == 1) {
                    cells <- matchingcells
                }
                else {
                    stop("function aborted by user")
                }
            }
        }

        # check metadata ----------------------------------------------------------

        if (is.null(group.by)) {
            group.by.vec <- X@active.ident
        } else {
            if (group.by %in% colnames(X@meta.data)) {
                group.by.vec <- X@meta.data[[group.by]]
                names(group.by.vec) <- rownames(X@meta.data)
            }
            else {
                stop(glue(
                    "{group.by} column not present in seurat object meta.data"
                ))
            }
        }
        fc <-
            list(
                dims = dims,
                features = features,
                cells = cells,
                group.by = group.by,
                group.by.vec = group.by.vec
            )
        return(fc)
    }

#' @rdname checkCelliDArg
checkCelliDArg.SingleCellExperiment <-
    function(X, reduction, dims, features = NULL, cells = NULL, group.by = NULL) {
        # dimension check ---------------------------------------------------------
        if (length(reduction) != 1) {
            stop("Please specify a singe reduction name")
        }
        dimSlot <- reducedDim(X, reduction)
        cellEmb <- dimSlot
        featureEmb <- attr(dimSlot, "genesCoordinates")
        if (is.null(dimSlot)) {
            stop(glue("specified reduction {reduction} does not exist"))
        }
        if (is.null(attr(dimSlot, "mcaFlag"))) {
            stop("selected reduction is not an mca")
        }
        if (length(dims) == 1) {
            dims <- seq(dims)
        }
        if (!all(dims %in% seq(ncol(dimSlot)))) {
            stop(
                glue(
                    "specified dimension not covered in {reduction}, only {ncol(dimSlot)} dimension computed"
                )
            )
        }
        # feature check -----------------------------------------------------------
        if (is.null(features)) {
            features <- rownames(featureEmb)
        }
        else {
            Tfeatures <- rownames(featureEmb)
            featureCheck <- features %in% Tfeatures
            if (!all(featureCheck)) {
                missingFeatures <- features[!featureCheck]
                matchingFeatures <- features[featureCheck]
                message(
                    glue(
                        "\nfeatures {paste(missingFeatures, collapse = \",\")} not in chosen reduction {reduction}\n"
                    )
                )
                featuresMenu <-
                    menu(
                        choices = c("continue without", "cancel"),
                        graphics = FALSE
                    )
                if (featuresMenu == 1) {
                    features <- matchingFeatures
                }
                else {
                    stop("function aborted by user")
                }
            }
        }
        # cells check -------------------------------------------------------------
        if (is.null(cells)) {
            cells <- rownames(cellEmb)
        } else {
            Tcells <- rownames(cellEmb)
            cellCheck <- cells %in% Tcells
            if (!all(cellCheck)) {
                missingcells <- cells[!cellCheck]
                matchingcells <- cells[cellCheck]
                message(glue(
                    "\ncells {missingcells} not in reduction {reduction}"
                ))
                cellsMenu <-
                    menu(
                        choices = c("continue without", "cancel"),
                        graphics = FALSE
                    )
                if (cellsMenu == 1) {
                    cells <- matchingcells
                }
                else {
                    stop("function aborted by user")
                }
            }
        }

        # check metadata ----------------------------------------------------------
        if (is.null(group.by)) {
            group.by <- NULL
            group.by.vec <- NULL
        } else {
            if (group.by %in% colnames(colData(X))) {
                group.by.vec <- colData(X)[[group.by]]
                names(group.by.vec) <- rownames(colData(X))
            }
            else {
                stop(
                    glue(
                        "{group.by} column not present in SingleCellExperiment object colData"
                    )
                )
            }
        }
        fc <-
            list(
                dims = dims, features = features, cells = cells, group.by.vec = group.by.vec, group.by = group.by
            )
        return(fc)
    }
