#' Run Gene Set Enrichment Analysis on cells
#'
#' @description Calculate cells gene specificty ranking and then perform geneset enrichment analysis (fgsea) on it.
#' However, due to the very long running time of gene set enrichment analysis, we recommend the usage of RunCellHGT.
#'
#' @param X Seurat or SingleCellExperiment object
#' @param pathways List of gene sets to check
#' @param reduction Which dimensionality reduction to use, must be based on MCA.
#' @param dims A vector of integers indicating which dimensions to use with reduction embeddings and loadings for distance calculation.
#' @param features Character vector of feature names to subset feature coordinates. If not specified will take all features available from specified reduction Loadings.
#' @param cells Character vector of cell names to subset cell coordinates. If not specified will take all features available from specified reduction Embeddings
#' @param nperm Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param gseaParam GSEA parameter value, all gene-level statis are raised to the power of 'gseaParam' before calculation of GSEA enrichment scores
#' @param n.core A single integer to specify the number of core for parallelisation.
#'
#' @return A data.table with geneset enrichment analysis statistics.
#' @export
#'
#' @examples
#' seuratPbmc <- RunMCA(seuratPbmc, nmcs = 5)
#' GSEAResults <- RunCellGSEA(seuratPbmc, Hallmark, dims = 1:5)
RunCellGSEA <- function(X, pathways, reduction, dims, features, cells, nperm, minSize, maxSize, gseaParam, n.core) {
    UseMethod("RunCellGSEA", X)
}

#' @rdname RunCellGSEA
#' @export
RunCellGSEA.Seurat <-
    function(X,
             pathways,
             reduction = "mca",
             dims = seq(50),
             features = NULL,
             cells = NULL,
             nperm = 1000,
             minSize = 10,
             maxSize = 500,
             gseaParam = 0,
             n.core = 1) {
        # Check -------------------------------------------------------------------
        check <-
            checkCelliDArg(
                X = X,
                reduction = reduction,
                dims = dims,
                features = features,
                cells = cells
            )
        features <- check$features
        dims <- check$dims
        cells <- check$cells

        # Parallel ----------------------------------------------------------------
        BPPARAM <- BiocParallel::bpparam()
        BPPARAM$workers <- n.core
        BPPARAM$tasks <- as.integer(n.core * 5)
        BPPARAM$progressbar <- TRUE

        # Ranking -----------------------------------------------------------------
        CellGeneDist <-
            GetCellGeneDistance(
                X,
                reduction = reduction,
                dims      = dims,
                features  = features,
                cells     = cells
            )
        CellGeneDist <- as.data.table(CellGeneDist, "features")
        features <- CellGeneDist$features
        message("starting GSEA")
        gseaResults <- bplapply(CellGeneDist[, -1], fgseaCelliDPar,
            pathways = pathways,
            features = features,
            nperm = nperm,
            minSize = minSize,
            maxSize = maxSize,
            BPPARAM = BPPARAM
        )
        gseaResults <- rbindlist(gseaResults, idcol = "cell")
        return(gseaResults)
    }

#' @rdname RunCellGSEA
#' @export
RunCellGSEA.SingleCellExperiment <-
    function(X,
             pathways,
             reduction = "mca",
             dims = seq(50),
             features = NULL,
             cells = NULL,
             nperm = 1000,
             minSize = 10,
             maxSize = 500,
             gseaParam = 0,
             n.core = 1) {
        # Check -------------------------------------------------------------------
        check <-
            checkCelliDArg(
                X = X,
                reduction = reduction,
                dims = dims,
                features = features,
                cells = cells
            )
        features <- check$features
        dims <- check$dims
        cells <- check$cells

        # Parallel ----------------------------------------------------------------
        BPPARAM <- BiocParallel::bpparam()
        BPPARAM$workers <- n.core
        BPPARAM$tasks <- as.integer(n.core * 5)
        BPPARAM$progressbar <- TRUE

        # Ranking -----------------------------------------------------------------
        CellGeneDist <-
            GetCellGeneDistance(
                X,
                reduction = reduction,
                dims = dims,
                features = features,
                cells = cells
            )
        CellGeneDist <- as.data.table(CellGeneDist, "features")
        features <- CellGeneDist$features
        message("starting GSEA")
        bplapply(CellGeneDist[, -1], fgseaCelliDPar,
            pathways = pathways,
            features = features,
            nperm = nperm,
            minSize = minSize,
            maxSize = maxSize,
            BPPARAM = BPPARAM
        )
        gseaResults <- rbindlist(gseaResults, idcol = "cell")
        return(gseaResults)
    }


#' Run GSEA on cluster/groups
#'
#' @description Calculate group gene specificty ranking and then perform geneset enrichment analysis on it.
#'
#' @param X pathways List of gene sets to check
#' @param pathways reduction Which dimensionality reduction to use, must be based on MCA.
#' @param group.by dims A vector of integers indicating which dimensions to use with reduction embeddings and loadings for distance calculation.
#' @param reduction features Character vector of feature names to subset feature coordinates. If not specified will take all features available from specified reduction Loadings.
#' @param dims cells Character vector of cell names to subset cell coordinates. If not specified will take all features available from specified reduction Embeddings
#' @param features cells Character vector of cell names to subset cell coordinates. If not specified will take all features available from specified reduction Embeddings
#' @param nperm nperm Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
#' @param minSize minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param gseaParam gseaParam GSEA parameter value, all gene-level statis are raised to the power of 'gseaParam' before calculation of GSEA enrichment scores
#'
#' @return A data.table with geneset enrichment analysis statistics.
#' @export
#'
#' @examples
#' seuratPbmc <- RunMCA(seuratPbmc, nmcs = 5)
#' GSEAResults <- RunGroupGSEA(seuratPbmc, Hallmark, group.by = "seurat_clusters", dims = 1:5)
RunGroupGSEA <-
    function(X, pathways, group.by, reduction, dims, features, nperm, minSize, maxSize, gseaParam) {
        UseMethod("RunGroupGSEA", X)
    }


#' @rdname RunGroupGSEA
#' @export
RunGroupGSEA.Seurat <-
    function(X, pathways, group.by = NULL, reduction = "mca", dims = seq(50), features = NULL, nperm = 1000, minSize = 10, maxSize = 500, gseaParam = 0) {
        check <-
            checkCelliDArg(
                X = X,
                group.by = group.by,
                reduction = reduction,
                dims = dims,
                features = features,
                cells = NULL
            )
        features <- check$features
        dims <- check$dims
        group.by <- check$group.by
        message("\nranking genes \n")
        # Ranking -----------------------------------------------------------------
        GroupGeneRanking <-
            GetGroupGeneRanking(
                X = X,
                group.by = group.by,
                reduction = reduction,
                dims = dims,
                features = features
            )
        GroupGeneRanking <-
            lapply(GroupGeneRanking, function(i) {
                  (i + 0.1)^-1
              })
        gseaResultsAll <- lapply(GroupGeneRanking, function(x,
                                                            pathways,
                                                            nperm,
                                                            minSize,
                                                            maxSize) {
            stats <- x
            gseaResults <- CelliD::fgseaCelliD(
                pathways = pathways,
                stats = stats,
                nperm = nperm,
                minSize = minSize,
                maxSize = maxSize
            )
        },
        pathways = pathways,
        nperm = nperm,
        minSize = minSize,
        maxSize = maxSize
        )
        gseaResultsAll <-
            rbindlist(gseaResultsAll, idcol = "cell")
        setnames(gseaResultsAll, "cell", "group")
        return(gseaResultsAll)
    }

#' @rdname RunGroupGSEA
#' @export
RunGroupGSEA.SingleCellExperiment <-
    function(X, pathways, group.by, reduction = "MCA", dims = seq(50), features = NULL, nperm = 1000, minSize = 10, maxSize = 500, gseaParam = 0) {
        check <-
            checkCelliDArg(
                X = X,
                group.by = group.by,
                reduction = reduction,
                dims = dims,
                features = features,
                cells = NULL
            )
        features <- check$features
        dims <- check$dims
        group.by <- check$group.by
        message("\nranking genes \n")
        # Ranking -----------------------------------------------------------------
        GroupGeneRanking <-
            GetGroupGeneRanking(
                X = X,
                group.by = group.by,
                reduction = reduction,
                dims = dims,
                features = features
            )
        GroupGeneRanking <-
            lapply(GroupGeneRanking, function(i) {
                  (i + 0.1)^-1
              })
        gseaResultsAll <- lapply(GroupGeneRanking, function(x,
                                                            pathways,
                                                            features,
                                                            nperm,
                                                            minSize,
                                                            maxSize) {
            stats <- x
            gseaResults <- CelliD::fgseaCelliD(
                pathways = pathways,
                stats = stats,
                nperm = nperm,
                minSize = minSize,
                maxSize = maxSize,
                .progress = TRUE
            )
        })
        gseaResultsAll <-
            rbindlist(gseaResultsAll, idcol = "cell")
        setnames(gseaResultsAll, "cell", "group")
        return(gseaResultsAll)
    }


#' Get Matrix from Enrichment Results
#' 
#' Extract enrcihment score Matrix from RunGSEA functions.
#'
#' @param X an enrichment results obtained by RunGroupGSEA or RunCellGSEA
#' @param metric a character indicating which metric to use as value of matrix (ES, NES, padj, pval)
#'
#' @return A matrix of geneset enrichment metric with cell/group at columns and pathways/genesets at rows
#' @export
#'
#' @examples
#' seuratPbmc <- RunMCA(seuratPbmc, nmcs = 5)
#' GSEAResults <- RunGroupGSEA(seuratPbmc, Hallmark, group.by = "seurat_clusters", dims = 1:5)
#' GSEAMatrix <- GetGSEAMatrix(GSEAResults)
GetGSEAMatrix <- function(X, metric = "ES") {
    names(X)[1] <- "cell"
    dt <- dcast(X, pathway ~ cell, value.var = metric)
    gsea_mat <- as.matrix(dt[, -1])
    rownames(gsea_mat) <- dt$pathway
    return(gsea_mat)
}


# fgsea -------------------------------------------------------------------

# Extracted from fgsea
calcGseaStatCumulativeBatch <-
    getFromNamespace("calcGseaStatCumulativeBatch", "fgsea")

#' Slight change in fgsea for ram and speed efficiency in CelliD
#'
#' @param pathways List of gene sets to check
#' @param stats Named vector of gene-level stats. Names should be the same as in 'pathways'
#' @param nperm Number of permutations to do. Minimal possible nominal p-value is about 1/nperm
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param gseaParam GSEA parameter value, all gene-level stats are raised to the power of 'gseaParam' before calculation of GSEA enrichment scores
#'
#' @return A table with GSEA results. Each row corresponds to a tested pathway.
#' The columns are the following:
#' \itemize{
#'   \item pathway -- name of the pathway as in `names(pathway)`;
#'   \item pval -- an enrichment p-value;
#'   \item padj -- a BH-adjusted p-value;
#'   \item ES -- enrichment score, same as in Broad GSEA implementation;
#'   \item NES -- enrichment score normalized to mean enrichment of random samples of the same size;
#'   \item nMoreExtreme` -- a number of times a random gene set had a more
#'      extreme enrichment score value;
#'   \item size -- size of the pathway after removing genes not present in `names(stats)`.
#'   \item leadingEdge -- vector with indexes of leading edge genes that drive the enrichment, see \url{http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading}.
#' }
#' @export
#' 
#' @examples 
#' seuratPbmc <- RunMCA(seuratPbmc, nmcs = 5)
#' ranking <- GetCellGeneRanking(seuratPbmc, reduction = "mca", dims = 1:5)
#' fgseaCelliD(pathways = Hallmark, stats = ranking[[1]])
fgseaCelliD <-
    function(pathways, stats, nperm = 1000, minSize = 10, maxSize = 500, gseaParam = 0) {
        # remove note
        . <- NULL
        ties <- sum(duplicated(stats[stats != 0]))
        if (any(duplicated(names(stats)))) {
            warning("There are duplicate gene names, fgsea may produce unexpected results")
        }
        granularity <- 1000
        seeds <- 1
        minSize <- max(minSize, 1)
        stats <- sort(stats, decreasing = TRUE)
        stats <- abs(stats)^gseaParam
        pathwaysFiltered <- lapply(pathways, function(p) {
            as.vector(na.omit(fastmatch::fmatch(p, names(stats))))
        })
        pathwaysSizes <-
            vapply(X = pathwaysFiltered, FUN = length, FUN.VALUE = as.integer(0))
        toKeep <-
            which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize)
        m <- length(toKeep)
        if (m == 0) {
            return(
                data.table(
                    pathway = character(),
                    pval = numeric(),
                    padj = numeric(),
                    ES = numeric(),
                    NES = numeric()
                )
            )
        }
        pathwaysFiltered <- pathwaysFiltered[toKeep]
        pathwaysSizes <- pathwaysSizes[toKeep]
        K <- max(pathwaysSizes)
        gseaStatRes <-
            vapply(
                X = pathwaysFiltered,
                FUN = calcGseaStat,
                FUN.VALUE = 0,
                stats = stats,
                scoreType = "std"
            )
        pathwayScores <- gseaStatRes
        universe <- seq_along(stats)
        leEs <- rep(0, m)
        geEs <- rep(0, m)
        leZero <- rep(0, m)
        geZero <- rep(0, m)
        leZeroSum <- rep(0, m)
        geZeroSum <- rep(0, m)
        if (m == 1) {
            for (i in seq_len(nperm)) {
                randSample <- sample.int(length(universe), K)
                randEsP <- calcGseaStat(
                    stats = stats,
                    selectedStats = randSample,
                    gseaParam = gseaParam,
                    scoreType = "std"
                )
                leEs <- leEs + (randEsP <= pathwayScores)
                geEs <- geEs + (randEsP >= pathwayScores)
                leZero <- leZero + (randEsP <= 0)
                geZero <- geZero + (randEsP >= 0)
                leZeroSum <- leZeroSum + pmin(randEsP, 0)
                geZeroSum <- geZeroSum + pmax(randEsP, 0)
            }
        } else {
            aux <- calcGseaStatCumulativeBatch(
                stats = stats,
                gseaParam = gseaParam,
                pathwayScores = pathwayScores,
                pathwaysSizes = pathwaysSizes,
                iterations = nperm,
                scoreType = "std",
                seed = 1
            )
            leEs <- get("leEs", aux)
            geEs <- get("geEs", aux)
            leZero <- get("leZero", aux)
            geZero <- get("geZero", aux)
            leZeroSum <- get("leZeroSum", aux)
            geZeroSum <- get("geZeroSum", aux)
        }
        counts <- data.table(
            pathway = seq_len(m),
            leEs = leEs,
            geEs = geEs,
            leZero = leZero,
            geZero = geZero,
            leZeroSum = leZeroSum,
            geZeroSum = geZeroSum
        )

        leEs <- leZero <- geEs <- geZero <- leZeroSum <- geZeroSum <- NULL
        pathway <-
            padj <- pval <- ES <- NES <- geZeroMean <- leZeroMean <- NULL
        nMoreExtreme <- nGeEs <- nLeEs <- size <- NULL
        leadingEdge <- NULL
        pvals <- counts[, list(
            pval = min(
                (1 + sum(leEs)) / (1 + sum(leZero)),
                (1 + sum(geEs)) / (1 + sum(geZero))
            ),
            leZeroMean = sum(leZeroSum) / sum(leZero),
            geZeroMean = sum(geZeroSum) / sum(geZero),
            nLeEs = sum(leEs),
            nGeEs = sum(geEs)
        ), by = .(pathway)]
        pvals[, `:=`(padj, p.adjust(pval, method = "BH"))]
        pvals[, `:=`(ES, pathwayScores[pathway])]
        pvals[, `:=`(NES, ES / ifelse(ES > 0, geZeroMean, abs(leZeroMean)))]
        pvals[, `:=`(leZeroMean, NULL)]
        pvals[, `:=`(geZeroMean, NULL)]
        pvals[, `:=`(nLeEs, NULL)]
        pvals[, `:=`(nGeEs, NULL)]
        pvals[, `:=`(pathway, names(pathwaysFiltered)[pathway])]
        pvals <- pvals[]
        return(pvals)
    }

fgseaCelliDPar <- function(x, pathways, features, nperm, minSize, maxSize) {
    x <- (x + 0.1)^-1
    names(x) <- features
    gseaResults <- CelliD::fgseaCelliD(
        pathways = pathways,
        stats = x,
        nperm = nperm,
        minSize = minSize,
        maxSize = maxSize
    )
}
