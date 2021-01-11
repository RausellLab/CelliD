#' Run HyperGeometric Test on cells
#'
#' @param X Seurat or SingleCellExperiment object with mca performed
#' @param pathways geneset to perform hypergeometric test on (named list of genes)
#' @param reduction name of the MCA reduction
#' @param n.features integer of top n features to consider for hypergeometric test
#' @param features vector of features to calculate the gene ranking by default will take everything in the selected mca reduction.
#' @param dims MCA dimensions to use to compute n.features top genes.
#' @param minSize minimum number of overlapping genes in geneset and
#' @param log.trans if TRUE tranform the pvalue matrix with -log10 and convert it to sparse matrix
#' @param p.adjust if TRUE apply Benjamini Hochberg correctionto p-value
#' @importFrom stats phyper
#'
#' @return a matrix of benjamini hochberg adjusted pvalue pvalue or a sparse matrix of (-log10+1) benjamini hochberg adjusted pvalue
#' @export
#'
#' @examples
#' seuratPbmc <- RunMCA(seuratPbmc, nmcs = 5)
#' seuratPbmc <- RunCellHGT(X = seuratPbmc, pathways = Hallmark, dims = 1:5)
RunCellHGT <- function(X, pathways, reduction, n.features, features, dims, minSize, log.trans, p.adjust) {
    UseMethod("RunCellHGT", X)
}


#' @rdname RunCellHGT
#' @export
RunCellHGT.SingleCellExperiment <-
    function(X, pathways, reduction = "MCA", n.features = 200, features = NULL, dims = seq(50), minSize = 10, log.trans = TRUE, p.adjust = TRUE) {
        DT <- GetCellGeneDistance(X, dims = dims, features = features, reduction = reduction)
        message("ranking genes")
        features <- rownames(DT)
        cells <- colnames(DT)

        # Target ------------------------------------------------------------------
        i <- pbapply::pbapply(DT, 2, order)[seq(n.features), ]
        j <- rep(seq(ncol(DT)), each = n.features)
        TargetMatrix <- sparseMatrix(i, j, x = 1, dims = c(length(features), length(cells)), dimnames = list(features, cells))

        # Geneset -----------------------------------------------------------------
        pathways <- lapply(pathways, function(x) x[x %fin% features])
        pathways <- pathways[sapply(pathways, function(x) length(x) >= minSize)]
        message("calculating number of success\n")
        PathwayMat <- pbapply::pbsapply(pathways, function(x) which(features %fin% x), simplify = FALSE)
        PathwayLen <- unlist(lapply(PathwayMat, length))
        j <- rep(seq(length(PathwayMat)), times = PathwayLen)
        PathwayMatrix <- sparseMatrix(unlist(PathwayMat), j, x = 1, dims = c(length(features), length(PathwayMat)), dimnames = list(features, names(PathwayMat)))

        # Hypergeo ----------------------------------------------------------------
        q <- as.data.frame((t(TargetMatrix) %*% PathwayMatrix) - 1)
        m <- sapply(pathways, function(x) sum(x %fin% features))
        n <- sapply(m, function(x) length(features) - x)
        k <- n.features
        message("performing hypergeometric test\n")
        A <- pbapply::pbmapply(
            FUN = function(q, m, n, k) {
                listhyper <- phyper(seq(-1, max(q)), m, n, k, lower.tail = FALSE)[q + 2]
                return(listhyper)
            },
            q = q,
            m = m,
            n = n,
            k = k
        )
        rownames(A) <- rownames(q)
        A <- t(A)
        if (p.adjust) {
            A <- apply(A, 2, function(x) p.adjust(x, "BH"))
        }
        if (log.trans) {
            A <- as.sparse(-log10(A))
        }
        return(A)
    }

#' @rdname RunCellHGT
#' @export
RunCellHGT.Seurat <-
    function(X, pathways, reduction = "mca", n.features = 200, features = NULL, dims = seq(50), minSize = 10, log.trans = TRUE, p.adjust = TRUE) {
        DT <- GetCellGeneDistance(X, dims = dims, features = features, reduction = reduction)
        message("ranking genes")
        features <- rownames(DT)
        cells <- colnames(DT)

        # Target ------------------------------------------------------------------
        i <- pbapply::pbapply(DT, 2, order)[seq(n.features), ]
        j <- rep(seq(ncol(DT)), each = n.features)
        TargetMatrix <- sparseMatrix(i, j, x = 1, dims = c(length(features), length(cells)), dimnames = list(features, cells))

        # Geneset -----------------------------------------------------------------
        pathways <- lapply(pathways, function(x) x[x %fin% features])
        pathways <- pathways[sapply(pathways, function(x) length(x) >= minSize)]
        message("calculating number of success\n")
        PathwayMat <- pbapply::pbsapply(pathways, function(x) which(features %fin% x), simplify = FALSE)
        PathwayLen <- unlist(lapply(PathwayMat, length))
        j <- rep(seq(length(PathwayMat)), times = PathwayLen)
        PathwayMatrix <- sparseMatrix(unlist(PathwayMat), j, x = 1, dims = c(length(features), length(PathwayMat)), dimnames = list(features, names(PathwayMat)))

        # Hypergeo ----------------------------------------------------------------
        q <- as.data.frame((t(TargetMatrix) %*% PathwayMatrix) - 1)
        m <- sapply(pathways, function(x) sum(x %fin% features))
        n <- sapply(m, function(x) length(features) - x)
        k <- n.features
        message("performing hypergeometric test\n")
        A <- pbapply::pbmapply(
            FUN = function(q, m, n, k) {
                listhyper <- phyper(seq(-1, max(q)), m, n, k, lower.tail = FALSE)[q + 2]
                return(listhyper)
            },
            q = q,
            m = m,
            n = n,
            k = k
        )
        rownames(A) <- rownames(q)
        A <- t(A)
        if (p.adjust) {
            A <- apply(A, 2, function(x) p.adjust(x, "BH"))
        }
        if (log.trans) {
            A <- as.sparse(-log10(A))
        }
        return(A)
    }
