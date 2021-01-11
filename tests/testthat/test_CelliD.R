example_mat <- as.matrix(GetAssayData(seuratPbmc, assay = "RNA", slot = "counts"))
colnames(example_mat) <- paste0("cell", seq(50))
rownames(example_mat) <- paste0("gene", seq(2000))

test_that("MCA cell and gene coordinates match when gene is expressed in only one cell Seurat ver", {
    test_mat <- example_mat
    # Create dummy gene expressed in one cell
    test_mat["gene1",] <- 0
    test_mat["gene1","cell1"] <- 100
    
    # Create Seurat object
    seurat <- CreateSeuratObject(test_mat)
    seurat <- Seurat::NormalizeData(seurat)
    seurat <- RunMCA(seurat, nmcs = 10)
    A <- Seurat::Loadings(seurat, "mca")["gene1",]
    B <- Seurat::Embeddings(seurat, "mca")["cell1",]
    testthat::expect_equal(A, B, tolerance = 0.1)
})


test_that("MCA cell and gene coordinates match when gene is expressed in only one cell SCE ver", {
    test_mat <- example_mat
    # Create dummy gene expressed in one cell
    test_mat["gene1",] <- 0
    test_mat["gene1","cell1"] <- 1
    
    # Create SCE object
    sce <- SingleCellExperiment(assays = list(counts = test_mat))
    sce <- scater::logNormCounts(sce)
    sce <- RunMCA(sce, nmcs = 10)
    A <- reducedDim(sce, "MCA")["cell1",]
    B <- attributes(reducedDim(sce, "MCA"))$genesCoordinates["gene1",]
    testthat::expect_equal(A, B, tolerance = 0.1)
})


# MCA centroid coordinates ------------------------------------------------

test_that("MCA centroids coordinates match when gene is expressed in 5 cells Seurat ver", {
    test_mat <- example_mat
    # Create dummy gene expressed in one cell
    test_mat["gene1",] <- 0
    test_mat["gene1",1:10] <- 10000
    
    # Create Seurat object
    seurat <- CreateSeuratObject(test_mat)
    seurat <- Seurat::NormalizeData(seurat)
    seurat <- RunMCA(seurat, nmcs = 10)
    A <- Seurat::Loadings(seurat, "mca")["gene1",]
    B <- colMeans(Seurat::Embeddings(seurat, "mca")[1:10,])
    testthat::expect_equal(A, B, tolerance = 0.1)
})


test_that("MCA centroids coordinates match when gene is expressed in 5 cells SCE ver", {
    test_mat <- example_mat
    # Create dummy gene expressed in one cell
    test_mat["gene1",] <- 0
    test_mat["gene1",1:10] <- 10000
    
    # Create SCE object
    sce <- SingleCellExperiment(assays = list(counts = test_mat))
    sce <- scater::logNormCounts(sce)
    sce <- RunMCA(sce, nmcs = 10)
    A <- colMeans(reducedDim(sce, "MCA")[1:10,])
    B <- attributes(reducedDim(sce, "MCA"))$genesCoordinates["gene1",]
    testthat::expect_equal(A, B, tolerance = 0.1)
})

# MCA centroid coordinates ------------------------------------------------

test_that("Testing GetCellGeneSet SCE ver", {
    set.seed(85)
    test_mat <- example_mat
    # Create dummy gene expressed in one cell that express only that gene
    test_mat["gene1",] <- 0
    test_mat["gene1","cell1"] <- 10000
    test_mat[-1,-1] <- (test_mat[-1,-1] + 1)
    # Create SCE object
    sce <- SingleCellExperiment(assays = list(counts = test_mat))
    sce <- scater::logNormCounts(sce)
    sce <- RunMCA(sce, nmcs = 10)
    Dist <- GetCellGeneDistance(sce, dims = 1:10)
    testthat::expect_equal(Dist["gene1","cell1"], 0, tolerance = 0.1)
    GS <- GetCellGeneSet(sce, dims = 1:10)
    Top1 <- GS$cell1[[1]]
    testthat::expect_identical(Top1, "gene1")
})

test_that("Testing GetCellGeneSet Seurat ver", {
    set.seed(85)
    test_mat <- example_mat
    # Create dummy gene expressed in one cell that express only that gene
    test_mat["gene1",] <- 0
    test_mat["gene1","cell1"] <- 10000
    test_mat[-1,-1] <- (test_mat[-1,-1] + 1)
    # Create SCE object
    seurat <- CreateSeuratObject(test_mat)
    seurat <- Seurat::NormalizeData(seurat)
    seurat <- RunMCA(seurat, nmcs = 10)
    Dist <- GetCellGeneDistance(seurat, dims = 1:10)
    testthat::expect_equal(Dist["gene1","cell1"], 0, tolerance = 0.1)
    GS <- GetCellGeneSet(seurat, dims = 1:10)
    Top1 <- GS$cell1[[1]]
    testthat::expect_identical(Top1, "gene1")
})


#  CheckCelliDArg ---------------------------------------------------------

test_that("Testing efficacy CheckCelliDArg Seurat ver", {
    seurat <- CreateSeuratObject(example_mat)
    seurat <- Seurat::NormalizeData(seurat)
    # features and cells to choose
    genesA <- sample(rownames(seurat), 1000)
    cells <- sample(colnames(seurat), 25)
    
    # MCA Arg
    seurat <- RunMCA(seurat, nmcs = 10, features = genesA)
    expect_identical(rownames(Loadings(seurat, reduction = "mca")), genesA)
    
    # Cell Gene distance Arg
    genesB <- sample(genesA, 500)
    Dist <- GetCellGeneDistance(seurat, dims = 1:10, features = genesB, cells = cells)
    expect_identical(rownames(Dist), genesB)
    expect_identical(colnames(Dist), cells)
    # Cell GeneSet Arg
    GS <- GetCellGeneSet(seurat, dims = 1:10, features = genesB, cells = cells)
    expect_true(all(unique(unlist(GS)) %in% genesB))
    expect_identical(names(GS), cells)
    genes2 <- 
    path1 <- sample(genesB, 20)
    path2 <- sample(genesB, 20)
    pathways <- list(path1 = path1, path2= path2)
    HGT <- RunCellHGT(seurat, pathways = pathways, dims = 1:10, n.features = 50, features = genesB)
})

test_that("Testing efficacy CheckCelliDArg SCE ver",{
    sce <- SingleCellExperiment(assays = list(counts = example_mat))
    sce <- scater::logNormCounts(sce)
    # features and cells to choose
    genesA <- sample(rownames(sce), 1000)
    cells <- sample(colnames(sce), 25)
    
    # MCA Arg
    sce <- RunMCA(sce, nmcs = 10, features = genesA)
    expect_identical(rownames(attributes(reducedDim(sce, "MCA"))$genesCoordinates), genesA)
    
    # Cell Gene distance Arg
    genesB <- sample(genesA, 500)
    Dist <- GetCellGeneDistance(sce, dims = 1:10, features = genesB, cells = cells)
    expect_identical(rownames(Dist), genesB)
    expect_identical(colnames(Dist), cells)
    # Cell GeneSet Arg
    GS <- GetCellGeneSet(sce, dims = 1:10, features = genesB, cells = cells)
    expect_true(all(unique(unlist(GS)) %in% genesB))
    expect_identical(names(GS), cells)
    genes2 <- 
        path1 <- sample(genesB, 20)
    path2 <- sample(genesB, 20)
    pathways <- list(path1 = path1, path2= path2)
    HGT <- RunCellHGT(sce, pathways = pathways, dims = 1:10, n.features = 50, features = genesB)
})
