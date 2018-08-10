#library(ascend)
devtools::load_all("~/CodeRepositories/ascend-dev/")
#devtools::load_all("~/CodeRepositories/ascend")

# BiocParallel configuration
library(BiocParallel)
ncores <- 3
register(MulticoreParam(workers = ncores, progressbar=TRUE), default = TRUE)

# Set path to data
em_set <- loadCellRanger("data/")
training_data <- readRDS(system.file("exdata", "human_cycle_markers.rds", 
                                     package = "scran"))
ensembl_set <- convertGeneID(em_set, new.annotation = "ensembl_gene_id")
em_set <- scranCellCycle(ensembl_set, training_set = training_data)
em_set <- convertGeneID(em_set, new.annotation = "gene_id")
em_set <- normaliseBatches(em_set)
filtered_set <- filterByOutliers(em_set, 
                                 control.threshold = 3,
                                 cell.threshold = 3)
filtered_set <- filterByControl(filtered_set, control = "Mt", pct.threshold = 20)
filtered_set <- filterByControl(filtered_set, control = "Rb", pct.threshold = 50)
filtered_set <- scranNormalise(filtered_set, quickCluster = FALSE, min.mean = 1e-5)
norm_set <- runPCA(filtered_set)
clustered_set <- runCORE(norm_set, conservative = FALSE, remove.outlier = TRUE)
clustered_set <- runTSNE(clustered_set, PCA = TRUE, seed = 1)
tsne_plot <- plotTSNE(clustered_set, Dim1 = 1, Dim2 = 2, group = "cluster")
cycle_plot <- plotTSNE(clustered_set, Dim1 = 1, Dim2 = 2, group = "phase")
batch_plot <- plotTSNE(clustered_set, Dim1 = 1, Dim2 = 2, group = "batch")

deseq1 <- runDESeq(clustered_set, group = "cluster", condition.a = 1, condition.b = c(2, 3), ngenes = 1500)
deseq2 <- runDESeq(clustered_set, group = "cluster", condition.a = 2, condition.b = c(1, 3), ngenes = 1500)
deseq3 <- runDESeq(clustered_set, group = "cluster", condition.a = 3, condition.b = c(1, 2), ngenes = 1500)