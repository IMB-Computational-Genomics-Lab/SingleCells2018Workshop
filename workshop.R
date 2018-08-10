devtools::load_all("~/CodeRepositories/ascend-dev")
#library(ascend)
library(BiocParallel)

ncores <- parallel::detectCores() - 1
register(MulticoreParam(workers = ncores, progressbar=TRUE), default = TRUE)

em_set <- loadCellRanger("data/")
colInfo <- colInfo(em_set)
batch1_cells <- colInfo$cell_barcode[which(colInfo$batch == 1)]
batch2_cells <- colInfo$cell_barcode[which(colInfo$batch == 2)]
batch1_cells <- sample(batch1_cells, length(batch2_cells), replace = FALSE)

cell_list <- c(batch1_cells, batch2_cells)
small_set <- em_set[, cell_list]
saveRDS(small_set, "SmallRawRetinaSet.rds")

qc_set <- plotGeneralQC(small_set)
filtered_set <- filterByOutliers(small_set, cell.threshold = 3, control.threshold = 3)
filtered_set <- filterByControl(filtered_set, control = "Mt", pct.threshold = 20)
filtered_set <- filterByControl(filtered_set, control = "Rb", pct.threshold = 50)
filtered_set <- filterLowAbundanceGenes(filtered_set, pct.threshold = 1)

filtered_log <- progressLog(filtered_set)
filter_qc <- plotGeneralQC(filtered_set)

norm_set <- normaliseByRLE(filtered_set)
norm_set <- excludeControl(norm_set, control = c("Mt", "Rb"))
norm_set <- calculateCV(norm_set)
pca_set <- runPCA(norm_set)

# Cluster
clustered_set <- runCORE(pca_set, remove.outlier = TRUE, conservative = FALSE)
tsne_set <- runTSNE(clustered_set, seed = 1)
  
# Generate heat map for clusters
object <- clustered_set
cluster_analysis <- clusterAnalysis(object)
cluster_list <- 1:cluster_analysis$nClusters
clustering_results <- lapply(cluster_list, function(n){
  # Determine other clusters to compare against
  other_clusters <- cluster_list[which(cluster_list != n)]
  
  # Perform clustering analysis
  de_results <- runDiffExpression(object, group = "cluster", condition.a = as.character(n), condition.b = as.character(other_clusters))
  
  # Retain significant results only and log2foldchange greater than two
  sig_results <- subset(de_results, de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 2)
  sig_results$cluster <- n
  return(sig_results)
})

# Compile all results together
all_de_results <- dplyr::bind_rows(clustering_results)

# Order by largest fold changes
all_de_results <- all_de_results[order(abs(all_de_results$log2FoldChange), decreasing = TRUE), ]
heatmap_genes <- unique(all_de_results$id)[1:30]

# Extract logged matrix
log_counts <- SingleCellExperiment::logcounts(object)
counts <- log_counts[heatmap_genes, ]

# Scale data
scaled_counts <- scale(t(counts))
count_matrix <- as.matrix(scaled_counts)
rownames(count_matrix) <- scaled_counts