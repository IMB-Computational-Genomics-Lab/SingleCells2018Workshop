library(ascend)
em_set <- LoadCellRanger("~/Data/IPSCRetina_scRNA_Aggr_V2/", "GRCh38p7")
em_set1 <- SubsetBatch(em_set, "1")
em_set2 <- SubsetBatch(em_set, "2")

matrix_1 <- GetExpressionMatrix(em_set1, "data.frame")
matrix_2 <- GetExpressionMatrix(em_set2, "data.frame")
cell_1 <- GetCellInfo(em_set1)
cell_2 <- GetCellInfo(em_set2)
gene_1 <- GetGeneInfo(em_set1)
gene_2 <- GetGeneInfo(em_set2)

gene_1 <- gene_1[, -3]
gene_2 <- gene_2[, -3]

library(data.table)
fwrite(matrix_1, "iPSC_RGCscRNASeq_Sample1_RawCounts.tsv", sep = "\t", row.names = TRUE)
fwrite(matrix_2, "iPSC_RGCscRNASeq_Sample2_RawCounts.tsv", sep = "\t", row.names = TRUE)
fwrite(cell_1, "iPSC_RGCscRNASeq_Sample1_CellInfo.tsv", sep = "\t", row.names = FALSE)
fwrite(cell_2, "iPSC_RGCscRNASeq_Sample2_CellInfo.tsv", sep = "\t", row.names = FALSE)
fwrite(gene_1, "iPSC_RGCscRNASeq_Sample1_GeneInfo.tsv", sep = "\t", row.names = FALSE)
fwrite(gene_2, "iPSC_RGCscRNASeq_Sample2_GeneInfo.tsv", sep = "\t", row.names = FALSE)