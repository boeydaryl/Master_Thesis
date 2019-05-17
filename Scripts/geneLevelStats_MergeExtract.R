#Script to merge and extract p1p2p3vsIM90 genes in common, extracting l2fc and padj from p1
library(plyr)

#### Loading the relevant CSV files #####

common_genes <- read.csv("/Volumes/scRNAseq_1/Compiled_Analysis/p1-5nGSE113957expnGSE66053/subsetnp1p2p3nGSE113957nGSE66053.csv", header = TRUE)
p1_geneLevelStats <- read.csv("/Volumes/scRNAseq_1/DESeq runs/p1p2p3p4p5_vs_GSE113957_exp/Analysis/P1p2p3p4p5vsGSE113957_geneLevelStats.csv", header = TRUE, row.names = "X")

#### Merging from p1_geneLevelStats to common_genes ####

common_geneLevelStats <- match_df(p1_geneLevelStats, common_genes, on = "ensembl")
write.csv(common_geneLevelStats, file = "/Volumes/scRNAseq_1/Compiled_Analysis/p1-5nGSE113957expnGSE66053/subsetnp1p2p3nGSE113957nGSE66053_geneLevelStats.csv")

common_test <- intersect(common_genes$ensembl, p1_geneLevelStats$ensembl)
write.csv(common_test, file = "/Volumes/scRNAseq_1/Compiled_Analysis/p1-5nGSE113957expnGSE66053/subsetnp1p2p3nGSE113957nGSE66053.csv")
csv_file <- common_genes[common_test,]
csv_file2 <- p1_geneLevelStats[common_test,]