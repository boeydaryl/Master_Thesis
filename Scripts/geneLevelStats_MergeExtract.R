#Script to merge and extract p1p2p3vsIM90 genes in common, extracting l2fc and padj from p1
library(plyr)

#### Loading the relevant CSV files #####

common_genes <- read.csv("/Volumes/scRNAseq_1/Compiled_Analysis/p4p5nGSE113957nGSE66053/p4p5nGSE113957nGSE66053.csv", header = TRUE)
p1_geneLevelStats <- read.csv("/Volumes/scRNAseq_1/DESeq runs/p4p5_vs_GSE66053/Analysis/p4p5vsGSE66053_geneLevelStats.csv", header = TRUE, row.names = "X")

#### Merging from p1_geneLevelStats to common_genes ####

common_geneLevelStats <- match_df(p1_geneLevelStats, common_genes, on ="ensembl")
write.csv(common_geneLevelStats, file = "/Volumes/scRNAseq_1/Compiled_Analysis/p4p5nGSE113957nGSE66053/geneLevelStats_from_p4p5vsGSE66053.csv")