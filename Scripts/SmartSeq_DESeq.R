#Deseq on SmartSeq2 Counts
library(DESeq2)
Smartseq_counts <- read.csv("/Volumes/scRNAseq_1/SS2_15_0150-0151/SmartSeq_Expression/counts.tab", header = TRUE, sep = "\t")
patient1layout <- read.csv("/Volumes/scRNAseq_1/SS2_15_0150-0151/patient_1_layout.csv", header = TRUE)
