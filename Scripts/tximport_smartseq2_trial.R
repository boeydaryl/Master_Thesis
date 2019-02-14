# Script to import SmartSeq2 transcript abundance files through DESeq analysis

###### Necessary libraries to run script #######
library(DESeq2)
library(data.table)
#library(readr)
#library(tximport)

#### Directory for abundance files ######
dir1 <- "/Volumes/scRNAseq_1/SS2_15_0150-0151/SmartSeq_Expression/Separated_counts/"

#### To load sample metadata and search for related count files #####
colData <- read.csv("/Volumes/scRNAseq_1/SS2_15_0150-0151/patient_1_layout.csv", header = TRUE, row.names = 1)
patient1dir <- "/Volumes/scRNAseq_1/SS2_15_0150-0151/SmartSeq_Expression/patient1.csv"
cts <- data.table(read.csv(patient1dir, sep = "\t"))
unique_cts <- cts[!duplicated(cts$gene)]
write.csv(unique_cts, file = "/Volumes/scRNAseq_1/SS2_15_0150-0151/SmartSeq_Expression/patient1_unique.csv", row.names = FALSE)
cts <- read.csv("/Volumes/scRNAseq_1/SS2_15_0150-0151/SmartSeq_Expression/patient1_unique.csv", row.names = "gene")
colData <- colData[,c("Condition", "Cell.type")]

#### Runs DESeq and saves results as res ######
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = ~ Condition)
dds
dds <- DESeq(dds)
res <- results(dds)
res
write.csv(res, "/Volumes/scRNAseq_1/SS2_15_0150-0151/SmartSeq_Expression/p1_analysis/deseq.csv")
