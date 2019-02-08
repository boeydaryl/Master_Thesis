# Script to import SmartSeq2 transcript abundance files through DESeq analysis

###### Necessary libraries to run script #######
library(DESeq2)
#library(readr)
#library(tximport)

#### Directory for abundance files ######
dir1 <- "/Volumes/scRNAseq_1/SS2_15_0150-0151/SmartSeq_Expression/Separated_counts/"

#### To load sample metadata and search for related count files #####
samples <- read.csv("/Volumes/scRNAseq_1/SS2_15_0150-0151/patient_1_layout.csv", header = TRUE)
patient1dir <- "/Volumes/scRNAseq_1/SS2_15_0150-0151/SmartSeq_Expression/patient1.csv"
cts <- as.matrix(read.csv(patient1dir, sep = "\t"))



#### Runs DESeq and saves results as res ######
dds <- DESeq(ddsTxi)
res <- results(dds)
res
