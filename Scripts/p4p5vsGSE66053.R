# Script to import salmon transcript abundance files through DESeq analysis

###### Necessary libraries to run script #######
library(DESeq2)
library(readr)
library(tximport)
library(GenomicFeatures)

#### Directory for abundance files ######
dir1 <- "/Volumes/scRNAseq_1/DESeq runs/p4p5_vs_GSE66053/"

#### To load sample metadata and search for related salmon files #####
samples <- read.csv(file.path(dir1, "p4p5vsGSE66053.csv"), header = TRUE)
files <- file.path(dir1, samples$Well.number,"quant.sf")

#### Annotation from transcripts to genes #####
txdb <- makeTxDbFromGFF("/Volumes/scRNAseq_1/Homo_sapiens/ENSEMBL/Homo_sapiens.GRCh38.95.gtf")
k <- keys(txdb, keytype ="TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
names(files) <- paste0(samples$Well.number, "_sample", 1:171)
all(file.exists(files))

#### Actual importation step ####
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi,colData = samples, design = ~ Condition)

#### Runs DESeq and saves results as res ######
dds <- DESeq(ddsTxi)
res <- results(dds)
res

write.csv(res, file = 
            "/Volumes/scRNAseq_1/DESeq runs/p4p5_vs_GSE66053/Analysis/p4p5vsGSE66053_deseq.csv")
