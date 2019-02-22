# Script to import salmon transcript abundance files through DESeq analysis

###### Necessary libraries to run script #######
library(DESeq2)
library(readr)
library(tximport)

#### Directory for abundance files ######
dir1 <- "/Volumes/scRNAseq_1/DESeq runs/p1p2p3vsGSE113957+GSE66053/"

#### To load sample metadata and search for related salmon files #####
samples <- read.csv(file.path(dir1, "p1p2p3vsGSE66053+GSE113957.csv"), header = TRUE)
files <- file.path(dir1, samples$Well.number,"quant.sf")

#### Annotation from transcripts to genes #####
txdb <- makeTxDbFromGFF("/Volumes/scRNAseq_1/Homo_sapiens/ENSEMBL/Homo_sapiens.GRCh38.95.gtf")
k <- keys(txdb, keytype ="TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
names(files) <- paste0(samples$Well.number, "_sample", 1:225)
all(file.exists(files))

#### Actual importation step ####
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi,colData = samples, design = ~ Condition)

#### Runs DESeq and saves results as res ######
dds <- DESeq(ddsTxi)
res <- results(dds)
res

write.csv(res, file = 
            "/Volumes/scRNAseq_1/DESeq runs/p1p2p3vsGSE113957+GSE66053/Analysis/p1p2p3vsGSE113957+GSE66053_deseq.csv")
