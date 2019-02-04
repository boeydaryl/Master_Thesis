library(DESeq2)
library(readr)
library(tximport)
library(GenomicFeatures)
dir1 <- "/Volumes/scRNAseq_1/SS2_15_0150-0151/"
samples <- read.csv(file.path(dir1, "patient_2_layout.csv"), header = TRUE)
dir2 <- "/Volumes/scRNAseq_1/SS2_15_0150-0151/salmon_quant_2/"
files <- file.path(dir2, samples$Well.number,"quant.sf")
txdb <- makeTxDbFromGFF("/Volumes/scRNAseq_1/Homo_sapiens/ENSEMBL/Homo_sapiens.GRCh38.95.gtf")
k <- keys(txdb, keytype ="TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
names(files) <- paste0(samples$Well.number, "_sample", 1:94)
all(file.exists(files))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi,colData = samples, design = ~ Condition)
dds <- DESeq(ddsTxi)
res <- results(dds)
res
res_p2 <- res[ ! (is.na(res$pvalue) | is.na(res$padj)), ]
geneLevelStats_p2 <- as.data.frame(res_p2[,c("log2FoldChange","padj")])
