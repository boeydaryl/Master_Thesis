library(DESeq2)
library(readr)
library(tximport)
dir1 <- "/Volumes/scRNAseq_1/DESeq runs/Patient_1_analysis/"
samples <- read.csv(file.path(dir1, "patient_1_layout.csv"), header = TRUE)
dir2 <- "/Volumes/scRNAseq_1/SS2_15_0150-0151/Patient_transcipt_quant_salmon/salmon_quant_1/"
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
res_p1 <- res[ ! (is.na(res$pvalue) | is.na(res$padj)), ]
geneLevelStats_p1 <- as.data.frame(res_p1[,c("log2FoldChange","padj")])