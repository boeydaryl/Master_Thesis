####Script to integrate DESeq and GSA Workflow########

###### To be Loaded Prior to running #####
library(DESeq2)
library(readr)
library(tximport)
library(GenomicFeatures)
library(biomaRt)
library(EnsDb.Hsapiens.v86)
library(piano)
library(EnhancedVolcano)

#### Directory for abundance files ######
args <- commandArgs(TRUE)
dir1 <- args[1]

#### To load sample metadata and search for related salmon files #####
samples <- read.csv(file.path(dir1, "layout.csv"), header = TRUE)
files <- file.path(dir1, samples$Well.number,"quant.sf")

#### Annotation from transcripts to genes #####
txdb <- makeTxDbFromGFF("/Volumes/scRNAseq_1/Homo_sapiens/ENSEMBL/Homo_sapiens.GRCh38.95.gtf")
k <- keys(txdb, keytype ="TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
names(files) <- paste0(samples$Well.number, "_sample", 1:nrow(samples))
all(file.exists(files))

#### Actual importation step ####
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi,colData = samples, design = ~ Condition)

#### Runs DESeq and saves results as res ######
dds <- DESeq(ddsTxi)
res <- results(dds)
res

write.csv(res, file = "deseq_output.csv")

##### To prepare geneLevelStats file from DESeq Res #####
res <- res[ ! (is.na(res$pvalue) | is.na(res$padj)), ]
geneLevelStats <- as.data.frame(res[,c("log2FoldChange","padj")])
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),mart=mart)
geneLevelStats <- merge(x=geneLevelStats, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
colnames(geneLevelStats) <- c("ensembl","log2fc","padj","gene")
rownames(geneLevelStats) <- geneLevelStats[,1]
geneLevelStats <- geneLevelStats[order(geneLevelStats$padj),]

write.csv(geneLevelStats, file = "geneLevelStats.csv")

#### GSA Run #######
gsc <- loadGSC("/Volumes/scRNAseq_1/Gene_Sets/KEGG/c2.cp.kegg.v6.2.symbols.gmt.txt", type="gmt")
padj <- geneLevelStats$padj
log2fc <- geneLevelStats$log2fc
names(padj) <- names(log2fc) <- geneLevelStats$gene
gsaRes <- runGSA(padj, log2fc, gsc=gsc)

#### PCA of samples ####
rld <- vst(dds)
ramp <- 1:3/3
cols <- c(rgb(ramp, 0, 0),
          rgb(0, ramp, 0),
          rgb(0, 0, ramp),
          rgb(ramp, 0, ramp))
print ( plotPCA( rld, intgroup = c("Condition")))

#### Diagnostic DESeq plots ####
hist (res$padj, breaks=20, col="grey")
plotMA(res, ylim = c(-1, 1))
plotDispEsts(dds)
EnhancedVolcano(geneLevelStats,
                
                lab = geneLevelStats$gene,
                
                x = "log2fc",
                
                y = "padj")

GSAheatmap(gsaRes, cutoff = 10, adjusted = FALSE, ncharLabel = 40,
           cellnote = "pvalue", columnnames = "full", colorkey = TRUE,
           colorgrad = NULL, cex = NULL)

