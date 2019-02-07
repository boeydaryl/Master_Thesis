#Piano GSA Script for KEGG Pathways

###### To be Loaded Prior to running #####
library(biomaRt)
library(piano)


##### To prepare geneLevelStats file from DESeq Res #####
res <- res[ ! (is.na(res$pvalue) | is.na(res$padj)), ]
geneLevelStats <- as.data.frame(res[,c("log2FoldChange","padj")])
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),mart=mart)
geneLevelStats <- merge(x=geneLevelStats, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
colnames(geneLevelStats) <- c("ensembl","log2fc","padj","gene")
rownames(geneLevelStats) <- geneLevelStats[,1]
geneLevelStats <- geneLevelStats[order(geneLevelStats$padj),]

#### Continue from here if geneLevelStats file exists #####
geneLevelStats <- read.csv("/Volumes/scRNAseq_1/P5_vs_IM90/Analysis/p5vsIM90_geneLevelStat.csv")

#### To load appropriate Gene set under gsc #######
#gsc <- loadGSC("/Volumes/scRNAseq_1/Gene_Sets/KEGG/c2.cp.kegg.v6.2.symbols.gmt.txt", type="gmt")
gsc <- loadGSC("/Volumes/scRNAseq_1/Gene_Sets/GSEA/Hallmark_sets/h.all.v6.2.symbols.gmt.txt", type="gmt")
padj <- geneLevelStats$padj
log2fc <- geneLevelStats$log2fc
names(padj) <- names(log2fc) <- geneLevelStats$gene
gsaRes <- runGSA(padj, log2fc, gsc=gsc)


##### for network plot, remember to zoom plot first! #######
networkPlot(gsaRes, "distinct", "both", adjusted=T, ncharLabel=Inf, significance=0.1,
nodeSize=c(3,20), edgeWidth=c(1,5), overlap=10,
scoreColors=c("red", "orange", "yellow", "blue", "lightblue", "lightgreen"))


###### for heatmap #########
GSAheatmap(gsaRes, cutoff = 5, adjusted = FALSE, ncharLabel = 25,
           cellnote = "pvalue", columnnames = "full", colorkey = TRUE,
           colorgrad = NULL, cex = NULL)


##### For GSA Summary Table and Boxplot #####
#GSAsummaryTable(gsaRes, save=T, file="/Volumes/scRNAseq_1/P5_vs_IM90/Analysis/p5_vs_im90gsares_KEGG.txt")
geneSetSummary(gsaRes, "HALLMARK_GLYCOLYSIS")
boxplot(list(-log10(geneLevelStats$padj),
             -log10(geneSetSummary(gsaRes,"HALLMARK_APOPTOSIS")$geneLevelStats)),
        names=c("all","HALLMARK_APOPTOSIS"))