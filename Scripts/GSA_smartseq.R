#Piano GSA Script for downstream analysis of DESeq (SmartSeq2-RSEM)

###### To be Loaded Prior to running #####
library(biomaRt)
library(piano)


##### To prepare geneLevelStats file from DESeq Res #####
res <- res[ ! (is.na(res$pvalue) | is.na(res$padj)), ]
geneLevelStats <- as.data.frame(res[,c("log2FoldChange","padj")])
#setDT(geneLevelStats, keep.rownames = TRUE)[]
#ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),mart=mart)
#geneLevelStats <- merge(x=geneLevelStats, y=ensembl2name, by.x=0, by.y=2)
colnames(geneLevelStats) <- c("log2fc","padj")
#rownames(geneLevelStats) <- geneLevelStats[,4]
geneLevelStats <- geneLevelStats[order(geneLevelStats$padj),]

#### Continue from here if geneLevelStats file exists #####
geneLevelStats <- read.csv("/Volumes/scRNAseq_1/Compiled_Analysis/p1p2p3vsIM90/p0.01/commonGeneLevelStats_exP1.csv")

#### To load appropriate Gene set under gsc #######
gsc <- loadGSC("/Volumes/scRNAseq_1/Gene_Sets/KEGG/c2.cp.kegg.v6.2.symbols.gmt.txt", type="gmt")
padj <- geneLevelStats$padj
log2fc <- geneLevelStats$log2FoldChange
names(padj) <- names(log2fc) <- geneLevelStats$rn
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
#geneSetSummary(gsaRes, "KEGG_RIBOSOME")
GSS_specific <- geneSetSummary(gsaRes, "KEGG_TIGHT_JUNCTION")
write.csv(GSS_specific[["geneLevelStats"]], 
          file = "/Volumes/scRNAseq_1/Compiled_Analysis/p1p2p3vsIM90/p0.01/GSA/KEGG_TightJcn.csv")
boxplot(list(-log10(geneLevelStats$padj),
             -log10(geneSetSummary(gsaRes,"HALLMARK_APOPTOSIS")$geneLevelStats)),
        names=c("all","HALLMARK_APOPTOSIS"))