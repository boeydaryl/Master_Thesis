#Piano GSA Script
library(biomaRt)
library(piano)
res <- res[ ! (is.na(res$pvalue) | is.na(res$padj)), ]
geneLevelStats <- as.data.frame(res[,c("log2FoldChange","padj")])
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),mart=mart)
"geneLevelStats <- merge(x=geneLevelStats, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
colnames(geneLevelStats) <- c("ensembl","log2fc","padj","gene")
rownames(geneLevelStats) <- geneLevelStats[,1]
geneLevelStats <- geneLevelStats[order(geneLevelStats$padj),]"
#geneLevelStats <- read.csv("/Volumes/scRNAseq_1/P5_vs_IM90/Analysis/p5vsIM90_geneLevelStat.csv")
gsc <- loadGSC("/Volumes/scRNAseq_1/Gene_Sets/GSEA/Hallmark_sets/h.all.v6.2.symbols.gmt.txt", type="gmt")
padj <- geneLevelStats$padj
log2fc <- geneLevelStats$log2fc
names(padj) <- names(log2fc) <- geneLevelStats$gene
gsaRes <- runGSA(padj, log2fc, gsc=gsc)
"networkPlot(gsaRes, "distinct", "both", adjusted=T, ncharLabel=Inf, significance=0.005,
            nodeSize=c(3,20), edgeWidth=c(1,5), overlap=10,
            scoreColors=c("red", "orange", "yellow", "blue", "lightblue", "lightgreen"))"
#GSAsummaryTable(gsaRes, save=T, file="/Volumes/scRNAseq_1/P4_vs_IM90/Analysis/p4_vs_im90gsares.txt")
geneSetSummary(gsaRes, "HALLMARK_APOPTOSIS")
boxplot(list(-log10(geneLevelStats$padj),
             -log10(geneSetSummary(gsaRes,"HALLMARK_APOPTOSIS")$geneLevelStats)),
        names=c("all","HALLMARK_APOPTOSIS"))
geneSetSummary(gsaRes, "HALLMARK_GLYCOLYSIS")
boxplot(list(-log10(geneLevelStats$padj),
             -log10(geneSetSummary(gsaRes,"HALLMARK_GLYCOLYSIS")$geneLevelStats)),
        names=c("all","HALLMARK_GLYCOLYSIS"))
geneSetSummary(gsaRes, "HALLMARK_MTORC1_SIGNALING")
boxplot(list(-log10(geneLevelStats$padj),
             -log10(geneSetSummary(gsaRes,"HALLMARK_MTORC1_SIGNALING")$geneLevelStats)),
        names=c("all","HALLMARK_MTORC1_SIGNALING"))
geneSetSummary(gsaRes, "HALLMARK_INTERFERON_ALPHA_RESPONSE")
boxplot(list(-log10(geneLevelStats$padj),
             -log10(geneSetSummary(gsaRes,"HALLMARK_INTERFERON_ALPHA_RESPONSE")$geneLevelStats)),
        names=c("all","HALLMARK_INTERFERON_ALPHA_RESPONSE"))
geneSetSummary(gsaRes, "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
boxplot(list(-log10(geneLevelStats$padj),
             -log10(geneSetSummary(gsaRes,"HALLMARK_OXIDATIVE_PHOSPHORYLATION")$geneLevelStats)),
        names=c("all","HALLMARK_OXIDATIVE_PHOSPHORYLATION"))