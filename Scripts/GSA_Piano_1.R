#Piano GSA Script
library(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),mart=mart)
geneLevelStats_p2 <- merge(x=geneLevelStats_p2, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
colnames(geneLevelStats_p2) <- c("ensembl","log2fc","padj","gene")
rownames(geneLevelStats_p2) <- geneLevelStats_p2[,1]
geneLevelStats_p2 <- geneLevelStats_p2[order(geneLevelStats_p2$padj),]
gsc <- loadGSC("/Volumes/scRNAseq_1/Gene_Sets/GSEA/Hallmark_sets/h.all.v6.2.symbols.gmt.txt", type="gmt")
padj <- geneLevelStats_p2$padj
log2fc <- geneLevelStats_p2$log2fc
names(padj) <- names(log2fc) <- geneLevelStats_p2$gene
gsaRes <- runGSA(padj, log2fc, gsc=gsc)
networkPlot(gsaRes, "distinct", "both", adjusted=T, ncharLabel=Inf, significance=0.05,
            nodeSize=c(3,20), edgeWidth=c(1,5), overlap=10,
            scoreColors=c("red", "orange", "yellow", "blue", "lightblue", "lightgreen"))