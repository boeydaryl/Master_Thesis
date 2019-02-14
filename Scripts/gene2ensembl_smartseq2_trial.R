# Script for converting gene symbols in combined RPKM file to Ensembl ID

library(biomaRt)
library(data.table)

rpkms_dir <- "/Volumes/scRNAseq_1/SS2_15_0150-0151/SmartSeq_Expression/rpkms.tab"
rpkms <- data.table(read.csv(rpkms_dir, sep = "\t"))
unique_rpkms <- rpkms[!duplicated(rpkms$gene)]
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
gene2ensembl <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),mart=mart)
labelled_rpkms <- merge(x=unique_rpkms, y=gene2ensembl, by.x="gene", by.y="external_gene_name", all.x=TRUE)

write.csv(labelled_rpkms, file = "/Volumes/scRNAseq_1/SS2_15_0150-0151/SmartSeq_Expression/labelled_RPKMS.tab", sep = "\t")