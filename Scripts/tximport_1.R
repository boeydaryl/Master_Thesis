dir <- "/Volumes/JetStream/Thesis_Raw_Data/"
samples <- read.table(file.path(dir, "samples_1.txt"), header = TRUE)
files <- file.path(dir, samples$run, "quant.sf")
names(files) <- paste0("sample", 1:3)
all(file.exists(files))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi,colData = samples, design = ~ 1)
dds <- DESeq(ddsTxi)
res <- results(dds)
res