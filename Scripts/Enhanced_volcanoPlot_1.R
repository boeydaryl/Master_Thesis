#Enhanced Volcano Basic Script
library(EnhancedVolcano)
EnhancedVolcano(geneLevelStats,
                
                lab = geneLevelStats$gene,
                
                x = "log2fc",
                
                y = "padj")
