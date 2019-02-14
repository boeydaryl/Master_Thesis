#Enhanced Volcano Basic Script
library(EnhancedVolcano)
EnhancedVolcano(geneLevelStats,
                
                lab = geneLevelStats$rn,
                
                x = "log2fc",
                
                y = "padj")
