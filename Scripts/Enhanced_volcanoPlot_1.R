#Enhanced Volcano Basic Script
library(EnhancedVolcano)
EnhancedVolcano(geneLevelStats_p3,
                
                lab = geneLevelStats_p3$gene,
                
                x = "log2fc",
                
                y = "padj")
