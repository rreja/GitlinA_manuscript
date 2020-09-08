# Volcano plot for all timepoints, comparing WT vs N4bp1 KO
library(EnhancedVolcano)

## Untreated/ 0h
lfc_0h = read.table('data/lfc-pval-KO_vs_WT_0h.txt',sep="\t",header = T,stringsAsFactors = F)
# Adjust point size
lfc_0h$ptsize = 0.8
lfc_0h$ptsize <- ifelse(((lfc_0h$lfc.0h >1 |lfc_0h$lfc.0h < -1) & (lfc_0h$p.0h < 0.05)) , 2, 0.8)
pdf('figures/VolcanoPlot_untreated.pdf')
EnhancedVolcano(lfc_0h,lab = lfc_0h$symbol, x = 'lfc.0h', y = 'p.0h', 
                xlim = c(-5, 5), selectLab = c('Il6','Csf3','Cxcl1','Ccl3','Ccl4','Tnf'),
                ylim = c(0,5),
                DrawConnectors = T, 
                widthConnectors = 1.0, 
                colConnectors = 'black',
                gridlines.major = F,
                gridlines.minor = F,
                FCcutoff = 1,
                col = c("grey30", "grey30", "royalblue", "red2"),
                transcriptPointSize = lfc_0h$ptsize,
                colAlpha = 1.0)
dev.off()

##  4h
lfc_4h = read.table('data/lfc-pval-KO_vs_WT_4h.txt',sep="\t",header = T,stringsAsFactors = F)
lfc_4h$ptsize = 0.8
lfc_4h$ptsize <- ifelse(((lfc_4h$lfc.4h >1 |lfc_4h$lfc.4h < -1) & (lfc_4h$p.4h < 0.05)) , 2, 0.8)
pdf('figures/VolcanoPlot_4h.pdf')
EnhancedVolcano(lfc_4h,lab = lfc_4h$symbol, x = 'lfc.4h', y = 'p.4h', 
                xlim = c(-5, 5), selectLab = c('Il6','Csf3','Cxcl1','Ccl3','Ccl4','Tnf'),
                ylim = c(0,5),
                DrawConnectors = T, 
                widthConnectors = 1.0, 
                colConnectors = 'black',
                gridlines.major = F,
                gridlines.minor = F,
                FCcutoff = 1,
                col = c("grey30", "grey30", "royalblue", "red2"),
                transcriptPointSize = lfc_4h$ptsize,
                colAlpha = 1.0)

dev.off()

##  8h
lfc_8h = read.table('data/lfc-pval-KO_vs_WT_8h.txt',sep="\t",header = T,stringsAsFactors = F)
# Adjust point size
lfc_8h$ptsize = 0.8
lfc_8h$ptsize <- ifelse(((lfc_8h$lfc.8h >1 |lfc_8h$lfc.8h < -1) & (lfc_8h$p.8h < 0.05)) , 2, 0.8)
pdf('figures/VolcanoPlot_8h.pdf')
EnhancedVolcano(lfc_8h,lab = lfc_8h$symbol, x = 'lfc.8h', y = 'p.8h', 
                xlim = c(-5, 5), selectLab = c('Il6','Csf3','Cxcl1','Ccl3','Ccl4','Tnf'),
                ylim = c(0,10),
                DrawConnectors = T, 
                widthConnectors = 1.0, 
                colConnectors = 'black',
                gridlines.major = F,
                gridlines.minor = F,
                FCcutoff = 1,
                col = c("grey30", "grey30", "royalblue", "red2"),
                transcriptPointSize = lfc_8h$ptsize,
                colAlpha = 1.0)

dev.off()
##  16h
lfc_16h = read.table('data/lfc-pval-KO_vs_WT_16h.txt',sep="\t",header = T,stringsAsFactors = F)
lfc_16h$ptsize = 0.8
lfc_16h$ptsize <- ifelse(((lfc_16h$lfc.16h >1 |lfc_16h$lfc.16h < -1) & (lfc_16h$p.16h < 0.05)) , 2, 0.8)
pdf('figures/VolcanoPlot_16h.pdf')
EnhancedVolcano(lfc_16h,lab = lfc_16h$symbol, x = 'lfc.16h', y = 'p.16h', 
                xlim = c(-7, 7), selectLab = c('Il6','Csf3','Cxcl1','Ccl3','Ccl4','Tnf'),
                ylim = c(0,10),
                DrawConnectors = T, 
                widthConnectors = 1.0, 
                colConnectors = 'black',
                gridlines.major = F,
                gridlines.minor = F,
                FCcutoff = 1,
                col = c("grey30", "grey30", "royalblue", "red2"),
                transcriptPointSize = lfc_16h$ptsize,
                colAlpha = 1.0)
dev.off()
