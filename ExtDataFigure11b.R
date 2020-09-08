## Make 4-way plots for 8h timeoint - TKO vs DKO
table = read.table('data/NGS3269_lfc_pval_LPS8h_comparisons.txt',sep="\t",header=T,stringsAsFactors = F)
merged = table[,c('symbol','MLKLkoC8koN4ko.v.MLKLkoC8ko.pval','MLKLkoC8koN4ko.v.MLKLkoC8ko.lfc','MLKLkoC8ko.v.MLKLko.pval','MLKLkoC8ko.v.MLKLko.lfc','MLKLkoC8ko.v.WT8h.lfc','MLKLkoC8ko.v.WT8h.pval')]
## Define Casp8 down genes - Genes that go down in both of these comparisons: MLKL-KO-Casp8-KO-LPS8h vs. MLKL-KO-LPS8h and
## MLKL-ko-Casp8-ko-LPS8h vs WT-LPS8h
casp8_genes_down = merged[(merged$MLKLkoC8ko.v.MLKLko.lfc < 0 & merged$MLKLkoC8ko.v.MLKLko.pval < 0.05),]
casp8_genes_down = casp8_genes_down[(casp8_genes_down$MLKLkoC8ko.v.WT8h.lfc < 0 & casp8_genes_down$MLKLkoC8ko.v.WT8h.pval < 0.05),]
# Total 310 Casp8 dependent downregulated genes
genes_to_color2 = casp8_genes_down
## Define genes that go up in the Triple KO (TKO)
up_in_TKO = merged[(merged$MLKLkoC8koN4ko.v.MLKLkoC8ko.lfc > 0 & merged$MLKLkoC8koN4ko.v.MLKLkoC8ko.pval < 0.05),]
# Define genes that are restored by N4bp1
up_in_TKO_restored_by_N4bp1 = up_in_TKO[up_in_TKO$symbol %in% casp8_genes_down$symbol,]
genes_to_color = up_in_TKO_restored_by_N4bp1

## Color plotting order grey -> blue -> red
merged$order <- ifelse(merged$symbol %in% genes_to_color2$symbol,2,1)
merged[merged$symbol %in% genes_to_color$symbol, 'order'] <- 3
merged = merged[order(merged$order),]


## Create 4 - way plot
red = "#e6550d"
blue = "#068ae5"
## Changing the opacity for color grey
grey = "#d6d4d44d" 

col = rep(grey, nrow(merged))
col[which(merged$symbol %in% genes_to_color2$symbol)] = blue # Genes restored by N4BP1
col[which(merged$symbol %in% genes_to_color$symbol)] = red # N4bp1 independent genes


## Label the following genes
glist = c('Tnf', 'Cxcl1', 'Ccl3', 'Ccl4', 'Csf3','Il6')
df = merged[merged$symbol %in% glist,]


## To add regression line
y = lm(merged$MLKLkoC8ko.v.MLKLko.lfc  ~ merged$MLKLkoC8koN4ko.v.MLKLkoC8ko.lfc, data = merged)
## Add corr coeff
corr = cor.test(merged$MLKLkoC8ko.v.MLKLko.lfc,merged$MLKLkoC8koN4ko.v.MLKLkoC8ko.lfc)

# Plotting
pch = rep(19, nrow(merged))
cex = rep(0.4, nrow(merged))
cex[which(merged$symbol %in% genes_to_color2$symbol)] = 0.6
cex[which(merged$symbol %in% genes_to_color$symbol)] = 0.5
pdf('figures/Ext.data.Figure11b.pdf')
xyplot(merged$MLKLkoC8ko.v.MLKLko.lfc ~ merged$MLKLkoC8koN4ko.v.MLKLkoC8ko.lfc, pch=pch, ylab="MLKLkoC8ko.vs.MLKLko.LPS.8h",
       xlab="MLKLkoC8koN4ko.vs.MLKLkoC8ko.LPS.8h", main = paste0('R = ',round(corr$estimate,2)), xlim = c(-6,6), ylim = c(-6,6),
       par.settings = list(plot.symbol = list(col = col,cex = cex,pch=pch)),
       panel=function(...) { panel.xyplot(...)
         panel.abline(h=1, lty = "dashed", col = "#63636388")
         panel.abline(h=-1, lty = "dashed", col = "#63636388")
         panel.abline(v=1, lty = "dashed", col = "#63636388")
         panel.abline(v=-1, lty = "dashed", col = "#63636388")
         panel.abline(h=0, lty = "solid", col = "#63636388")
         panel.abline(v=0, lty = "solid", col = "#63636388")
         panel.abline(a = 0, b = 1, lty = "solid", col = "#63636388")
         panel.text(x= df$MLKLkoC8koN4ko.v.MLKLkoC8ko.lfc, y = df$MLKLkoC8ko.v.MLKLko.lfc, labels = df$symbol, col = 'black', pos=3)
       })
dev.off()
