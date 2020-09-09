# 4-way plot comparing 4hr time point with untreated samples for WT and N4bp1 KO
library(lattice)

## NFKB gene list from MsigDB
nfkb = read.table('data/NFkB_msigDB_geneset.txt',sep="\t",header=F,stringsAsFactors = F)

## list of mouse cytokines from here - https://reader.elsevier.com/reader/sd/pii/S0167488914001967?token=33286FB86A9CDB8DF2BB6856278AF41D55D99CF9E6646ADB3AB3E4E2E3949AC96ADD2359A8C580D7C13663C39B106EDE
cytokines = c('Il1a','Il1b','Il6','Il8','Il12b','Tnf','Csf3','Ccl1','Ccl2','Ccl3',
              'Ccl4','Ccl5','Ccl6','Ccl7','Ccl8','Ccl9','Ccl11','Ccl12',
              'Ccl19','Ccl20','Ccl21a','Ccl21b','Ccl21c','Ccl22','Ccl24',
              'Ccl25','Ccl26','Ccl27a','Ccl27b','Ccl28','Xcl1','Cxcl1','Cxcl2',
              'Cxcl3','Cxcl5','Cxcl9','Cxcl10','Cxcl11','Cxcl12','Cxcl13','Cxcl14',
              'Cxcl15','Cxcl16','Cxcl17','Cx3cl1')


table = read.table('data/Partek_analysis_NGS2979_glist.txt',sep="\t",header=T,stringsAsFactors = F)
colnames(table) = c('Symbol','ID','Total_reads',
                    'KO2hr_untreat.pval','KO2hr_untreat.adjpval','KO2hr_untreat.lfc',
                    'KO4hr_untreat.pval','KO4hr_untreat.adjpval','KO4hr_untreat.lfc',
                    'KO8hr_untreat.pval','KO8hr_untreat.adjpval','KO8hr_untreat.lfc',
                    'KO16hr_untreat.pval','KO16hr_untreat.adjpval','KO16hr_untreat.lfc',
                    'WT2hr_untreat.pval','WT2hr_untreat.adjpval','WT2hr_untreat.lfc',
                    'WT4hr_untreat.pval','WT4hr_untreat.adjpval','WT4hr_untreat.lfc',
                    'WT8hr_untreat.pval','WT8hr_untreat.adjpval','WT8hr_untreat.lfc',
                    'WT16hr_untreat.pval','WT16hr_untreat.adjpval','WT16hr_untreat.lfc')


## 4-way plot for 4hr time point
merged = table[,c('Symbol','KO4hr_untreat.adjpval','KO4hr_untreat.lfc','WT4hr_untreat.adjpval','WT4hr_untreat.lfc')]
## Let's filter to only keep significant hits
merged = merged[(merged$KO4hr_untreat.adjpval < 0.05 | merged$WT4hr_untreat.adjpval <0.05), ]

# Set plotting order
merged$order <- ifelse(toupper(merged$Symbol) %in% nfkb$V1,2,1)
merged[merged$Symbol %in% cytokines, 'order'] <- 3
merged = merged[order(merged$order),]


## Create 4 - way plot
red = "#e6550d"
blue = "#068ae5"
grey = "#d6d4d44d" 

col = rep(grey, nrow(merged))
# NFKB genes in blue, cytokines in red
col[which(toupper(merged$Symbol) %in% nfkb$V1)] = blue
col[which(merged$Symbol %in% cytokines)] = red

## To add regression line
y = lm(merged$WT4hr_untreat.lfc  ~ merged$KO4hr_untreat.lfc, data = merged)
## Add corr coeff
corr = cor.test(merged$WT4hr_untreat.lfc,merged$KO4hr_untreat.lfc)

# Plotting
pch = rep(19, nrow(merged))
cex = rep(0.4, nrow(merged))
cex[which(merged$Symbol %in% cytokines)] = 0.6
cex[which(toupper(merged$Symbol) %in% nfkb$V1)] = 0.5

## Label the following genes
glist = c('Il6','Csf3','Cxcl1','Ccl3','Ccl4','Tnf')
df = merged[merged$Symbol %in% glist,]

pdf('figures/4wayPlots_4.pdf')
xyplot(merged$WT4hr_untreat.lfc ~ merged$KO4hr_untreat.lfc, pch=pch, ylab="WT-4hr vs Untreated",
       xlab="KO-4hrs vs Untreated", main = paste0('R = ',round(corr$estimate,4)),xlim = c(-10,15), ylim = c(-10,15),
       par.settings = list(plot.symbol = list(col = col,cex = cex,pch=pch)),
       panel=function(...) { panel.xyplot(...)
         panel.abline(h=1, lty = "dashed", col = "#63636388")
         panel.abline(h=-1, lty = "dashed", col = "#63636388")
         panel.abline(v=1, lty = "dashed", col = "#63636388")
         panel.abline(v=-1, lty = "dashed", col = "#63636388")
         panel.abline(h=0, lty = "solid", col = "#63636388")
         panel.abline(v=0, lty = "solid", col = "#63636388")
         panel.abline(a = 0, b = 1, lty = "solid", col = "#63636388")
         panel.text(x= df$KO4hr_untreat.lfc, y = df$WT4hr_untreat.lfc, labels = df$Symbol, col = 'black', pos=3)
       })

dev.off()
