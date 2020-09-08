## Values for Pie charts in Ext. Data Fig 11c
## and up in TKO vs DKO at 8 hrs
table = read.table('data/NGS3269_lfc_pval_LPS8h_comparisons.txt',sep="\t",header=T,stringsAsFactors = F)
merged = table[,c('symbol','MLKLkoC8koN4ko.v.MLKLkoC8ko.pval','MLKLkoC8koN4ko.v.MLKLkoC8ko.lfc','MLKLkoC8ko.v.MLKLko.pval','MLKLkoC8ko.v.MLKLko.lfc','MLKLkoC8ko.v.WT8h.lfc','MLKLkoC8ko.v.WT8h.pval','WT8h.v.untreat.lfc','WT8h.v.untreat.pval')]
## Define Casp8 down genes - Genes that go down in both of these comparisons: MLKL-KO-Casp8-KO-LPS8h vs. MLKL-KO-LPS8h and
## MLKL-ko-Casp8-ko-LPS8h vs WT-LPS8h
casp8_genes_down = merged[(merged$MLKLkoC8ko.v.MLKLko.lfc < 0 & merged$MLKLkoC8ko.v.MLKLko.pval < 0.05),]
casp8_genes_down = casp8_genes_down[(casp8_genes_down$MLKLkoC8ko.v.WT8h.lfc < 0 & casp8_genes_down$MLKLkoC8ko.v.WT8h.pval < 0.05),]
# Total 310 Casp8 dependent downregulated genes
## Genes up in the triple KO
up_in_TKO = merged[(merged$MLKLkoC8koN4ko.v.MLKLkoC8ko.lfc > 0 & merged$MLKLkoC8koN4ko.v.MLKLkoC8ko.pval < 0.05),]

# Casp8 dependent genes are resucned by N4bp1 = 105
overlap = casp8_genes_down[casp8_genes_down$symbol %in% up_in_TKO$symbol,]

# Casp8 depdenet genes are NOT rescued by N4bp1 = 205
casp8_down_up_in_TKO = casp8_genes_down[casp8_genes_down$symbol %in% up_in_TKO$symbol,]
casp8_down_NOTup_in_TKO = casp8_genes_down[!(casp8_genes_down$symbol %in% casp8_down_up_in_TKO$symbol),]

# Total LPS incued genes = 2590
lps_induced_genes = merged[(merged$WT8h.v.untreat.lfc > 0 & merged$WT8h.v.untreat.pval < 0.05),]

# 97/105 Casp8 depen genes that are rescued by N4bp1 are LPS induced
casp8_down_up_in_TKO_lps = casp8_down_up_in_TKO[casp8_down_up_in_TKO$symbol %in% lps_induced_genes$symbol,]

# 117/205 Casp8 depen genes that are NOT rescued by N4bp1 are LPS induced
casp8_down_NOTup_in_TKO_lps = casp8_down_NOTup_in_TKO[casp8_down_NOTup_in_TKO$symbol %in% lps_induced_genes$symbol,]

## Total LPS-induced = 97+117 = 214, N4bp1-dependent = 97/214 = 45%
## Total LPS-noninduced = (205-117) + (105-97) = 96, N4bp1-dependent = 8/96 = 8%


