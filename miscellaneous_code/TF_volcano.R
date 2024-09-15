library(EnhancedVolcano)
library("plotly")

setwd("D:\\xiaodong\\FOXA2\\TCseq\\by_merged_narrow_peak\\NKX2\\GSE74685_FHCRC2016")
res <- read.table("TF_GSE74685_patient_all_columns.txt", header=TRUE)

selected_genes = list('NKX2-1', 'FOXA2', 'ASCL1', 'AR', 'FOXA1')
pdf(file = "TF_GSE74685_patient_volcano.pdf")
EnhancedVolcano(res, lab = res$GENE_NAME, x = "log2FoldChange", y = "padj", 
title = 'GSE74685_patient CRPC vs NEPC', ylab = bquote(~-Log[10] ~ italic(adjP)),
    pCutoff = 0.05,
    FCcutoff = 1,
    legendLabels=c('NC','Log2FC','adjP', 'adjP & Log2FC'),
    #shape = c(1, 4, 23, 25),    
    col = c('black', 'green', 'purple','red3'),
    colAlpha = 4/5,
    labSize = 3.0,
    selectLab = selected_genes,
    labCol = 'black', labFace = 'bold',
    pointSize = c(ifelse(res$GENE_NAME %in% selected_genes, 5, 1)),
    drawConnectors = TRUE, widthConnectors = 1, colConnectors = 'black') + coord_flip()
dev.off()    
    
    
    




library(EnhancedVolcano)
library("plotly")

setwd("D:\\xiaodong\\FOXA2\\TCseq\\by_merged_narrow_peak\\NKX2\\GSE126078_patient")
res <- read.table("TF_GSE126078_patient_all_columns.txt", header=TRUE)

selected_genes = list('NKX2-1', 'FOXA2', 'ASCL1', 'AR', 'FOXA1')
pdf(file = "TF_GSE126078_patient_volcano.pdf")
EnhancedVolcano(res, lab = res$GENE_NAME, x = "log2FoldChange", y = "padj", 
title = 'GSE126078_patient CRPC vs NEPC', ylab = bquote(~-Log[10] ~ italic(adjP)),
    pCutoff = 0.05,
    FCcutoff = 1,
    legendLabels=c('NC','Log2FC','adjP', 'adjP & Log2FC'),
    #shape = c(1, 4, 23, 25),    
    col = c('black', 'green', 'purple','red3'),
    colAlpha = 4/5,
    labSize = 3.0,
    selectLab = selected_genes,
    labCol = 'black', labFace = 'bold',
    pointSize = c(ifelse(res$GENE_NAME %in% selected_genes, 5, 1)),
    drawConnectors = TRUE, widthConnectors = 1, colConnectors = 'black') + coord_flip()
dev.off()    
    
    
    



library(EnhancedVolcano)
library("plotly")

setwd("D:\\xiaodong\\FOXA2\\TCseq\\by_merged_narrow_peak\\NKX2\\beltran_rnaseq2016")
res <- read.table("TF_beltran_CRPC_NEPC_all_columns.txt", header=TRUE)

selected_genes = list('NKX2-1', 'FOXA2', 'ASCL1', 'AR', 'FOXA1')
pdf(file = "TF_beltran_CRPC_NEPC_volcano.pdf")
EnhancedVolcano(res, lab = res$GENE_NAME, x = "log2FoldChange", y = "padj", 
title = 'Beltran_patient CRPC vs NEPC', ylab = bquote(~-Log[10] ~ italic(adjP)),
    pCutoff = 0.05,
    FCcutoff = 1,
    legendLabels=c('NC','Log2FC','adjP', 'adjP & Log2FC'),
    #shape = c(1, 4, 23, 25),    
    col = c('black', 'green', 'purple','red3'),
    colAlpha = 4/5,
    labSize = 3.0,
    selectLab = selected_genes,
    labCol = 'black', labFace = 'bold',
    pointSize = c(ifelse(res$GENE_NAME %in% selected_genes, 5, 1)),
    drawConnectors = TRUE, widthConnectors = 1, colConnectors = 'black') + coord_flip()
dev.off()    
    
    
    

