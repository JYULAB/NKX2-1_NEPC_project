library(tidyverse)
library(copykat)
library(scCustomize)

sink("output.txt")
color_clones <- c("#FB8072", "#80B1D3")

combined_matrix <- read_delim("../combined_cna.txt")
# testing with small version
#combined_matrix <- combined_matrix[sample(nrow(combined_matrix), 100), c(1:3, sample(ncol(combined_matrix), 100))]
#hcc <- hclust(parallelDist::parDist(t(combined_matrix[,4:ncol(combined_matrix)]),threads =4, method = "euclidean"), method = "ward.D2")
#tree_2 <- cutree(hcc, 2)

tree_2 <- readRDS("tree2.RDS")

my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
chr <- as.numeric(combined_matrix$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)
days <- str_extract(colnames(combined_matrix)[4:ncol(combined_matrix)], "[:alnum:]+(?=-)")
days <- factor(days, levels = c("d2", "d14", "d21"))


days_color <- ColorBlind_Pal()[days]
clones_color_2 <- color_clones[tree_2]

print(table(days))
print(table(days_color))

row_annotations_2 <- rbind(clones_color_2, days_color)
rownames(row_annotations_2) <- c("Clones", "Days")


col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))


jpeg("heatmap_2_clones.jpeg", res = 100, height = 3750, width = 4000)
heatmap.3(t(combined_matrix[,4:ncol(combined_matrix)]),dendrogram="r", 
          distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), 
          hclustfun = function(x) hclust(x, method="ward.D2"),
          ColSideColors=chr1, RowSideColors=row_annotations_2, RowSideColorsSize = 2,
          Colv=NA, Rowv=TRUE,
          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          keysize=1, density.info="none", trace="none",
          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1, 
          symm=F,symkey=F,symbreaks=T,cex=10, cex.main=4, margins=c(10,5),
          labRow = FALSE, labCol = FALSE)
dev.off()

sink()
