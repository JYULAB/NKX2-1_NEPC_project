####GSVA####

library(GSVA)
library(Biobase)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggrepel)

#Import Data

##Please note: the following steps should first be taken
# 1) Align with STAR
# 2) Combine raw count outputs from LNCaP/FOXA2, cell line data, and patient tumor data and batch correct with Combat-seq

df <- as.data.frame(read_xlsx("all_mr_cl_new_raw_combat.xlsx"))


df <- df %>% distinct(X, .keep_all=TRUE)
df <- df[complete.cases(df),]
rownames(df)<- df[,1]
df <- df[,-1]
df <- as.data.frame(sapply(df,as.numeric), row.names = rownames(df)) 
df <- df[rowSums(df[])>0,]

####Find gene expression activity scores ####

df <- data.matrix(df, rownames.force = TRUE)

#gene sets
AR <- (read.delim("AR_genes_Nelson.txt", header=TRUE, sep="\t", na.strings="NA"))
NE <- (read.delim("NE_genes_Nelson.txt", header=TRUE, sep="\t", na.strings="NA"))
NE <- unlist(NE) #need to make character vector or else gsva doesn't work
AR <- unlist(AR)
geneSets <- list(NE = NE, AR = AR)



gsva_es <- gsva(df, geneSets, method="gsva", kcdf="Poisson" ) #Poisson for integer counts (raw counts), Gaussian for FPKM
write.table(gsva_es, file = "Nelson_GSVA_scores_wnewdata.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE, na = "")


####Scatterplot####
gsva_es <- read.delim("Beltran_GSVA_scores_wnewdata.txt", header=TRUE, sep="\t", na.strings="NA")
scores <- as.data.frame(gsva_es)
scores <- t(scores)


#metadata
metadata <- read.delim("metadata_5.txt", header=TRUE, sep="\t", na.strings="NA")
rownames(metadata) <- metadata$sample
scores <- merge(scores, metadata, by = 0)
rownames(scores) <- scores[,1]
scores <- scores[,-1]

#Remove cell lines that are not new data or NCI-H660
scores <- scores[-c(1:4, 9, 28:30, 34:36),]
# Hide all of the text labels.
scores$select <- ""
# Let's just label these items.
ix_label <- c(1,3, 5, 8, 11, 14, 17, 20, 23:25, 26, 29, 42, 54, 105:106, 117:120) 
scores$select[ix_label] <- scores$cluster2[ix_label]


require('RColorBrewer') #color palette for the scatter plot
pick.col <- c("#6699FF", "#CC9900", "#FF66CC", "#006633",
              "#000000",
              "#9E0142", "#CC3333",  "#00CC33", "#FFCC00", "#FF9900","#999999", "#CC6666")
colII <- colorRampPalette(pick.col)(7)[factor(7, levels = metadata)]

ggplot(scores, aes(x=NE, y=AR, color = cluster)) + 
  geom_point(aes(shape=cluster, color=cluster), size=2) +
  scale_shape_manual(values=c(16,16,16,16,16, 
                              17,17,17,17,17,16,17))+
  theme_minimal() +
  scale_color_manual(values = pick.col) + 
  labs(x="NE Signature", y="AR Signature")+
  ylim(-0.8, 0.7)+
  xlim(-0.75, 0.9)+
  geom_text_repel(aes(label=select), point.size = 7,
                  max.overlaps = Inf, show.legend = FALSE,
                  direction="both")

