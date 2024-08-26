#01/10/2022

library(unixtools)

set.tempdir("/projects/b1042/YuLab/irina/ATAC-seq/w_bowtie/tmp")
tempdir()
dir.exists(tempdir())
testing <- c(1,2,3)
save(testing, file = "data/testing.RData")


library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)
library(RColorBrewer)
library(cluster)
library(Rsubread)
library(ggfortify)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(ChIPseeker)
library(dplyr)
set.seed(7712)
source("/projects/p20023/Viriya/analysis/foxa2/atac_seq/FOXA2_OE_timecourse/bedops_merge/get_cluster_info.R")
source("/projects/p20023/Viriya/analysis/foxa2/atac_seq/FOXA2_OE_timecourse/bedops_merge/make_func_plots.R")

# FeatureCounts () ----------
#merged <- read_tsv("/projects/b1042/YuLab/irina/ATAC-seq/w_bowtie/both_merged/FOXA2_Cejas_summmit_extend_merge.bed", col_names = FALSE)
#colnames(merged) <- c("chr", "start", "end", "peak_name", "score")
# unique(merged$peak_name) %>% length() # not all the peak names from bedops are unique
#rownames(merged) <- paste0("peaks", rownames(merged))

#saf <- data.frame(GeneID = rownames(merged),
#                  Chr = merged$chr,
#                  Start = merged$start,
#                  End = merged$end,
#                  Strand = rep(".", nrow(merged)))

#f1 <- list.files("/projects/b1042/YuLab/irina/ATAC-seq/w_bowtie/Cejas_bams/unsorted", pattern = ".bam$", full.names = TRUE)
#f2 <- list.files("/projects/b1126/Viriya/ENCODE_atac-seq/m863-m872/bams", pattern = ".bam$", full.names = TRUE)
#f <- c(f1,f2)
#feat_counts <- featureCounts(f, isGTFAnnotationFile = FALSE, annot.ext = saf, isPairedEnd = TRUE)
#df_feat_counts <- feat_counts$counts
#write_csv(as.data.frame(df_feat_counts), "df_feat_counts.csv")

#for_colnames <- c()
#for_colnames[1:28] <- colnames(df_feat_counts[1:28]) %>% substr(1,11)
#for_colnames[29:39] <- colnames(df_feat_counts[29:39]) %>% substr(1,4)

#colnames(df_feat_counts) <- for_colnames
#save(df_feat_counts, file = "data/df_feat_counts.RData")
#write_csv(as.data.frame(df_feat_counts), "featurecounts.csv")

## read in ------------

testing <- c(1,2,3)

df_feat_counts <- read.csv("featurecounts_1.csv", )
rownames(df_feat_counts) <- df_feat_counts[,1]
df_feat_counts <- df_feat_counts[,-1]
samples <- read.csv("sample.csv")
colnames(samples) <- c("SampleID","Condition","Replicate","Group","Batch","bamReads","peaks","PeakCaller")
sampleTable <- samples[, c(2,3)]
rownames(sampleTable) <- samples[, "SampleID"]
sampleTable <- arrange(sampleTable, Condition, Replicate)
colnames(sampleTable) <- tolower(colnames(sampleTable))
counts <- as.data.frame(df_feat_counts)
counts <- relocate(counts, rownames(sampleTable))



## DESeq2 ----------------------------------------------------------------------------------------------------------------------------
#dds <- DESeqDataSetFromMatrix(counts, colData = sampleTable, design = ~condition)


#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep, ]
#dds <- DESeq(dds, test = "LRT", reduced = ~1)
#save(dds, file = "dds_1.RData")
load(file="dds_1.RData")
res <- results(dds, alpha = 0.05)
save(res, file = "res_1.RData")

rlog_res <- rlog(dds, blind = FALSE)
save(rlog_res, file = "rlog_res_1.RData")

#load("data/rlog_res.RData")
depeaks <- res %>% subset(padj < 0.05)
save(depeaks, file = "depeaks_1.RData")

rlog_de_all <- assay(rlog_res) %>% subset(rownames(assay(rlog_res)) %in% rownames(depeaks))
load("rlog_de_all_1.RData")
rlog_de_all <- removeBatchEffect(assay(rlog_res))
save(rlog_de_all, file = "rlog_de_all_1.RData")

pca_rlog <- t(rlog_de_all) %>% prcomp()
PCi<-data.frame(pca_rlog$x,group=sampleTable$group)

pick.col <- c("#666666", "#FF66CC", "#CC9900","#669900","#000000",
	"#CC6666", "#9E0142", "#CC3333", "#00CC33", "#FFCC00","#FF9900")

colors <- factor(PCi$group, labels = pick.col,
			levels = c("ef1_cfce", "LuCaP_CRPC", "LuCaP_NEPC", "MCC", "NCI-H660",
					"D0", "D2", "D7", "D14", "D21", "D28"))
PCi$select <- ""
# Let's just label these items.
ix_label <- c(1,3,5,6,12,22,28,30,31,33,35) 
PCi$select[ix_label] <- as.character(PCi$group[ix_label])

library("RColorBrewer")


plt_pca <- ggplot(PCi, aes(x=PC1, y=PC2, color = group)) + 
		geom_point(aes(shape=group, color=group), size=3) +
  		scale_shape_manual(values=c(16,16,16,16,16, 
                              		17,17,17,17,17,17))+
  		#theme_minimal() +
  		scale_color_manual(values = pick.col) + 
  		labs(title = "ATAC-Seq Time Course", subtitle = paste(num_rows, " peaks, all DE peaks", sep = ""), 
					x="PC1 (28.7%)", y="PC2 (14.78%)") +
		geom_text_repel(aes(label=select), point.size = 7,
                		max.overlaps = 8, show.legend = FALSE,
                		direction="both")
plt_pca
plotPCA(rlog_res)
ggsave("output/pca_all_de_peaks.jpeg", plt_pca)

plt_pca <- autoplot(pca_rlog, data = sampleTable, colour = "condition", label = FALSE,  size = 3) +
  labs(title = "ATAC-Seq Time Course", subtitle = "All DE peaks, Rlog")
plt_pca
plotPCA(rlog_res)
ggsave("output/pca_all_de_peaks.jpeg", plt_pca)

plt_pca1 <- autoplot(pca_rlog, data = sampleTable, colour = "condition", label = FALSE,  size = 3) +
  labs(title = "ATAC-Seq Time Course", subtitle = "All DE peaks, Rlog")
pdf("output/all_de_peaks/pca_all_de_peaks.pdf")
plt_pca1
dev.off()