# Finding differentially accessible peaks
# This was ran on R/3.6.3
# Fig 3A and B

library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)
library(RColorBrewer)
library(Rsubread)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)

set.seed(7712)

source("../../bedops_merge/make_func_plots.R")

# FeatureCounts (copied over, not run. will use featurecounts.csv) ----------
merged <- read_tsv("FOXA2_summmit_extend_merge.bed", col_names = FALSE)
colnames(merged) <- c("chr", "start", "end", "peak_name", "score")
# unique(merged$peak_name) %>% length() # not all the peak names from bedops are unique
rownames(merged) <- paste0("peaks", rownames(merged))

saf <- data.frame(GeneID = rownames(merged),
                  Chr = merged$chr,
                  Start = merged$start,
                  End = merged$end,
                  Strand = rep(".", nrow(merged)))

f <- list.files("/projects/b1126/Viriya/ENCODE_atac-seq/m863-m872/bams", pattern = ".bam$", full.names = TRUE)
feat_counts <- featureCounts(f, isGTFAnnotationFile = FALSE, annot.ext = saf, isPairedEnd = TRUE)
df_feat_counts <- feat_counts$counts
colnames(df_feat_counts) <- colnames(df_feat_counts) %>% substr(1,4)
write_csv(as.data.frame(df_feat_counts), "featurecounts.csv")

## read in ------------

samples <- read.csv("sample.csv")
sampleTable <- samples[, c(2,3)]
rownames(sampleTable) <- samples[, "SampleID"]
sampleTable <- arrange(sampleTable, Condition, Replicate)
colnames(sampleTable) <- tolower(colnames(sampleTable))

counts <- as.data.frame(df_feat_counts)
counts <- relocate(counts, rownames(sampleTable))

## DESeq2 ----------------------------------------------------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(counts, colData = sampleTable, design = ~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds <- DESeq(dds, test = "LRT", reduced = ~1)
res <- results(dds, alpha = 0.0001)

rlog_res <- rlog(dds, blind = FALSE)
depeaks <- res %>% subset(padj < 0.0001)

rlog_de_all <- assay(rlog_res) %>% subset(rownames(assay(rlog_res)) %in% rownames(depeaks))
rlog_de_scaled <- t(scale(t(rlog_de_all)))
## Heatmap parameters -----------------------------------------------------------------------------------------------------------------

mat_colors <- list(
  condition = c((brewer.pal(6, "Set1"))),
  replicate = c("#7FC97F", "#BEAED4"))
names(mat_colors$condition) <- unique(sampleTable$condition)
names(mat_colors$replicate) <- unique(sampleTable$replicate)

col_anno <- HeatmapAnnotation(df = sampleTable,
                              which = 'col',
                              col = mat_colors)

## All DE Peaks -------

hmap_6 <- Heatmap(rlog_de_scaled,
                  name = "scaled",
                  
                  # Row Params
                  row_km = 6,
                  row_km_repeats = 3,
                  show_row_names = FALSE,
                  row_title_rot=0,
                  cluster_row_slices = FALSE,
                  
                  border = TRUE,
                  show_row_dend = FALSE,
                  
                  # Column Params
                  cluster_columns = FALSE,
                  column_title = paste0("Differentially Variable Peaks \npadj 0.0001, number of peaks: ", nrow(rlog_de_scaled)),
                  top_annotation = col_anno)

ht_6 <- draw(hmap_6, merge_legend = TRUE)
ht_6


# functions -------
get_cluster_peak_names <- function(mat, ht, res) {
  n <- length(row_order(ht))
  clust <- list()
  if (is.matrix(mat)) {
    for(i in 1:n) {
      df <- rownames(mat[row_order(ht)[[i]],])
      clust[[i]] <- df
    }
  } else {
    for (i in 1:n) {
      df <- res@rownames[mat[row_order(ht)[[i]]]]
      clust[[i]] <- df
    }
  }
  
  clust
}

peak_coordinates_to_df <- function(full_mat, list_of_peaks) {
  n <- length(list_of_peaks)
  cl_list_df <- list()
  if (length(full_mat) == 5 ) {
    for (i in 1:n) {
      df <- full_mat[list_of_peaks[[i]],]
      df[[4]] <- "+"
      cl_list_df[[i]] <- df
    }
  }
  if (length(full_mat) == 3) {
    for (i in 1:n) {
      df <- full_mat[list_of_peaks[[i]],]
      df[[4]] <- "+"
      cl_list_df[[i]] <- df
    }
  }
  cl_list_df
}

# make functional analysis plots and annotation files

functional_plots <- function(peaks, tagMatrix, path) {
  plotAvgProf(tagMatrix, xlim=c(-3000, 3000))
  ggsave(paste0(path, "avg_profiles.jpg"), last_plot(), width = 5, height = 5)
  
  peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,
                         tssRegion=c(-3000, 3000), verbose=FALSE)
  plotAnnoBar(peakAnnoList)
  ggsave(paste0(path, "feat_dist.jpg"), last_plot(), width = 5, height = 5)
  
  plotDistToTSS(peakAnnoList)
  ggsave(paste0(path, "TSS_dist.jpg"), last_plot(), width = 5, height = 5)
  
  genes <- lapply(peakAnnoList, function(x) as.data.frame(x)$geneId)
  names(genes) = grep("[0-9]", names(genes))
  
  annot <- lapply(peakAnnoList, function(x) x@anno)
  names(annot) <- rep(1:length(annot))
  annot_edb <- lapply(genes, function(x) {
    AnnotationDbi::select(EnsDb.Hsapiens.v86,
                          keys = x,
                          columns = c("GENENAME"),
                          keytype = "ENTREZID")
  })
  
  annot_edb <- lapply(annot_edb, function(x) mutate(x, entrezID = as.character(x$ENTREZID)))
  
  clust_anno <- lapply(1:length(annot), function(x)
    annot_save(x, annot, annot_edb, path)
  )
  compGO <- compareCluster(geneCluster = genes,
                           fun = "enrichGO",
                           ont = "BP",
                           OrgDb = org.Hs.eg.db)
  plt_dot <- dotplot(compGO, title = "GO enrichment")
  ggsave(paste0(path, "go_dot_plot.jpg"), plt_dot, width = 10, height = 10)
  
  jpeg(paste0(path, "tag_heatmap.jpg"))
  tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color=NULL)
  dev.off()
}

annot_save <- function(x, annot, annot_edb, path) {
  annot[[x]] <- left_join(as.data.frame(annot[[x]]), as.data.frame(annot_edb[[x]]), by=c("geneId"="entrezID"))
  write_delim(data.frame(annot[[x]]),
              file = paste0(path, "genes_clust_annot_", x, ".tsv"),
              delim = "\t")
}

plot_things <- function(rlog_scaled, ht, merged, path) {
  peak_name_by_clust <- get_cluster_peak_names(rlog_scaled, ht)
  peak_coordinates_by_clust <- peak_coordinates_to_df(merged, peak_name_by_clust)
  if(length(peak_coordinates_by_clust) == 4) {
    names(peak_coordinates_by_clust) <-  c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")
  }
  if(length(peak_coordinates_by_clust) == 6) {
    names(peak_coordinates_by_clust) <-  c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6")
  }
  
  lapply(1:length(peak_coordinates_by_clust), function(x) 
    write.table(data.frame(peak_coordinates_by_clust[[x]]),
                file = paste0(path, "homer/", names(peak_coordinates_by_clust)[x], ".tsv"), 
                sep = "\t", quote = FALSE, col.names = FALSE))
  
  peaks <- lapply(peak_coordinates_by_clust, makeGRangesFromDataFrame)
  peaks_tag <- lapply(peaks, getTagMatrix, windows = promoter)
  anno <- functional_plots(peaks, peaks_tag, path)
}
# -----
peak_name_by_clust_6 <- get_cluster_peak_names(rlog_de_scaled, ht_6)
peak_coordinates_by_clust_6 <- peak_coordinates_to_df(merged, peak_name_by_clust_6)
names(peak_coordinates_by_clust_6) <-  c("c1", "c2", "c3", "c4", "c5", "c6")
lapply(1:length(peak_coordinates_by_clust_6), function(x) 
  write.table(data.frame(peak_coordinates_by_clust_6[[x]]),
              file = paste0("output/", names(peak_coordinates_by_clust_6)[x],".bed"), 
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE))

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
k6_peaks <- lapply(peak_coordinates_by_clust_6, makeGRangesFromDataFrame)
k6_peaks_tag <- lapply(k6_peaks, getTagMatrix, windows = promoter)
k6_anno <- functional_plots(k6_peaks, k6_peaks_tag, "output/")

