# LRT for Timecourse FOXA2 OE RNA-Seq DE Analysis
# Fig. 2A

library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
set.seed(7712)
library(clusterProfiler)

## Read in counts and sample table ----------------------------------------------------------------------------------------------------

ff <- list.files(path = "/projects/b1126/Viriya/rna-seq/mr283_mr300", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE)
counts.files <- lapply(ff, read.table, skip = 4)
counts <- as.data.frame(counts.files)
counts <- as.data.frame(sapply(counts.files, function(x) x[ , 2 ]))

f <- str_extract(ff, "/mr[:digit:]*R") %>% str_extract("mr[:digit:]*")
colnames(counts) <- f
row.names(counts) <- counts.files[[1]]$V1
# write.csv(counts, "/projects/p20023/Viriya/analysis/foxa2/rna-seq/raw_counts.csv", row.names = TRUE)

counts <- read.csv("/projects/p20023/Viriya/analysis/foxa2/rna-seq/raw_counts.csv", row.names = 1)
sampleTable <- data.frame(condition = factor(rep(c("plv_2d", "FOXA2_2d", "FOXA2_1w", "FOXA2_2w", "FOXA2_3w", "FOXA2_4w"), each = 3), 
                                             levels = c("plv_2d", "FOXA2_2d", "FOXA2_1w", "FOXA2_2w", "FOXA2_3w", "FOXA2_4w")),
                          replicate = factor(rep(seq(1,3))))

rownames(sampleTable) <- colnames(counts)

sampleTable_pretty <- data.frame(condition = factor(rep(c("Control 2 days", "2 days", "1 week", "2 weeks", "3 weeks", "4 weeks"), each = 3), 
                                             levels = c("Control 2 days", "2 days", "1 week", "2 weeks", "3 weeks", "4 weeks")),
                          replicate = factor(rep(seq(1,3))))

rownames(sampleTable_pretty) <- colnames(counts)

## DESeq2 -----------------------------------------------------------------------------------------------------------------------------

dds_lrt <- DESeqDataSetFromMatrix(counts, sampleTable, design = ~condition)
dds_lrt$condition <- relevel(dds_lrt$condition, ref = "plv_2d")
keep <- rowSums(counts(dds_lrt)) > 1
dds_lrt <- dds_lrt[keep,]
dds_lrt <- DESeq(dds_lrt, test = "LRT", reduced = ~1) # DESeq2 with LTR test compared to just the reduced model
res_lrt <- results(dds_lrt, alpha = 0.0001) # selects for FDR 0.0001

## Heatmap parameters -----------------------------------------------------------------------------------------------------------------
mat_colors <- list(
  replicate = c(brewer.pal(3, "Accent")),
  condition = c(brewer.pal(6, "Set1")))
names(mat_colors$replicate) <- unique(sampleTable$replicate)
names(mat_colors$condition) <- unique(sampleTable$condition)

col_anno <- HeatmapAnnotation(df = sampleTable,
                              which = 'col',
                              col = mat_colors
)
## Cluster functions ------------------------------------------------------------------------------------------------------------------

# get the gene names from each cluster in heatmap draw object. returns a list of vectors
get_cluster_genes <- function(mat, ht, n, res) {
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

# given a matrix of values, will divide the values into a list of df based on gene name list
cluster_values_to_df <- function(full_mat, n, list_of_genes) {
  cl_list_df <- list()
  for (i in 1:n) {
    df <- full_mat[list_of_genes[[i]],]
    cl_list_df[[i]] <- df 
  }
  cl_list_df
}

## Rlog Transformation and Z-score Scaling --------------------------------------------------------------------------------------------

degenes <- res_lrt %>% subset(padj < 0.0001)
ltr_rlog <- rlog(dds_lrt)
rlog_de <- assay(ltr_rlog) %>% subset(rownames(ltr_rlog) %in% rownames(degenes))
rlog_de_scaled <- t(scale(t(rlog_de)))

## Heatmap ----------------------------------------------------------------------------------------------------------------------------

hmap_k6 <- Heatmap(rlog_de_scaled,
                             name = "scaled",
                             
                             # Row Params
                             row_km = 6,
                             show_row_names = FALSE,
                             row_title_rot=0,
                             cluster_row_slices = FALSE,
                             # row_gap = unit(2, "mm"),
                             row_km_repeats = 50,
                             border = TRUE,
                             
                             # Column Params
                             cluster_columns = FALSE,
                             column_title = "Rlog Transformed Expression for all DE genes",
                             top_annotation = col_anno)

pdf("heatmaps/rna_all_genes_k6.pdf")
ht_k6 <- draw(hmap_k6, merge_legend = TRUE)
dev.off()

k6_genes <- get_cluster_genes(rlog_de_scaled, ht_k6, 6)
saveRDS(k6_genes, "k6_genes.RDS")
k6_values <- cluster_values_to_df(assay(ltr_rlog), 6, k6_genes)
names(k6_values) <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6")
saveRDS(ht_k6, "ht_k6.rds")

## GO --------------------------------------------------------------------------------------------------------------------------------

geneUniverse <- rownames(res_lrt) # all genes that were tested for significance
get_go_results <- function(genes_list) {
  result <- list()
  for ( i in 1:length(genes_list)) {
    res <- enrichGO(gene = genes_list[[i]],
                    universe = geneUniverse,
                    keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db,
                    ont = "BP")
    result[[i]] <- res
  }
  result
}

k6_go_results <- get_go_results(k6_genes)

for (i in 1:length(k6_genes)) {
  plt <- dotplot(k6_go_results[[i]], showCategory = 10)
  ggsave(paste0("plots/go_dotplot", str_split(substitute(k6_genes), "_", simplify = TRUE)[1],"_c",i,".jpeg"), plt)
}
