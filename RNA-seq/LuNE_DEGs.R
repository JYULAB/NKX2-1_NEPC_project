# LuNE RNA-Seq
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)
set.seed(7712)

# CCS1477 Treatment DEGs Fig 7F -----------------
# mr399 to mr 404
ff <- list.files(path = "/projects/b1126/Viriya/rna-seq/mr399_mr404/", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE)
lune_counts.files <- lapply(ff, read.table, skip = 4)
lune_counts <- as.data.frame(lune_counts.files)
lune_counts <- as.data.frame(sapply(lune_counts.files, function(x) x[ , 2 ]))

f <- str_extract(ff, "/mr[:digit:]*R") %>% str_extract("mr[:digit:]*")
colnames(lune_counts) <- f
row.names(lune_counts) <- lune_counts.files[[1]]$V1

lune_sampleTable <- data.frame(condition = factor(rep(c("DMSO", "CCS"), each = 3)),
                               replicate = factor(rep(seq(1,3))))
rownames(lune_sampleTable) <- colnames(lune_counts)
dds_lune <- DESeqDataSetFromMatrix(lune_counts, lune_sampleTable, design = ~condition)
dds_lune$condition <- relevel(dds_lune$condition, ref = "DMSO")
keep_lune <- rowSums(counts(dds_lune)) > 10
dds_lune <- dds_lune[keep_lune,]
dds_lune <- DESeq(dds_lune)
res_lune <- results(dds_lune, alpha = 0.05)

sig <- subset(res_lune, padj < 0.05 & abs(log2FoldChange) > 1)
norm_counts_lune <- counts(dds_lune, normalized = T) %>% as.data.frame()
file <- merge(as.data.frame(sig), norm_counts_lune, by = 0)
colnames(file)[1] <- "Genes"
write_csv(file, "DMSO_vs_CCS1477_2fold_padj_05.csv")

rlog_lune_full <- rlog(dds_lune, blind = FALSE)
rlog_lune <- assay(rlog_lune_full) %>% subset(rownames(rlog_lune_full) %in% rownames(sig))
rlog_de_scaled <- t(scale(t(rlog_lune)))

hmap_lune <- Heatmap(rlog_de_scaled,
                        name = "scaled",
                        
                        # Row Params
                        show_row_names = FALSE,
                        row_title_rot= 0,
                        cluster_row_slices = FALSE,
                        # row_gap = unit(2, "mm"),
                        border = TRUE,
                        
                        # Column Params
                        cluster_columns = FALSE,
                        column_title = "DMSO vs CCS1477 2 Fold Padj 0.05, n = 590")

ht_lune <- draw(hmap_lune, merge_legend = TRUE)

# LuNE KO DEGs --------------

ff <- list.files(path = "/projects/b1126/Viriya/rna-seq/mr405_mr430", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE)
ff <- ff[1:24] # mr405 to mr428
counts.files <- lapply(ff, read.table, skip = 4)
counts <- as.data.frame(counts.files)
counts <- as.data.frame(sapply(counts.files, function(x) x[ , 2 ]))

f <- str_extract(ff, "/mr[:digit:]*R") %>% str_extract("mr[:digit:]*")
colnames(counts) <- f
row.names(counts) <- counts.files[[1]]$V1

counts_nk <- counts[,1:12] # NKX2-1 experiments
counts_cbp <- counts[,13:24] # CBP/P300 experiments

# FOXA2 and NXK2-1 regulated genes Fig 7E ---------

sampletable_nk <- data.frame(condition = factor(rep(c("V2", "sgNKX2-1", "sgFOXA2", "sgFOXA2 + NKX2-1"), each = 3), 
                                                levels = c("V2", "sgNKX2-1", "sgFOXA2", "sgFOXA2 + NKX2-1")),
                             replicate = factor(rep(seq(1,3))))
rownames(sampletable_nk) <- colnames(counts_nk)

dds_nk <- DESeqDataSetFromMatrix(counts_nk, sampletable_nk, design = ~condition)
dds_nk$condition <- relevel(dds_nk$condition, ref = "V2")
keep_nk <- rowSums(counts(dds_nk)) > 10
dds_nk <- dds_nk[keep_nk,]
dds_nk <- DESeq(dds_nk)
rld_nk <- rlog(dds_nk, blind = FALSE)

res_list_nk <- list(results(dds_nk, contrast = c("condition", "sgNKX2.1", "V2"), alpha = 0.05),
                    results(dds_nk, contrast = c("condition", "sgFOXA2", "V2"), alpha = 0.05),
                    results(dds_nk, contrast = c("condition", "sgFOXA2...NKX2.1", "V2"), alpha = 0.05))
                      
res_list_nk_fc <- lapply(res_list_nk, function(x) subset(x, padj < 0.05 & abs(log2FoldChange) > 1))
rld_nk_2v0 <- subset(assay(rld_nk), rownames(assay(rld_nk)) %in% rownames(res_list_nk_fc[[2]]))

library(circlize)
col <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

# CBP Group of RNA-Seq -------
sampletable_cbp <- data.frame(condition = factor(rep(c("plko.1", "shP300", "shCBP", "shP300 + CBP"), each = 3), 
                                                 levels = c("plko.1", "shP300", "shCBP", "shP300 + CBP")),
                              replicate = factor(rep(seq(1,3))))
rownames(sampletable_cbp) <- colnames(counts_cbp)

dds_cbp <- DESeqDataSetFromMatrix(counts_cbp, sampletable_cbp, design = ~condition)
dds_cbp$condition <- relevel(dds_cbp$condition, ref = "plko.1")
keep_cbp <- rowSums(counts(dds_cbp)) > 10
dds_cbp <- dds_cbp[keep_cbp,]
dds_cbp <- DESeq(dds_cbp)

res_list_cbp  <- list(results(dds_cbp, contrast = c("condition", "shP300", "plko.1"), alpha = 0.05),
                      results(dds_cbp, contrast = c("condition", "shCBP", "plko.1"), alpha = 0.05),
                      results(dds_cbp, contrast = c("condition", "shP300...CBP", "plko.1"), alpha = 0.05))

# combining log2FoldChanges compared to control for both set of experiments ----
fc_nk <- lapply(res_list_nk, function(x) {
  as.data.frame(x) %>% select(log2FoldChange)
}) %>% as.data.frame()
fc_cbp <- lapply(res_list_cbp, function(x) {
  as.data.frame(x) %>% select(log2FoldChange)
}) %>% as.data.frame()

fc_nk_nk <- fc_nk %>% subset(rownames(fc_nk) %in% rownames(res_list_nk_fc[[3]]))
fc_cbp_nk <- fc_cbp %>% subset(rownames(fc_cbp) %in% rownames(res_list_nk_fc[[3]]))

fc_combined_nk <- cbind(fc_nk_nk, fc_cbp_nk)
colnames(fc_combined_nk) <- c("sgNKX2.1", "sgFOXA2", "sgFOXA1+NKX2.1", "shP300", "shCBP", "shP300+CBP")

hmap_fc_nk <- Heatmap(as.matrix(fc_combined_nk),
                      name = "Log2FC",
                      col = col,
                      
                      #Row Params
                      show_row_names = FALSE,
                      cluster_row_slices = FALSE,
                      row_km_repeats = 50,
                      border = TRUE,
                      
                      # Column Params
                      cluster_columns = FALSE,
                      column_names_rot = 45,
                      column_title = "Log2FC respective to control,
sgNKX2.1+FOXA2 vs V2 DEGs, n = 228")
ht_fc_nk <- draw(hmap_fc_nk)

pdf("final_lune_heatmaps.pdf")
ht_lune
ht_fc_nk
dev.off()

