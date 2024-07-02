library(tidyverse)
library(Signac)
library(Seurat)
library(patchwork)
set.seed(1234)
library(scCustomize)
library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

d14 <- readRDS("../single/simple_d14.RDS")
d2 <- readRDS("../single_d2/simple_d2.RDS")
d21 <- readRDS("../single_d21/simple_d21.RDS")

merge <- merge(d2, y = c(d14, d21), add.cell.ids = c("d2", "d14", "d21"), merge.data = TRUE)
# saveRDS(merge, "simple_merge.RDS")

# rpca --------
# merge <- readRDS("simple_merge.RDS")
merge_list <- SplitObject(merge, split.by = "orig.ident")
merge_list <- lapply(X = merge_list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = merge_list, nfeatures = 3000)
merge_list <- PrepSCTIntegration(object.list = merge_list, anchor.features = features)
merge_list <- lapply(X = merge_list, FUN = RunPCA, features = features)

anchors <- FindIntegrationAnchors(object.list = merge_list, normalization.method = "SCT",
                                  anchor.features = features, reduction = "rpca", k.anchor = 10)

rna_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
rna_integrated <- RunPCA(rna_integrated)
rna_integrated <- RunUMAP(rna_integrated, dims = 1:30)
DimPlot(rna_integrated, group.by = "orig.ident") + ggtitle("RNA")
DefaultAssay(rna_integrated) <- "RNA"
FeaturePlot(rna_integrated, features = c("AR", "KLK3", "TMPRSS2", "ALDH1A1", "FOLH1"))
FeaturePlot(rna_integrated, features = c("FOXA2", "NKX2-1", "POU3F2", "INSM1", "NCAM1", "ENO2"))

# finding clusters -------
DefaultAssay(rna_integrated) <- "integrated"
rna_integrated <- FindNeighbors(rna_integrated, dims = 1:30)
rna_integrated <- FindClusters(rna_integrated, resolution = 0.8)
DimPlot(rna_integrated, group.by = "integrated_snn_res.0.8", label = TRUE, label.box = TRUE) + ggtitle("resolution = 0.8")
rna_integrated <- FindClusters(rna_integrated, resolution = 0.4)
DimPlot(rna_integrated, group.by = "integrated_snn_res.0.4", label = TRUE, label.box = TRUE) + ggtitle("resolution = 0.4")
rna_integrated <- FindClusters(rna_integrated, resolution = 0.2)
DimPlot(rna_integrated, group.by = "integrated_snn_res.0.2", label = TRUE, label.box = TRUE) + ggtitle("resolution = 0.2")
saveRDS("rna_integrated_rpca_default.RDS")

# ARS NES ----
rna_integrated <- readRDS("rna_integrated_rpca_default.RDS")
Idents(rna_integrated) <- rna_integrated$orig.ident
Idents(rna_integrated) <- factor(Idents(rna_integrated), levels = c('d2', 'd14', 'd21'))
rna_integrated$orig.ident <- Idents(rna_integrated)

DefaultAssay(rna_integrated) <- "RNA"
our_AR <- c("AR", "KLK2", "KLK3", "TMPRSS2", "NKX3-1", "FKBP5", "PLPP1", "PMEPA1", "PART1", "ALDH1A3", "STEAP4")
our_NE <- c("SYP", "NCAM1", "ENO2", "INSM1", "SOX2", "NKX2-1")

rna_integrated <- AddModuleScore(rna_integrated, features = list(our_AR), name = "AR_signature", search = TRUE)
rna_integrated <- AddModuleScore(rna_integrated, features = list(our_NE), name = "NE_signature", search = TRUE)

# AR/ARS/NES classification --------

stat_ar <- ifelse(FetchData(rna_integrated, vars = "AR") > 0.6, "AR+", "AR-")
stat_ars <- ifelse(rna_integrated$AR_signature1 >= 0, "ARS+", "ARS-")
stat_nes <- ifelse(rna_integrated$NE_signature1 >= 0, "NES+", "NES-")

pos <- cbind(stat_ar, stat_ars) %>% as.data.frame()
pos$day <- rna_integrated$orig.ident
pos$stat <- paste0(pos$AR, "/", pos$stat_ars)
pos_table <- table(pos$day, pos$stat)

neg <- cbind(stat_ar, stat_nes) %>% as.data.frame()
neg$day <- rna_integrated$orig.ident
neg$stat <- paste0(neg$AR, "/", neg$stat_nes)
neg_table <- table(neg$day, neg$stat)

pos_neg_table <- cbind(pos_table[,c(3,4)], neg_table[,c(1,2)])
prop_sig <- pivot_longer(as_data_frame(prop.table(pos_neg_table, 1)), everything())
prop_sig$day <- factor(rep(c("d2", "d14", "d21"), each = 4), levels = c("d2", "d14", "d21"))
prop_sig$name <- factor(prop_sig$name, levels = c("AR+/ARS+", "AR+/ARS-", "AR-/NES-", "AR-/NES+"))

cell_stat <- data.frame(pos = pos$stat, neg = neg$stat)
rownames(cell_stat) <- rownames(pos)
cell_stat$stat <- ifelse(cell_stat$pos %in% c("AR+/ARS+", "AR+/ARS-", "AR-/NES-", "AR-/NES+"), cell_stat$pos, cell_stat$neg)
rna_integrated$status <- factor(cell_stat$stat, levels = c("AR+/ARS+", "AR+/ARS-", "AR-/NES-", "AR-/NES+"))
rna_integrated$status_day <- factor(paste0(rna_integrated$orig.ident, " ", rna_integrated$status))

# Figure 3C, D, E, F ------

pdf("rna_intergrated_output.pdf")
DimPlot_scCustom(rna_integrated, group.by = "orig.ident", colors_use = ColorBlind_Pal(), pt.size = 0.7) + ggtitle("RNA") &
  labs(x = "rna_UMAP 1", y = "rna_UMAP 2")

ggplot(prop_sig, aes(x = day, y = value, fill = name)) +
  geom_bar(width = 0.5, stat = "identity") +
  labs(x= "Timepoint",y = "Percentage", title) +
  scale_fill_manual(values = c("orange","#0072B2", "#009E73", "#CC79A7", "#F0E442"), name = "Status") +
  theme_classic()

FeaturePlot(rna_integrated, order = T, pt.size = 0.7, features = "AR_signature1") + 
  labs(x = "rna_UMAP 1", y = "rna_UMAP 2", title = "AR Signature")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 1), na.value = "lightgrey")
FeaturePlot(rna_integrated, order = T, pt.size = 0.7, features = "NE_signature1") +
  labs(x = "rna_UMAP 1", y = "rna_UMAP 2", title = "NE Signature")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 0.5), na.value = "lightgrey")

dev.off()

# Figure S3F ----
pdf("rna_single_genes_ordered.pdf")
FeaturePlot(rna_integrated, features = c("AR"), pt.size = 0.7, order = T) & 
  labs(x = "rna_UMAP 1", y = "rna_UMAP 2")
FeaturePlot(rna_integrated, features = c("KLK3"), pt.size = 0.7, order = T) & 
  labs(x = "rna_UMAP 1", y = "rna_UMAP 2")
FeaturePlot(rna_integrated, features = c("NKX3-1"), pt.size = 0.7, order = T) & 
  labs(x = "rna_UMAP 1", y = "rna_UMAP 2")

FeaturePlot(rna_integrated, features = c("SYP"), pt.size = 0.7, order = T) & 
  labs(x = "rna_UMAP 1", y = "rna_UMAP 2")
FeaturePlot(rna_integrated, features = c("NCAM1"), pt.size = 0.7, order = T) & 
  labs(x = "rna_UMAP 1", y = "rna_UMAP 2")
FeaturePlot(rna_integrated, features = c("POU3F2"), pt.size = 0.7, order = T) & 
  labs(x = "rna_UMAP 1", y = "rna_UMAP 2")
dev.off()

