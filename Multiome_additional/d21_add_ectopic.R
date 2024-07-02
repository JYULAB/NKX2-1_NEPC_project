library(tidyverse)
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(scCustomize)
set.seed(1234)
options(scipen = 999)
polychrome_pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
polychrome_pal <- polychrome_pal[c(1, 3, 6:36)]

counts_d21 <- Read10X_h5("/projects/b1042/YuLab/Viriya/multiome/d21_same/d21/outs/filtered_feature_bc_matrix.h5")
fragpath_d21 <- "/projects/b1042/YuLab/Viriya/multiome/d21_same/d21/outs/atac_fragments.tsv.gz"

edb <- EnsDb.Hsapiens.v86
seqlevelsStyle(edb) <- "UCSC"

annotation <- GetGRangesFromEnsDb(ensdb = edb)

d21 <- CreateSeuratObject(counts = counts_d21$`Gene Expression`, assay = "RNA", project = "d21")
d21[["ATAC"]] <- CreateChromatinAssay(counts = counts_d21$Peaks,
                                      sep = c(":", "-"),
                                      fragments = fragpath_d21,
                                      annotation = annotation)
DefaultAssay(d21) <- "ATAC"

d21 <- NucleosomeSignal(d21)
d21 <- TSSEnrichment(d21)

VlnPlot(
  object = d21,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 2,
  pt.size = 0
)
saveRDS(d21, "d21_raw.RDS")

v1 <- VlnPlot(d21, features = c("nCount_RNA"), pt.size = 0) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 3500) +
  geom_hline(yintercept = 80000)
v2 <- VlnPlot(d21, features = c("nCount_ATAC"), pt.size = 0) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 10000) +
  geom_hline(yintercept = 150000)
v1_1 <- VlnPlot(d21, features = c("nCount_RNA"), pt.size = 0, y.max = 90000) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 3500) +
  geom_hline(yintercept = 80000)
v2_1 <- VlnPlot(d21, features = c("nCount_ATAC"), pt.size = 0, y.max = 150000) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 10000) +
  geom_hline(yintercept = 150000)
v3 <- VlnPlot(d21, features = c("TSS.enrichment"), pt.size = 0) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 2)
v4 <- VlnPlot(d21, features = c("nucleosome_signal"), pt.size = 0) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 1.5)

v1 + v2 + v1_1 + v2_1 + v3 + v4 + plot_layout(nrow = 3, ncol = 2) + plot_annotation("FOXA2_same")

d21 <- subset(
  x = d21,
  subset = nCount_ATAC < 150000 &
    nCount_RNA < 80000 &
    nCount_ATAC > 10000 &
    nCount_RNA > 3500 &
    nucleosome_signal < 1.5 &
    TSS.enrichment > 2
)
peaks_d21 <- CallPeaks(d21, macs2.path = "/projects/b1126/Viriya/conda/macs/bin/macs2")
peaks_d21 <- keepStandardChromosomes(peaks_d21, pruning.mode = "coarse")
peaks_d21 <- subsetByOverlaps(x = peaks_d21, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts_d21 <- FeatureMatrix(
  fragments = Fragments(d21),
  features = peaks_d21,
  cells = colnames(d21)
)

saveRDS(macs2_counts_d21, "macs2_counts.RDS")

d21[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts_d21,
  fragments = fragpath_d21,
  annotation = annotation
)

DefaultAssay(d21) <- "RNA"
d21 <- SCTransform(d21)
d21 <- RunPCA(d21)
d21 <- FindNeighbors(d21, reduction = "pca", dims = 1:50)
d21 <- RunUMAP(d21, reduction = "pca", dims = 1:50, reduction.name = "umap_rna_alone")
# d21 <- FindClusters(d21)
DimPlot_scCustom(d21, reduction = "umap_rna_alone")
FeaturePlot_scCustom(d21, "FOXA2", pt.size = 0.7)
FeaturePlot_scCustom(d21, "FOXA2-polyA", pt.size = 0.7)

DefaultAssay(d21) <- "RNA"
d21 <- NormalizeData(d21)
pdf("same_FOXA2.pdf")
VlnPlot_scCustom(d21, "FOXA2")
FeaturePlot_scCustom(d21, "FOXA2", pt.size = 0.7)
dev.off()

expression <- FetchData(d21, c("FOXA2", "NKX2-1"), slot = "counts")
F_pos_N_pos <- ifelse((expression$FOXA2 > 0 & expression$`NKX2-1` > 0), TRUE, FALSE)
names(F_pos_N_pos) <- Cells(d21)
F_pos_N_neg <- ifelse(expression$FOXA2 > 0 & expression$`NKX2-1` == 0, TRUE, FALSE)
F_neg_N_pos <- ifelse(expression$FOXA2 == 0 & expression$`NKX2-1` > 0, TRUE, FALSE)

sum(expression$`NKX2-1` > 0)
sum(expression$FOXA2 > 0)

# expression[new_cells,]
d21$F_pos_N_pos <- as.factor(F_pos_N_pos) 
d21$F_pos_N_neg <- as.factor(F_pos_N_neg) 
result <-character(length(F_neg_N_pos))

for (i in 1:length(F_pos_N_pos)) {
  if(F_pos_N_pos[i]) {
    result[i] <- "F_pos_N_pos"
  } else if(F_pos_N_neg[i]) {
    result[i] <- "F_pos_N_neg"
  } else {
    result[i] <- "other"
  }
}

d21$FN_status <- factor(result)

# ALRA -------
library(SeuratWrappers)
DefaultAssay(d21) <- "RNA"
d21_alra <- RunALRA(d21)
# d21_alra <- readRDS("d21_alra.RDS")
DefaultAssay(d21_alra) <- "RNA"
FeaturePlot_scCustom(d21_alra, c("AR", "KLK3", "FOXA2", "NKX2-1"), reduction = "umap_rna_alone",
                     pt.size = 0.7, colors_use = viridis_light_high) + plot_annotation("RNA") &
  labs(x = "rna_UMAP 1", y = "rna_UMAP 2")

FeaturePlot_scCustom(d21_alra, c("SYP", "ASCL1", "ENO2", "NCAM1"), reduction = "umap_rna_alone",
                     pt.size = 0.7, colors_use = viridis_light_high) + plot_annotation("RNA") &
  labs(x = "rna_UMAP 1", y = "rna_UMAP 2")


DefaultAssay(d21_alra) <- "alra"
FeaturePlot_scCustom(d21_alra, c("AR", "KLK3", "FOXA2", "NKX2-1"), reduction = "umap_rna_alone",
                     pt.size = 0.7, colors_use = viridis_light_high) + plot_annotation("ALRA") &
  labs(x = "rna_UMAP 1", y = "rna_UMAP 2")

FeaturePlot_scCustom(d21_alra, c("SYP", "ASCL1", "ENO2", "NCAM1"), reduction = "umap_rna_alone",
                     pt.size = 0.7, colors_use = viridis_light_high) + plot_annotation("ALRA") &
  labs(x = "rna_UMAP 1", y = "rna_UMAP 2")

FeaturePlot_scCustom(d21_alra, c("TMPRSS2", "NEAT1", "HOXB13", "VIM"), reduction = "umap_rna_alone",
                     pt.size = 0.7, colors_use = viridis_light_high) + plot_annotation("RNA") &
  labs(x = "rna_UMAP 1", y = "rna_UMAP 2")

VlnPlot_scCustom(d21, features = c("AR", "KLK3", "FOXA2", "NKX2-1"), group.by = "orig.ident", num_columns = 4) + 
  plot_annotation("D21 RNA")
VlnPlot_scCustom(d21_alra, features = c("AR", "KLK3", "FOXA2", "NKX2-1"), group.by = "orig.ident", num_columns = 4) + 
  plot_annotation("D21_alra")


expression_alra <- FetchData(d21_alra, c("SYP", "FOXA2", "ENO2", "NCAM1", "NKX2-1"), assay = "alra")

F_pos_N_pos_alra <- ifelse(expression_alra$FOXA2 > 0 & expression_alra$`NKX2-1` > 0, TRUE, FALSE)
F_pos_N_neg_alra <- ifelse(expression_alra$FOXA2 > 0 & expression_alra$`NKX2-1` == 0, TRUE, FALSE)
sum(expression_alra$`NKX2-1` > 0)
sum(expression_alra$FOXA2 > 0)

d21_alra$F_pos_N_pos <- as.factor(F_pos_N_pos_alra) 
d21_alra$F_pos_N_neg <- as.factor(F_pos_N_neg_alra)

result_alra <-character(length(F_pos_N_pos_alra))

for (i in 1:length(F_pos_N_pos_alra)) {
  if(F_pos_N_pos_alra[i]) {
    result_alra[i] <- "F_pos_N_pos"
  } else if(F_pos_N_neg_alra[i]) {
    result_alra[i] <- "F_pos_N_neg"
  } else {
    result_alra[i] <- "other"
  }
}

d21_alra$FN_status <- factor(result_alra)
saveRDS(d21_alra@meta.data, "d21_alrea_meta_data.RDS")

F_pos_N_pos_alra %>% table()

d21_alra <- AddModuleScore(d21_alra, features = list(our_AR), name = "AR_signature", search = TRUE, assay = "alra")
d21_alra <- AddModuleScore(d21_alra, features = list(our_NE), name = "NE_signature", search = TRUE, assay = "alra")
FeaturePlot_scCustom(d21_alra, c("AR", "KLK3", "FOXA2", "NKX2-1"), reduction = "umap_rna_alone", pt.size = 1)
FeaturePlot_scCustom(d21_alra, "AR_signature1", pt.size = 1)
FeaturePlot_scCustom(d21_alra, "NE_signature1", pt.size = 1)
DimPlot_scCustom(d21, group.by = "F_pos_N_pos", pt.size = 1)
DimPlot_scCustom(d21, group.by = "F_pos_N_neg", pt.size = 1)
Cell_Highlight_Plot(d21_alra, list(F_pos_N_pos = rownames(d21_alra@meta.data %>% filter(F_pos_N_pos == TRUE))), pt.size = 2)
Cell_Highlight_Plot(d21_alra, list(F_pos_N_neg = rownames(d21_alra@meta.data %>% filter(F_pos_N_neg == TRUE))), pt.size = 2)

VlnPlot_scCustom(d21_alra, features = c("AR_signature1", "NE_signature1"), group.by = c("F_pos_N_pos")) + 
  plot_annotation("FOXA2 postive NKX2-1 Positive") + theme(legend.position = "right")
VlnPlot_scCustom(d21_alra, features = c("AR_signature1", "NE_signature1"), group.by = c("F_pos_N_neg")) + 
  plot_annotation("FOXA2 postive NKX2-1 Negative") + theme(legend.position = "right")

VlnPlot_scCustom(d21_alra, features = c("AR_signature1", "NE_signature1"), group.by = c("FN_status"), colors_use = polychrome_pal) + 
  plot_annotation("FOXA2 NKX2-1 Status") + theme(legend.position = "right")

pdf("d21_alra_FN_0524.pdf")
FeaturePlot_scCustom(d21_alra, order = T, pt.size = 0.7, features = "AR_signature1") + 
  labs(x = "atac_UMAP 1", y = "atac_UMAP 2", title = "AR Signature (RNA)")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 1), na.value = "lightgrey")
FeaturePlot_scCustom(d21_alra, order = T, pt.size = 0.7, features = "AR_signature1") + 
  labs(x = "atac_UMAP 1", y = "atac_UMAP 2", title = "AR Signature (RNA)")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 2.5), na.value = "lightgrey")

FeaturePlot(d21_alra, order = T, pt.size = 0.7, features = "NE_signature1") +
  labs(x = "atac_UMAP 1", y = "atac_UMAP 2", title = "NE Signature (RNA)")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 0.5), na.value = "lightgrey")
FeaturePlot(d21_alra, order = T, pt.size = 0.7, features = "NE_signature1") +
  labs(x = "atac_UMAP 1", y = "atac_UMAP 2", title = "NE Signature (RNA)")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 0.8), na.value = "lightgrey")

Cell_Highlight_Plot(d21_alra, list(F_pos_N_pos = rownames(d21_alra@meta.data %>% filter(F_pos_N_pos == TRUE))), pt.size = 0.7)
Cell_Highlight_Plot(d21_alra, list(F_pos_N_neg = rownames(d21_alra@meta.data %>% filter(F_pos_N_neg == TRUE))), pt.size = 0.7)
VlnPlot_scCustom(d21_alra, features = c("AR_signature1", "NE_signature1"), group.by = c("FN_status"), colors_use = polychrome_pal) + 
  plot_annotation("FOXA2 NKX2-1 Status") + theme(legend.position = "right")
VlnPlot_scCustom(d21_alra, features = c("AR_signature1"), group.by = c("FN_status"), colors_use = polychrome_pal)
VlnPlot_scCustom(d21_alra, features = c("NE_signature1"), group.by = c("FN_status"), colors_use = polychrome_pal)

dev.off()

# others

d21_add <- readRDS("../single_d21/d21_simple.RDS")
FeaturePlot_scCustom(d21_add, "FOXA2", reduction = "umap_rna_alone",
                     pt.size = 0.7, colors_use = viridis_light_high)
new_cells <- Cells(d21)[!Cells(d21) %in% Cells(d21_add)]

d21_custom <- readRDS("../d21_custom/d21_custom.RDS")
expression_d21_custom <- FetchData(d21_custom, c("FOXA2", "FOXA2-polyA", "NKX2-1"), slot = "counts")
F_pos_P_pos <-  ifelse((expression_d21_custom$FOXA2 > 0 & expression_d21_custom$`FOXA2-polyA` > 0), TRUE, FALSE)
cells_18 <- expression_d21_custom[F_pos_P_pos,] %>% rownames()
expression[cells_18,]

our_AR <- c("AR", "KLK2", "KLK3", "TMPRSS2", "NKX3-1", "FKBP5", "PLPP1", "PMEPA1", "PART1", "ALDH1A3", "STEAP4")
our_NE <- c("SYP", "NCAM1", "ENO2", "INSM1", "SOX2", "NKX2-1")

d21 <- AddModuleScore(d21, features = list(our_AR), name = "AR_signature", search = TRUE)
d21 <- AddModuleScore(d21, features = list(our_NE), name = "NE_signature", search = TRUE)
FeaturePlot_scCustom(d21, c("AR", "KLK3", "FOXA2", "NKX2-1"), reduction = "umap_rna_alone", pt.size = 1)
FeaturePlot_scCustom(d21, "AR_signature1", pt.size = 1)
FeaturePlot_scCustom(d21, "NE_signature1", pt.size = 1)
DimPlot_scCustom(d21, group.by = "F_pos_N_pos", pt.size = 1)
DimPlot_scCustom(d21, group.by = "F_pos_N_neg", pt.size = 1)
Cell_Highlight_Plot(d21, list(F_pos_N_pos = rownames(d21@meta.data %>% filter(F_pos_N_pos == TRUE))), pt.size = 2)
Cell_Highlight_Plot(d21, list(F_pos_N_neg = rownames(d21@meta.data %>% filter(F_pos_N_neg == TRUE))), pt.size = 2)

VlnPlot_scCustom(d21, features = c("AR_signature1", "NE_signature1"), group.by = c("F_pos_N_pos")) + 
  plot_annotation("FOXA2 postive NKX2-1 Positive") + theme(legend.position = "right")
VlnPlot_scCustom(d21, features = c("AR_signature1", "NE_signature1"), group.by = c("F_pos_N_neg")) + 
  plot_annotation("FOXA2 postive NKX2-1 Negative") + theme(legend.position = "right")
VlnPlot_scCustom(d21, features = c("AR_signature1", "NE_signature1"), group.by = c("FN_status"), colors_use = polychrome_pal) + 
  plot_annotation("FOXA2 NKX2-1 Status") + theme(legend.position = "right")


cell_expression <- d21@assays$RNA@data[c("MIPOL1", "ETV1", "DGKB", "AR", "FOXA2",
                                            "SYP", "RB1", "TP53", "ENO2", "VIM", "HOXB13",
                                            "NKX2-1"),] %>% as.matrix() %>%  t()

condition <- d21@meta.data

plot_data <- cbind(condition, cell_expression)
saveRDS(plot_data, "data_violin_stats.RDS")
