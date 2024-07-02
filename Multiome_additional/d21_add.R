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

counts_d21 <- Read10X_h5("../../multiome/arc_ranger/d21_add/filtered_feature_bc_matrix.h5")
fragpath_d21 <- "../../multiome/arc_ranger/d21/atac_fragments.tsv.gz"

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

v <- VlnPlot(
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

v1 + v2 + v1_1 + v2_1 + v3 + v4 + plot_layout(nrow = 3, ncol = 2)

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

DefaultAssay(d21) <- "peaks"
d21 <- FindTopFeatures(d21, min.cutoff = 5)
d21 <- RunTFIDF(d21)
d21 <- RunSVD(d21)

DepthCor(d21)
DefaultAssay(d21) <- "peaks"

d21 <- RunUMAP(d21, reduction = "lsi", dims = 2:40, reduction.name = "umap_atac_alone")
d21 <- FindNeighbors(d21, reduction = "lsi", dims = 2:40)
d21 <- FindClusters(d21)
DimPlot(d21)

DefaultAssay(d21) <- "SCT"
d21 <- FindNeighbors(d21, reduction = "pca", dims = 1:50)
d21 <- RunUMAP(d21, reduction = "pca", dims = 1:50, reduction.name = "umap_rna_alone")
d21 <- FindClusters(d21)
DimPlot_scCustom(d21, reduction = "umap_rna_alone")

DefaultAssay(d21) <- "RNA"
d21 <- NormalizeData(d21)

pdf("single_d21.pdf")
FeaturePlot_scCustom(d21, c("nCount_RNA"), reduction = "umap_rna_alone")
FeaturePlot_scCustom(d21, c("nCount_ATAC"), reduction = "umap_atac_alone")
FeaturePlot_scCustom(original_d21, c("nCount_RNA"), reduction = "umap_rna_alone")
FeaturePlot_scCustom(original_d21, c("nCount_ATAC"), reduction = "umap_atac_alone")
FeaturePlot_scCustom(original_d2, c("nCount_RNA"), reduction = "umap_rna_alone")
FeaturePlot_scCustom(original_d2, c("nCount_ATAC"), reduction = "umap_atac_alone")
FeaturePlot_scCustom(d21, c("KLK3", "AR", "FOXA2", "NKX2-1"), reduction = "umap_atac_alone") + plot_annotation("ATAC")
FeaturePlot_scCustom(d21, c("KLK3", "AR", "FOXA2", "NKX2-1"), reduction = "umap_rna_alone") + plot_annotation("RNA")
FeaturePlot_scCustom(d21, c("SYP", "ASCL1", "ENO2", "NCAM1"), reduction = "umap_rna_alone") + plot_annotation("RNA")
FeaturePlot_scCustom(d21, c("SYP", "ASCL1", "ENO2", "NCAM1"), reduction = "umap_atac_alone") + plot_annotation("ATAC")

FeaturePlot(d21, c("KLK3", "AR", "FOXA2", "NKX2-1"), reduction = "umap_rna_alone", slot = "counts") + plot_annotation("RNA Raw Counts")
dev.off()
saveRDS(d21, "simple_d21.RDS")

DefaultAssay(d21) <- "RNA"
FeatureScatter(d21, "FOXA2", "NKX2-1", group.by = 'orig.ident')
FeatureScatter(d21, "KLK3", "AR", group.by = 'orig.ident')

expression <- FetchData(d21, c("SYP", "FOXA2", "ENO2", "NCAM1", "NKX2-1"), slot = "counts")
expression_norm <- FetchData(d21, c("SYP", "FOXA2", "ENO2", "NCAM1", "NKX2-1"))

F_pos_N_pos <- ifelse((expression$FOXA2 > 0 & expression$`NKX2-1` > 0), TRUE, FALSE)
F_pos_N_neg <- ifelse(expression$FOXA2 > 0 & expression$`NKX2-1` == 0, TRUE, FALSE)

sum(expression$`NKX2-1` > 0)
sum(expression$FOXA2 > 0)
expression[F_pos_N_pos,]

d21@meta.data[!d21@meta.data$seurat_clusters %in% c(2,7) & expression$`FOXA2` > 0,] %>% dim()

