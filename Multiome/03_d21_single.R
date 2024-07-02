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

counts_d21 <- Read10X_h5("../arc_ranger/d21/filtered_feature_bc_matrix.h5")
fragpath_d21 <- "../arc_ranger/d21/atac_fragments.tsv.gz"

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
VlnPlot(d21, features = c("nCount_RNA"), y.max = 35000)
# scale_y_continuous(breaks = seq(2, 10, by = 1))
scale_y_continuous(labels = scales::comma, breaks = seq(0, 160000, by = 20000)) +
  geom_vline(xintercept = 10000)
# scale_y_log10(, breaks = seq(0,200000, by = 10000))

d21 <- subset(
  x = d21,
  subset = nCount_ATAC < 75000 &
    nCount_RNA < 35000 &
    nCount_ATAC > 5000 &
    nCount_RNA > 5000 &
    nucleosome_signal < 1.5 &
    TSS.enrichment > 2
)

peaks_d21 <- CallPeaks(d21, macs2.path = "/home/vkp2256/.conda/envs/singlecell/bin/macs2")
peaks_d21 <- keepStandardChromosomes(peaks_d21, pruning.mode = "coarse")
peaks_d21 <- subsetByOverlaps(x = peaks_d21, ranges = blacklist_hg38_unified, invert = TRUE)

macs2_counts_d21 <- FeatureMatrix(
  fragments = Fragments(d21),
  features = peaks_d21,
  cells = colnames(d21)
)

saveRDS(macs2_counts_d21, "macs2_counts")

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
DimPlot(d21)

d21 <- NormalizeData(d21, assay = "RNA")
pdf("single_d21.pdf")
DefaultAssay(d21) <- "RNA"
DimPlot(d21, reduction = "umap_atac_alone", group.by = "peaks_snn_res.0.8")
DimPlot(d21, reduction = "umap_rna_alone", group.by = "SCT_snn_res.0.8", label = TRUE, label.box = TRUE)
FeaturePlot(d21, c("KLK3", "SYP", "FOXA2", "AR"), reduction = "umap_atac_alone")
FeaturePlot(d21, c("KLK3", "SYP", "FOXA2", "AR"), reduction = "umap_rna_alone")
FeaturePlot(d21, c("NKX2-1", "ASCL1", "ENO2", "NCAM1"), reduction = "umap_rna_alone")
FeaturePlot(d21, c("NKX2-1", "ASCL1", "ENO2", "NCAM1"), reduction = "umap_atac_alone")
dev.off()

saveRDS(d21, "simple_d21.RDS")
