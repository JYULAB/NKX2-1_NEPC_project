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

counts_d2 <- Read10X_h5("../arc_ranger/d2/filtered_feature_bc_matrix.h5")
fragpath_d2 <- "../arc_ranger/d2/atac_fragments.tsv.gz"

edb <- EnsDb.Hsapiens.v86
seqlevelsStyle(edb) <- "UCSC"

annotation <- GetGRangesFromEnsDb(ensdb = edb)

d2 <- CreateSeuratObject(counts = counts_d2$`Gene Expression`, assay = "RNA", project = "d2")
d2[["ATAC"]] <- CreateChromatinAssay(counts = counts_d2$Peaks,
                                     sep = c(":", "-"),
                                     fragments = fragpath_d2,
                                     annotation = annotation)
DefaultAssay(d2) <- "ATAC"

d2 <- NucleosomeSignal(d2)
d2 <- TSSEnrichment(d2)

VlnPlot(
  object = d2,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

VlnPlot(d2, features = c("TSS.enrichment"), y.max = 5) +
scale_y_continuous(labels = scales::comma, breaks = seq(0, 160000, by = 20000)) +
  geom_vline(xintercept = 10000)

d2 <- subset(
  x = d2,
  subset = nCount_ATAC < 150000 &
    nCount_RNA < 30000 &
    nCount_ATAC > 20000 &
    nCount_RNA > 5000 &
    nucleosome_signal < 2.5 &
    TSS.enrichment > 1
)

peaks_d2 <- CallPeaks(d2, macs2.path = "/home/vkp2256/.conda/envs/singlecell/bin/macs2")
peaks_d2 <- keepStandardChromosomes(peaks_d2, pruning.mode = "coarse")
peaks_d2 <- subsetByOverlaps(x = peaks_d2, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts_d2 <- FeatureMatrix(
  fragments = Fragments(d2),
  features = peaks_d2,
  cells = colnames(d2)
)

saveRDS(macs2_counts_d2, "macs2_counts")

d2[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts_d2,
  fragments = fragpath_d2,
  annotation = annotation
)

DefaultAssay(d2) <- "RNA"
d2 <- SCTransform(d2)
d2 <- RunPCA(d2)

DefaultAssay(d2) <- "peaks"
d2 <- FindTopFeatures(d2, min.cutoff = 5)
d2 <- RunTFIDF(d2)
d2 <- RunSVD(d2)

DepthCor(d2)
DefaultAssay(d2) <- "peaks"

d2 <- RunUMAP(d2, reduction = "lsi", dims = 2:40, reduction.name = "umap_atac_alone")
d2 <- FindNeighbors(d2, reduction = "lsi", dims = 2:40)
d2 <- FindClusters(d2)
DimPlot(d2)

DefaultAssay(d2) <- "SCT"
d2 <- FindNeighbors(d2, reduction = "pca", dims = 1:50)
d2 <- RunUMAP(d2, reduction = "pca", dims = 1:50, reduction.name = "umap_rna_alone")
d2 <- FindClusters(d2)


d2 <- NormalizeData(d2, assay = "RNA")
pdf("single_d2.pdf")
DefaultAssay(d2) <- "RNA"
DimPlot(d2, reduction = "umap_atac_alone", group.by = "peaks_snn_res.0.8")
DimPlot(d2, reduction = "umap_rna_alone", group.by = "SCT_snn_res.0.8")
FeaturePlot(d2, c("KLK3", "SYP", "FOXA2", "AR"), reduction = "umap_atac_alone")
FeaturePlot(d2, c("KLK3", "SYP", "FOXA2", "AR"), reduction = "umap_rna_alone")
FeaturePlot(d2, c("NKX2-1", "ASCL1", "ENO2", "NCAM1"), reduction = "umap")
dev.off()


saveRDS(d2, "simple_d2.RDS")

