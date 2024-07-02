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

counts <- Read10X_h5("../arc_ranger/d14/filtered_feature_bc_matrix.h5")
fragpath <- "../arc_ranger/d14/atac_fragments.tsv.gz"

edb <- EnsDb.Hsapiens.v86
seqlevelsStyle(edb) <- "UCSC"

annotation <- GetGRangesFromEnsDb(ensdb = edb)

d14 <- CreateSeuratObject(counts = counts$`Gene Expression`, assay = "RNA", project = "d14")
d14[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,
                                      sep = c(":", "-"),
                                      fragments = fragpath,
                                      annotation = annotation)
DefaultAssay(d14) <- "ATAC"

d14 <- NucleosomeSignal(d14)
d14 <- TSSEnrichment(d14)

VlnPlot(
  object = d14,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 2,
  pt.size = 0
)
saveRDS(d14, "d14_raw.RDS")
d14 <- readRDS("d14_raw.RDS")

VlnPlot(d14, features = c("nCount_ATAC"))
  # scale_y_continuous(breaks = seq(2, 10, by = 1))
  scale_y_continuous(labels = scales::comma, breaks = seq(0,90000, by = 5000)) + 
  scale_y_log10(labels = scales::comma, breaks = seq(0,90000, by = 5000))

d14 <- subset(
  x = d14,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 30000 &
    nCount_ATAC > 20000 &
    nCount_RNA > 5000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 3
)

peaks <- CallPeaks(d14, macs2.path = "/home/vkp2256/.conda/envs/singlecell/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(d14),
  features = peaks,
  cells = colnames(d14)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
d14[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

saveRDS(macs2_counts, file = "d14_macs2_counts.RDS")


#
DefaultAssay(d14) <- "RNA"
d14 <- SCTransform(d14)
d14 <- RunPCA(d14)

DefaultAssay(d14) <- "peaks"
d14 <- FindTopFeatures(d14, min.cutoff = 5)
d14 <- RunTFIDF(d14)
d14 <- RunSVD(d14)

saveRDS(d14, "simple_d14.RDS")
d14 <- readRDS("simple_d14.RDS")

FragmentHistogram(d14)

DepthCor(d14)
DefaultAssay(d14) <- "peaks"

d14 <- RunUMAP(d14, reduction = "lsi", dims = 2:40, reduction.name = "umap_atac_alone")
d14 <- FindNeighbors(d14, reduction = "lsi", dims = 2:40)
d14 <- FindClusters(d14)

DefaultAssay(d14) <- "SCT"
d14 <- FindNeighbors(d14, reduction = "pca", dims = 1:50)
d14 <- RunUMAP(d14, reduction = "pca", dims = 1:50, reduction.name = "umap_rna_alone")

DimPlot(d14, reduction = "umap_rna_alone", group.by = "SCT_snn_res.0.8")

d14 <- NormalizeData(d14, assay = "RNA")
pdf("single_d14.pdf")
DefaultAssay(d14) <- "RNA"
DimPlot(d14, reduction = "umap_atac_alone", group.by = "peaks_snn_res.0.8")
DimPlot(d14, reduction = "umap_rna_alone", group.by = "SCT_snn_res.0.8")
FeaturePlot_scCustom(d14, c("AR", "KLK3", "FOXA2", "NCAM1"), reduction = "umap_atac_alone") + plot_annotation("ATAC")
FeaturePlot_scCustom(d14, c( "AR", "KLK3", "FOXA2", "NCAM1"), reduction = "umap_rna_alone") + plot_annotation("RNA")
dev.off()

saveRDS(d14, "simple_d14.RDS")

RidgePlot(d14, "KLK3")
klk3_low <- subset(d14, KLK3 <= 2) %>% Cells()
klk3_high <- subset(d14, KLK3 > 2) %>% Cells()
write_csv(as.data.frame(klk3_low), "klk3_low.csv", col_names = F)
write_csv(as.data.frame(klk3_high), "klk3_high.csv", col_names = F)

klk3_status <- ifelse(Cells(d14) %in% klk3_high, "high", "low")
d14$klk3_status <- klk3_status
DimPlot_scCustom(d14, group.by = "klk3_status", pt.size = 0.8)

