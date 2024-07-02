# atac endogenous d21 add

library(tidyverse)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(future)
library(stringr)
library(patchwork)
library(scCustomize)
library(EnsDb.Hsapiens.v86)
set.seed(1234)

polychrome_pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
polychrome_pal <- polychrome_pal[c(1, 3, 6:36)]

# processing with d21 add endogenous -------
d21 <- readRDS("../single_d21/simple_d21.RDS")
d14 <- readRDS("../../multiome/single/simple_d14.RDS")
d2 <- readRDS("../../multiome/single_d2/simple_d2.RDS")

DefaultAssay(d14) <- "ATAC"
DefaultAssay(d2) <- "ATAC"
DefaultAssay(d21) <- "ATAC"

mac_d14 <- readRDS("../../multiome/single/d14_macs2_counts.RDS")
mac_d2 <- readRDS("../../multiome/single_d2/macs2_counts")
mac_d21 <- readRDS("../single_d21/macs2_counts.RDS")

peaks_d14 <- str_split(mac_d14@Dimnames[[1]], "-", simplify = TRUE) %>% as.data.frame()
colnames(peaks_d14) <- c("chr", "start", "end")

peaks_d2 <- str_split(mac_d2@Dimnames[[1]], "-", simplify = TRUE) %>% as.data.frame()
colnames(peaks_d2) <- c("chr", "start", "end")

peaks_d21 <- str_split(mac_d21@Dimnames[[1]], "-", simplify = TRUE) %>% as.data.frame()
colnames(peaks_d21) <- c("chr", "start", "end")

gr_d14 <- makeGRangesFromDataFrame(peaks_d14)
gr_d2 <- makeGRangesFromDataFrame(peaks_d2)
gr_d21 <- makeGRangesFromDataFrame(peaks_d21)

combined_peaks <- reduce(x = c(gr_d14, gr_d2, gr_d21))
peakwidths <- width(combined_peaks)
range(peakwidths)

tmp_frag <- Fragments(d14)
Fragments(d14) <- NULL
tmp_frag[[1]] <- UpdatePath(tmp_frag[[1]], "../../multiome/arc_ranger/d14/atac_fragments.tsv.gz")
Fragments(d14) <- tmp_frag

d14_counts <- FeatureMatrix(
  fragments = Fragments(d14),
  features = combined_peaks,
  cells = colnames(d14)
)
saveRDS(d14_counts, "d14_counts.RDS")

tmp_frag <- Fragments(d2)
Fragments(d2) <- NULL
tmp_frag[[1]] <- UpdatePath(tmp_frag[[1]], "../../multiome/arc_ranger/d2/atac_fragments.tsv.gz")
Fragments(d2) <- tmp_frag

d2_counts <- FeatureMatrix(
  fragments = Fragments(d2),
  features = combined_peaks,
  cells = colnames(d2)
)
saveRDS(d2_counts, "d2_counts.RDS")

d21_counts <- FeatureMatrix(
  fragments = Fragments(d21),
  features = combined_peaks,
  cells = colnames(d21)
)
saveRDS(d21_counts, "d21_counts.RDS")

d14_assay <- CreateChromatinAssay(d14_counts, fragments = Fragments(d14))
d14_new <- CreateSeuratObject(d14_assay, assay = "ATAC")
d14_new[["RNA"]] <- d14[["RNA"]]
d14_new[["SCT"]] <- d14[["SCT"]]
d14_new$dataset <- "d14"
saveRDS(d14_new, "d14.RDS")

d2_assay <- CreateChromatinAssay(d2_counts, fragments = Fragments(d2))
d2_new <- CreateSeuratObject(d2_assay, assay = "ATAC")
d2_new[["RNA"]] <- d2[["RNA"]]
d2_new[["SCT"]] <- d2[["SCT"]]
d2_new$dataset <- "d2"
saveRDS(d2_new, "d2.RDS")

d21_assay <- CreateChromatinAssay(d21_counts, fragments = Fragments(d21))
d21_new <- CreateSeuratObject(d21_assay, assay = "ATAC")
d21_new[["RNA"]] <- d21[["RNA"]]
d21_new[["SCT"]] <- d21[["SCT"]]
d21_new$dataset <- "d21"
saveRDS(d21_new, "d21.RDS")

d2_new <- RenameCells(d2_new, add.cell.id = "d2")
d21_new <- RenameCells(d21_new, add.cell.id = "d21")
d14_new <- RenameCells(d14_new, add.cell.id = "d14")

d2_new <- FindTopFeatures(d2_new, min.cutoff = 10)
d2_new <- RunTFIDF(d2_new)
d2_new <- RunSVD(d2_new)
saveRDS(d2_new, "d2.RDS")

d21_new <- FindTopFeatures(d21_new, min.cutoff = 10)
d21_new <- RunTFIDF(d21_new)
d21_new <- RunSVD(d21_new)
saveRDS(d21_new, "d21.RDS")

d14_new <- FindTopFeatures(d14_new, min.cutoff = 10)
d14_new <- RunTFIDF(d14_new)
d14_new <- RunSVD(d14_new)
saveRDS(d14_new, "d14.RDS")

# integration of timepoints -----
combined <- merge(
  x = d2_new,
  y = list(d14_new, d21_new)
)
saveRDS(combined, "combined_plain.RDS")
combined[["ATAC"]]
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
DimPlot(combined, group.by = "dataset") + ggtitle("ATAC Combined")
saveRDS(combined, "combined.RDS")

integration.anchors <- FindIntegrationAnchors(
  object.list = list(d2_new, d14_new, d21_new),
  anchor.features = rownames(d2_new),
  reduction = "rlsi",
  dims = 2:50
)

integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:50
)

integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
DimPlot_scCustom(integrated, group.by = "dataset", colors_use = ColorBlind_Pal(), pt.size = 1)
integrated <- FindNeighbors(integrated, reduction = "integrated_lsi", dims = 2:30)
integrated <- FindClusters(integrated, graph.name = "ATAC_snn")
integrated <- FindClusters(integrated, graph.name = "ATAC_snn", resolution = 0.3)
integrated$dataset <- factor(integrated$dataset, levels = c("d2", "d14", "d21"))
DimPlot_scCustom(integrated, group.by = "dataset", colors_use = ColorBlind_Pal(), pt.size = 1)
saveRDS(integrated, "atac_integrated_new.RDS")

cell_types <- c("0" = "Luminal", "1" = "NEPC", "2" = "Luminal", "3" = "NEPC", "4" = "NEPC")
Idents(integrated) <- "ATAC_snn_res.0.3"
integrated <- RenameIdents(integrated, cell_types)
integrated$cell_types <- Idents(integrated)
DimPlot_scCustom(integrated, group.by = "cell_types")

# RNA ----
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated)

# annotation
library(BSgenome.Hsapiens.UCSC.hg38)

DefaultAssay(integrated) <- "ATAC"
edb <- EnsDb.Hsapiens.v86
seqlevelsStyle(edb) <- "UCSC"

annotation <- GetGRangesFromEnsDb(ensdb = edb)
integrated@assays$ATAC@annotation <- annotation
integrated <- RegionStats(integrated, genome = BSgenome.Hsapiens.UCSC.hg38)

# RNA signature -----
DefaultAssay(integrated) <- "RNA"
our_AR <- c("AR", "KLK2", "KLK3", "TMPRSS2", "NKX3-1", "FKBP5", "PLPP1", "PMEPA1", "PART1", 
            "ALDH1A3", "STEAP4")
our_NE <- c("SYP", "NCAM1", "ENO2", "INSM1", "SOX2", "NKX2-1")

integrated <- AddModuleScore(integrated, features = list(our_AR), name = "AR_signature", search = TRUE)
integrated <- AddModuleScore(integrated, features = list(our_NE), name = "NE_signature", search = TRUE)

FeaturePlot(integrated, order = T, pt.size = 0.7, features = "AR_signature1") + 
  labs(x = "atac_UMAP 1", y = "atac_UMAP 2", title = "AR Signature (RNA)")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 1), na.value = "lightgrey")
FeaturePlot(integrated, order = T, pt.size = 0.7, features = "NE_signature1") +
  labs(x = "atac_UMAP 1", y = "atac_UMAP 2", title = "NE Signature (RNA)")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 0.5), na.value = "lightgrey")
VlnPlot_scCustom(integrated, "AR_signature1", group.by = "dataset")

Cell_Highlight_Plot(integrated, cells_highlight = list(AR_sig_pos = rownames(integrated@meta.data %>% filter(AR_signature1 > 0))))


pdf("scores.pdf")
FeaturePlot(integrated, order = T, pt.size = 0.7, features = "AR_signature1") + 
  labs(x = "atac_UMAP 1", y = "atac_UMAP 2", title = "AR Signature (RNA)")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 1), na.value = "lightgrey")
FeaturePlot(integrated, order = T, pt.size = 0.7, features = "NE_signature1") +
  labs(x = "atac_UMAP 1", y = "atac_UMAP 2", title = "NE Signature (RNA)")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 0.5), na.value = "lightgrey")

dev.off()

saveRDS(integrated, "atac_integrated_meta.RDS")
