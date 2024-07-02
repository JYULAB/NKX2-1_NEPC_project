library(Seurat)
library(Signac)
library(tidyverse)
library(patchwork)
library(scCustomize)
library(EnsDb.Hsapiens.v86)
set.seed(1234)

polychrome_pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
polychrome_pal <- polychrome_pal[c(1, 3, 6:36)]

# Reading in data from Cheng et. al -----
# using CellRanger aggr output
# to stay consistent with Jiaoti's paper, otherwise many more cells will be included
cheng <- Read10X("../../scrna_seq/integration/cellranger/cellranger/cheng_stat/outs/count/filtered_feature_bc_matrix/")
all_pat <- CreateSeuratObject(counts = cheng, project = "patient", min.cells = 3, min.features = 200)

sample <- rownames(all_pat@meta.data) %>% str_extract("[\\d]+")
stat <- read_csv("../../scrna_seq/integration/cellranger/cellranger/aggr_stat.csv") %>% dplyr::select(status) %>% deframe()
samples <- sapply(sample, switch, "1"=stat[1], "2"=stat[2], "3"=stat[3], "4"=stat[4], "5"=stat[5],
                  "6"=stat[6], "7"=stat[7], "8"=stat[8], "9"=stat[9], "10"=stat[10], "11"=stat[11])
all_pat$orig.ident <- as.factor(samples)
all_pat[["per_mt"]] <- PercentageFeatureSet(all_pat, pattern = "^MT-")

pat <- subset(all_pat, subset = nFeature_RNA <= 8000 & nFeature_RNA >= 500 & per_mt < 10) # matching the description in Cheng et. al.
pat <- NormalizeData(pat, normalization.method = "LogNormalize", scale.factor = 10000)
tvb <- if_else(str_sub(pat$orig.ident, -1) == "B", "Adjacent", "Tumor")
pat$tvb <- tvb
pat$dataset <- as.factor(rep("Cheng", 25425))

pat <- SCTransform(pat)
pat <- RunPCA(pat)
pat <- RunUMAP(pat, dims = 1:30, reduction = "pca")
pat <- FindNeighbors(pat, dim = 1:30)
pat <- FindClusters(pat, resolution = 0.3)
DimPlot_scCustom(pat, group.by = "SCT_snn_res.0.3")
DimPlot_scCustom(pat, group.by = "orig.ident")

# plotting markers mentioned in Cheng et. al. to assign cell types
DefaultAssay(pat) <- "RNA"
FeaturePlot_scCustom(pat, c("SYP"))
FeaturePlot_scCustom(pat, c("AR"))
FeaturePlot_scCustom(pat, c("ITGA6"))
FeaturePlot_scCustom(pat, c("KLK3"))
FeaturePlot_scCustom(pat, c("TMPRSS2"))
FeaturePlot_scCustom(pat, c("CLDN5"))
FeaturePlot_scCustom(pat, c("DCN"))
FeaturePlot_scCustom(pat, c("PTPRC"))

lineage <- c("0" = "Basal", "1" = "PSA High", "2" = "Basal", "3" = "PSA High",
             "4" = "mCRPC", "5" = "PSA Low", "6" = "PSA High","7" = "SNSC",
             "8" = "PSA Low", "9" = "CRPC AR High",
             "10" = "non-epi", "11" = "non-epi", "12" = "non-epi", "13" = "non-epi")
Idents(pat) <- "SCT_snn_res.0.3"
pat <- RenameIdents(pat, lineage)
pat$status <- Idents(pat)
DimPlot_scCustom(pat, group.by = "status", colors_use = polychrome_pal)

saveRDS(pat, "pat_SCT.RDS")
pat <- readRDS("pat_SCT.RDS")

# merging with scRNA-seq FOXA2 OE as a group -----

foxa2 <- readRDS("../../multiome/rna_integration_3/rna_integrated_rpca_default.RDS") # integrated scRNA-seq object
samples_list <- list(pat, foxa2)
samples_list <- lapply(samples_list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = samples_list, nfeatures = 3000)
samples_list <- PrepSCTIntegration(object.list = samples_list, anchor.features = features)
samples_list <- lapply(X = samples_list, FUN = RunPCA, features = features)
anchors_rpca <- FindIntegrationAnchors(object.list = samples_list, normalization.method = "SCT",
                                  anchor.features = features, reduction = "rpca",  k.anchor = 30)
integrated_rpca <- IntegrateData(anchorset = anchors_rpca, normalization.method = "SCT")

integrated_rpca <- RunPCA(integrated_rpca)
integrated_rpca <- RunUMAP(integrated_rpca, dims = 1:30)
DimPlot_scCustom(integrated_rpca, group.by = "dataset")
integrated_rpca <- FindNeighbors(integrated_rpca, dim = 1:30)
integrated_rpca <- FindClusters(integrated_rpca, resolution = 0.3)

dataset <- c(rep("Patient", 25425), rep("Day 2", 2110), rep("Day 14", 1739), rep("Day 21", 2344))
dataset <- factor(dataset, levels = c("Day 2", "Day 14", "Day 21", "Patient"))
names(dataset) <- Cells(integrated_rpca)
integrated_rpca$dataset <- dataset

new_clusts <- c(as.character(pat$status), rep("Day 2", 2110), rep("Day 14", 1739), rep("Day 21", 2344))
integrated_rpca$new_clusts <- new_clusts
Idents(integrated_rpca) <- new_clusts

DimPlot_scCustom(integrated_rpca, group.by = "new_clusts")

# subset for luminal clusters only ------
library(monocle3)
library(SeuratWrappers)
integrated_rpca$batch <- c(rep("patient", 25425), rep("cellline", 6193))
integrated_rpca$pat_clusters <- factor(pat$SCT_snn_res.0.3)
DimPlot_scCustom(integrated_rpca, group.by = "pat_clusters")


# seurat subset Fig S3D------
Idents(integrated_rpca) <- "new_clusts"
sub_integrated <- subset(integrated_rpca, idents = c("non-epi", "Basal"), invert = TRUE)
sub_integrated <- RunUMAP(sub_integrated, dims = 1:30)
saveRDS(sub_integrated, "sub_integrated.RDS")

# monocle on umap
library(SeuratWrappers)
library(monocle3)

# by sample
simple <- c("CRPC1" = "CRPC", "CRPC2" = "NEPC", "d14" = "D14", "d2" = "D0", "d21" = "D21", "CRPC3" = "CRPC", 
            "PC1_T" = "PC", "PC2_B" = "PC", "PC2_T" ="PC", "PC3_B" = "PC", "PC3_T" = "PC")
Idents(sub_integrated) <- "orig.ident"
sub_integrated <- RenameIdents(sub_integrated, simple)
sub_integrated$simplified <- factor(Idents(sub_integrated), levels = c("D0", "D14", "D21", "PC", "CRPC", "NEPC"))
DimPlot_scCustom(sub_integrated, group.by = "simplified", colors_use = ColorBlind_Pal())
saveRDS(sub_integrated, "sub_integrated.RDS")

# by cell types
simple_by_cell <- c("CRPC AR High" = "CRPC", "Day 14" = "D14", "Day 2" = "D0", "Day 21" = "D21",
                    "mCRPC" = "CRPC", "PSA High" = "PSA", "PSA Low" = "PSA", "SNSC" = "NEPC")

Idents(sub_integrated) <- "new_clusts"
sub_integrated <- RenameIdents(sub_integrated, simple_by_cell)
sub_integrated$simplified_ct <- factor(Idents(sub_integrated), levels = c("D0", "D14", "D21", "PSA", "CRPC", "NEPC"))
DimPlot_scCustom(sub_integrated, group.by = "simplified_ct", colors_use = polychrome_pal, pt.size = 0.7)

# by_combo
status <- cbind(sub_integrated$orig.ident, sub_integrated$new_clusts)
status <- as.data.frame(status)

status$annotation <- case_when(
  status$V1 == "CRPC2" & status$V2 == "PSA High" ~ "CRPC",
  status$V1 == "CRPC2" & status$V2 == "PSA Low" ~ "CRPC",
  status$V1 == "CRPC1" & status$V2 == "PSA High" ~ "CRPC",
  status$V1 == "CRPC1" & status$V2 == "PSA Low" ~ "CRPC",
  status$V1 == "CRPC3" & status$V2 == "PSA High" ~ "CRPC",
  status$V1 == "CRPC3" & status$V2 == "PSA Low" ~ "CRPC",
  # "PC1_T" ~ "PC", "PC2_B" ~ "PC", "PC2_T" ~"PC", "PC3_B" ~ "PC", "PC3_T" ~ "PC",
  .default = status$V2
)

anno <- factor(status$annotation)
names(anno) <- rownames(status)

anno <- case_when(
  anno == "CRPC AR High" ~ "CRPC",
  anno == "mCRPC" ~ "CRPC",
  anno == "PSA High" ~ "PC",
  anno == "PSA Low" ~ "PC",
  anno == "Day 2" ~ "D0",
  anno == "Day 14" ~ "D14",
  anno == "Day 21" ~ "D21",
  .default = anno
)
anno <- factor(anno, levels = c("D0", "D14", "D21", "PC", "CRPC", "SNSC"))
sub_integrated$anno <- anno

# pseudotime analysis with monocle3 Fig S3E --------
DefaultAssay(sub_integrated) <- "RNA"
lum <- as.cell_data_set(sub_integrated)
lum <- cluster_cells(lum)
plot_cells(lum)
lum <- learn_graph(lum, use_partition = FALSE)
plot_cells(lum) + scale_color_discrete(polychrome_pal)
lum <- order_cells(lum)
plot_cells(lum, color_cells_by = "pseudotime")
polychrome_pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
polychrome_pal <- polychrome_pal[c(1, 3, 6:36)]
saveRDS(lum, "lum_fin.RDS")


