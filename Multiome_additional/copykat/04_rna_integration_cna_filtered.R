library(tidyverse)
library(Seurat)
library(patchwork)
set.seed(1234)
options(scipen = 999)

stepped <- DiscretePalette_scCustomize(num_colors = 24, palette = "stepped")
stepped <- stepped[c(7,6,5, 16,15,14, 12,11,10,9, 24, 24)]

source("../../../../../multiome2/atac_integration/epi_loose/meta_highlight_manual.R")
source("../../../../../multiome2/atac_integration/epi_loose/Meta_Highlight_Plot_manual.R")

cna_filtered_mat <- read_delim("../combined_cna_1_heatmap.txt")
cna_d2 <- colnames(cna_filtered_mat)[colnames(cna_filtered_mat) %>% str_starts("d2-")] %>% str_remove("d2-")
cna_d14 <- colnames(cna_filtered_mat)[colnames(cna_filtered_mat) %>% str_starts("d14-")] %>% str_remove("d14-")
cna_d21 <- colnames(cna_filtered_mat)[colnames(cna_filtered_mat) %>% str_starts("d21-")] %>% str_remove("d21-")

cna_d21_filtered <- read_delim("../d21/D21_add_copykat_CNA_results.txt")

counts_d2 <- Read10X("../../../../../multiome/arc_ranger/d2/filtered_feature_bc_matrix")
d2_all <- CreateSeuratObject(counts = counts_d2$`Gene Expression`, assay = "RNA", project = "d2")
d2 <- d2_all %>% subset(cells = cna_d2)
d2 <- NormalizeData(d2)
d2[["percent.mt"]] <- PercentageFeatureSet(d2, pattern = "^MT-")


counts_d14 <- Read10X("../../../../../multiome/arc_ranger/d14/filtered_feature_bc_matrix")
d14_all <- CreateSeuratObject(counts = counts_d14$`Gene Expression`, assay = "RNA", project = "d14")
d14 <- d14_all %>% subset(cells = cna_d14)
d14 <- NormalizeData(d14)

d14[["percent.mt"]] <- PercentageFeatureSet(d14, pattern = "^MT-")


counts_d21 <- Read10X("../../../../../multiome/arc_ranger/d21_add/filtered_feature_bc_matrix")
d21_all <- CreateSeuratObject(counts = counts_d21$`Gene Expression`, assay = "RNA", project = "d21")
d21 <- d21_all %>% subset(cells = cna_d21)
d21 <- NormalizeData(d21)
d21[["percent.mt"]] <- PercentageFeatureSet(d21, pattern = "^MT-")

# merging ----
merged <- merge(d2, y = c(d14, d21), add.cell.ids = c("D2", "D14", "D21"), merge.data = TRUE)
merge_list <- SplitObject(merged, split.by = "orig.ident")
merge_list <- lapply(X = merge_list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = merge_list, nfeatures = 3000)
merge_list <- PrepSCTIntegration(object.list = merge_list, anchor.features = features)
merge_list <- lapply(X = merge_list, FUN = RunPCA, features = features)

anchors_rpca <- FindIntegrationAnchors(object.list = merge_list, normalization.method = "SCT",
                                       anchor.features = features, k.anchor = 10, reduction = "rpca")
rna_integrated <- IntegrateData(anchorset = anchors_rpca, normalization.method = "SCT")
rna_integrated <- RunPCA(rna_integrated)
rna_integrated <- RunUMAP(rna_integrated, dims = 1:30)
rna_integrated$day <- ifelse(rna_integrated$orig.ident == "d2", "d0", rna_integrated$orig.ident)
DimPlot_scCustom(rna_integrated, group.by = "day", colors_use = ColorBlind_Pal(), pt.size = 0.7)
table(rna_integrated$orig.ident)

tree_2 <- readRDS("../heatmap_all_filtered_by_heatmap/tree_2.RDS")
names(tree_2) <- names(tree_2) %>% str_replace("-", "_")
names(tree_2) <- toupper(names(tree_2))
names(tree_2) %in% Cells(rna_integrated) %>% table()
Cells(rna_integrated) %in% names(tree_2) %>% table()
rna_integrated$tree_2 <- tree_2
DimPlot_scCustom(rna_integrated, group.by = "tree_2", colors_use = c("#FB8072", "#80B1D3"), pt.size = 1)
rna_integrated$day <- ifelse(integrated$orig.ident == "d2", "d0", integrated$orig.ident)

rna_integrated$category <- factor(paste0(toupper(rna_integrated$day), " Clone ", rna_integrated$tree_2),
                              levels = c("D0 Clone 1", "D0 Clone 2", "D14 Clone 1", "D14 Clone 2", 
                                         "D21 Clone 1", "D21 Clone 2"))

table(rna_integrated$day, rna_integrated$tree_2)
DimPlot_scCustom(rna_integrated, group.by = "category", pt.size = 1)

DefaultAssay(rna_integrated) <- "RNA"
our_AR <- c("AR", "KLK2", "KLK3", "TMPRSS2", "NKX3-1", "FKBP5", "PLPP1", "PMEPA1", "PART1", 
            "ALDH1A3", "STEAP4")
our_NE <- c("SYP", "NCAM1", "ENO2", "INSM1", "SOX2", "NKX2-1")

rna_integrated <- AddModuleScore(rna_integrated, features = list(our_AR), name = "AR_signature", search = TRUE)
rna_integrated <- AddModuleScore(rna_integrated, features = list(our_NE), name = "NE_signature", search = TRUE)


table(integrated$day, integrated$tree_2)
VlnPlot_scCustom(rna_integrated, "AR_signature1", group.by = "category", stepped[c(1,2,5,6,9,10)], pt.size = 0)
VlnPlot_scCustom(rna_integrated, "NE_signature1", group.by = "category", stepped[c(1,2,5,6,9,10)], pt.size = 0)
Meta_Highlight_Plot(rna_integrated, meta_data_column = "category", pt.size = 3, meta_data_highlight = c("D0 Clone 2"))
Meta_Highlight_Plot(rna_integrated, meta_data_column = "category", pt.size = 3, meta_data_highlight = c("D14 Clone 2"))

clone2_d14_cells <- rna_integrated[[]] %>% filter(category == "D14 Clone 2") %>% rownames()
clone2_d2_cells <- rna_integrated[[]] %>% filter(category == "D0 Clone 2") %>% rownames()

pdf("copykat_CNA_filters.pdf")
gridExtra::grid.table(data.frame(table(rna_integrated$category)))

DimPlot_scCustom(rna_integrated, group.by = "day", colors_use = ColorBlind_Pal(), pt.size = 0.7)
FeaturePlot(rna_integrated, order = T, pt.size = 0.7, features = "AR_signature1") + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 1.5), na.value = "lightgrey")
FeaturePlot(rna_integrated, order = T, pt.size = 0.7, features = "NE_signature1") +
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 0.6), na.value = "lightgrey")
VlnPlot_scCustom(rna_integrated, "AR_signature1", group.by = "category", stepped[c(1,2,5,6,9,10)], pt.size = 0)
VlnPlot_scCustom(rna_integrated, "NE_signature1", group.by = "category", stepped[c(1,2,5,6,9,10)], pt.size = 0)


Meta_Highlight_Plot_man(rna_integrated, meta_data_column = "category", pt.size = 0.7,
                        meta_data_highlight = c("D0 Clone 1", "D0 Clone 2"), c("#FB8072", "#80B1D3"),
                        order = c( "D0 Clone 2", "D0 Clone 1"))
Meta_Highlight_Plot(rna_integrated, meta_data_column = "category", pt.size = 3, meta_data_highlight = c("D0 Clone 2"))
Meta_Highlight_Plot_man(rna_integrated, meta_data_column = "category", pt.size = 0.7,
                        meta_data_highlight = c("D14 Clone 1", "D14 Clone 2"), c("#FB8072", "#80B1D3"),
                        order = c( "D14 Clone 2", "D14 Clone 1"))
Meta_Highlight_Plot(rna_integrated, meta_data_column = "category", pt.size = 3, meta_data_highlight = c("D14 Clone 2"))
Meta_Highlight_Plot_man(rna_integrated, meta_data_column = "category", pt.size = 0.7,
                        meta_data_highlight = c("D21 Clone 1", "D21 Clone 2"), c("#80B1D3", "#FB8072"),
                        order = c( "D21 Clone 2", "D21 Clone 1"))


dev.off()
