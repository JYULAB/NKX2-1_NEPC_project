#enndogenous d21

library(tidyverse)
library(Signac)
library(Seurat)
library(patchwork)
set.seed(1234)
library(scCustomize)

library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

d14 <- readRDS("../../multiome/single/simple_d14.RDS")
d2 <- readRDS("../../multiome/single_d2/simple_d2.RDS")
d21 <- readRDS("../single_d21/simple_d21.RDS")

DefaultAssay(d2) <- "RNA"
d2 <- NormalizeData(d2)
DefaultAssay(d14) <- "RNA"
d14 <- NormalizeData(d14)

# merge <- merge(d2, y = c(d14, d21), add.cell.ids = c("d2", "d14", "d21"), merge.data = TRUE)
saveRDS(merge, "simple_merge.RDS")


merge_list <- SplitObject(merge, split.by = "orig.ident")
merge_list <- lapply(X = merge_list, FUN = SCTransform, return.only.var.genes = FALSE)
features <- SelectIntegrationFeatures(object.list = merge_list, nfeatures = 3000)
saveRDS(features, "features.RDS")
saveRDS(features, "features_3k.RDS")

merge_list <- PrepSCTIntegration(object.list = merge_list, anchor.features = features)

merge_list <- lapply(X = merge_list, FUN = RunPCA, features = features)
anchors <- FindIntegrationAnchors(object.list = merge_list, normalization.method = "SCT",
                                  anchor.features = features, reduction = "rpca")
rna_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
rna_integrated <- RunPCA(rna_integrated)
rna_integrated <- RunUMAP(rna_integrated, dims = 1:30)
rna_integrated <- FindNeighbors(rna_integrated, dims = 1:30)
rna_integrated <- FindClusters(rna_integrated, resolution = 0.8)
DimPlot_scCustom(rna_integrated, group.by = "orig.ident", colors_use = ColorBlind_Pal())
saveRDS(rna_integrated, "rna_integrated_3k.RDS")

DefaultAssay(rna_integrated) <- "RNA"
our_AR <- c("AR", "KLK2", "KLK3", "TMPRSS2", "NKX3-1", "FKBP5", "PLPP1", "PMEPA1", "PART1", 
            "ALDH1A3", "STEAP4")
our_NE <- c("SYP", "NCAM1", "ENO2", "INSM1", "SOX2", "NKX2-1")
rna_integrated <- AddModuleScore(rna_integrated, features = list(our_AR), name = "AR_signature", search = TRUE)
rna_integrated <- AddModuleScore(rna_integrated, features = list(our_NE), name = "NE_signature", search = TRUE)

FeaturePlot(rna_integrated, order = T, pt.size = 0.7, features = "AR_signature1") + 
  # labs(x = "atac_UMAP 1", y = "atac_UMAP 2", title = "AR Signature (RNA)")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 1), na.value = "lightgrey")
FeaturePlot(rna_integrated, order = T, pt.size = 0.7, features = "NE_signature1") +
  # labs(x = "atac_UMAP 1", y = "atac_UMAP 2", title = "NE Signature (RNA)")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 0.5), na.value = "lightgrey")
VlnPlot_scCustom(integrated, "AR_signature1", group.by = "dataset")

# copyKat 2 -----
copykat_2 <- readRDS("../../multiome_clean/copykat/unfiltered/heatmap/tree2.RDS")
names(copykat_2) <- names(copykat_2) %>% str_replace("-", "_")
stepped <- DiscretePalette_scCustomize(num_colors = 24, palette = "stepped")
stepped <- stepped[c(7,6,5, 16,15,14,  12,11,10,9, 24, 24)]

rna_integrated$rna_clones <- copykat_2
DimPlot_scCustom(rna_integrated, group.by = "rna_clones", colors_use = ColorBlind_Pal(), pt.size = 1)

rna_integrated$day <- ifelse(rna_integrated$orig.ident == "d2", "d0", rna_integrated$orig.ident)
rna_integrated$category <- factor(paste0(toupper(rna_integrated$day), " Clone ", rna_integrated$rna_clones),
                                  levels = c("D0 Clone 1", "D0 Clone 2", "D14 Clone 1", "D14 Clone 2", 
                                             "D21 Clone 1", "D21 Clone 2", "D21 Clone NA"))
DimPlot_scCustom(rna_integrated, group.by = "category", pt.size = 1, stepped[c(1,2,5,6,9,10,12)])

# pdfs --------
source("../atac_integration/epi_loose/Meta_Highlight_Plot_manual.R")

pdf("copykat_rpca.pdf")
gridExtra::grid.table(data.frame(table(rna_integrated$category)))

DimPlot_scCustom(rna_integrated, group.by = "day", colors_use = ColorBlind_Pal(), pt.size = 0.7)
FeaturePlot(rna_integrated, order = T, pt.size = 0.7, features = "AR_signature1") + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 1.5), na.value = "lightgrey")
FeaturePlot(rna_integrated, order = T, pt.size = 0.7, features = "NE_signature1") +
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 0.6), na.value = "lightgrey")
VlnPlot_scCustom(rna_integrated, "AR_signature1", group.by = "category", stepped[c(1,2,5,6,9,10,12)], pt.size = 0)
VlnPlot_scCustom(rna_integrated, "NE_signature1", group.by = "category", stepped[c(1,2,5,6,9,10,12)], pt.size = 0)

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




