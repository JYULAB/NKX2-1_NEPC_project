library(tidyverse)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(future)
library(stringr)
library(patchwork)
library(scCustomize)
set.seed(1234)

polychrome_pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
polychrome_pal <- polychrome_pal[c(1, 3, 6:36)]

integrated <- readRDS("../atac_integrated_meta.RDS")
epi_2 <- readRDS("../../../multiome_clean/epiAneufinder/endogenous_loose/analysis/depth2.RDS")
epi_clones_2 <- epi_2[[2]]$annot
names(epi_clones_2) <- paste0(str_extract(epi_2[[2]]$cell, "d\\d{1,}"), "_", str_extract(epi_2[[2]]$cell, "(?<=l-)[:graph:]*(?=-d)"))
integrated$atac_clones_2 <- epi_clones_2
table(integrated$dataset, integrated$atac_clones_2)
DimPlot_scCustom(integrated, group.by = "atac_clones_2", colors_use = c("#F8766D", "#00BFC4"), pt.size = 1)
prop_clones_2 <- prop.table(table(integrated$dataset, integrated$atac_clones_2), margin = 1)

table(integrated$dataset, integrated$atac_clones_2, useNA = "ifany")

d21_c2 <- integrated@meta.data %>% dplyr::filter(dataset == "d21" & atac_clones_2 == "Clone 2")
Cell_Highlight_Plot(integrated, cells_highlight = list(D21_C2 = rownames(d21_c2)), pt.size = 2)

integrated$cat <- factor(paste0(toupper(integrated$dataset), "_", integrated$atac_clones_2, "_", integrated$cell_types),
                         levels = c("D2_Clone 1_Luminal", "D2_Clone 1_NEPC", "D2_Clone 2_Luminal",
                                    "D14_Clone 1_Luminal", "D14_Clone 1_NEPC", "D14_Clone 2_Luminal",
                                    "D21_Clone 1_Luminal", "D21_Clone 1_NEPC", "D21_Clone 2_Luminal", "D21_Clone 2_NEPC", 
                                    "D21_NA_Luminal", "D21_NA_NEPC"))

integrated$category <- factor(paste0(toupper(integrated$dataset), " ", integrated$atac_clones_2),
                              levels = c("D2 Clone 1", "D2 Clone 2", "D14 Clone 1", "D14 Clone 2", 
                                         "D21 Clone 1", "D21 Clone 2", "D21 NA"))
DimPlot_scCustom(integrated, group.by = "category", pt.size = 1)

stepped <- DiscretePalette_scCustomize(num_colors = 24, palette = "stepped")
stepped <- stepped[c(7,6,5, 16,15,14,  12,11,10,9, 24, 24)]

library(gridExtra)

pdf("epiAneufinder_clones_2.pdf")
grid.table(t(table(integrated$dataset, integrated$atac_clones_2)))
ggplot(as.data.frame(prop_clones_2), aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = "fill") +
  labs(x= "Timepoint", y = "Proportion") +
  scale_fill_discrete(name = "Clones") +
  theme_classic(base_size = 18)

ggplot(as.data.frame(prop_clones_2), aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = "fill", width = 0.5) +
  labs(x= "Timepoint", y = "Proportion") +
  scale_fill_discrete(name = "Clones") +
  theme_classic(base_size = 18)
Cell_Highlight_Plot(integrated, cells_highlight = list(D21_C2 = rownames(d21_c2)), pt.size = 2, highlight_color = c("#00BFC4"))
DimPlot_scCustom(integrated, group.by = "dataset", colors_use = ColorBlind_Pal(), pt.size = 1)
DimPlot_scCustom(integrated, group.by = "atac_clones_2", colors_use = c("#F8766D", "#00BFC4"), pt.size = 1)
DimPlot_scCustom(integrated, group.by = "cell_types", pt.size = 1)
DimPlot_scCustom(integrated, group.by = "cat", pt.size = 1, colors_use = stepped)
grid::grid.newpage()
grid.table(data.frame(table(integrated$cat)))
VlnPlot(integrated, features = "AR_signature1", group.by = "cat", cols = stepped, pt.size = 0) + 
  scale_fill_manual(values = stepped, drop = FALSE)
VlnPlot_scCustom(integrated, features = "NE_signature1", group.by = "cat", colors_use = stepped, pt.size = 0) + 
  scale_fill_manual(values = stepped, drop = FALSE)
dev.off()


pdf("epiAneufinder_0528.pdf")

Meta_Highlight_Plot_man(integrated, meta_data_column = "category", pt.size = 0.7,
                        meta_data_highlight = c("D2 Clone 1", "D2 Clone 2"), c("#FB8072", "#80B1D3"),
    order = c( "D2 Clone 2", "D2 Clone 1"))
Meta_Highlight_Plot_man(integrated, meta_data_column = "category", pt.size = 0.7,
                        meta_data_highlight = c("D14 Clone 1", "D14 Clone 2"), c("#FB8072", "#80B1D3"),
                        order = c( "D14 Clone 2", "D14 Clone 1"))
Meta_Highlight_Plot_man(integrated, meta_data_column = "category", pt.size = 0.7,
                        meta_data_highlight = c("D21 Clone 1", "D21 Clone 2"), c("#80B1D3", "#FB8072"),
                        order = c( "D21 Clone 2", "D21 Clone 1"))

dev.off()



pdf("epiAneufinder_0524_signatures.pdf")
FeaturePlot_scCustom(integrated, order = T, pt.size = 0.7, features = "AR_signature1") + 
  labs(x = "atac_UMAP 1", y = "atac_UMAP 2", title = "AR Signature (RNA)")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 1.5), na.value = "lightgrey")
FeaturePlot(integrated, order = T, pt.size = 0.7, features = "NE_signature1") +
  labs(x = "atac_UMAP 1", y = "atac_UMAP 2", title = "NE Signature (RNA)")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 0.6), na.value = "lightgrey")

VlnPlot(integrated, features = "AR_signature1", group.by = "category", cols = stepped[c(1,2,5,6,9,10,12)], pt.size = 0) + 
  scale_fill_manual(values = stepped[c(1,2,5,6,9,10,12)], drop = FALSE)
VlnPlot_scCustom(integrated, features = "NE_signature1", group.by = "category", colors_use = stepped[c(1,2,5,6,9,10,12)], pt.size = 0) + 
  scale_fill_manual(values = stepped[c(1,2,5,6,9,10,12)], drop = FALSE)
dev.off()

