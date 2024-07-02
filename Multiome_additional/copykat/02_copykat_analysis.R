library(tidyverse)
library(Seurat)
library(scCustomize)

polychrome_pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
polychrome_pal <- polychrome_pal[c(1, 3, 6:36)]

integrated <- readRDS("../../rna_integration/rna_integrated.RDS")

d2_copykat <- read_delim("d2/D2_copykat_prediction.txt")
d14_copykat <- read_delim("d14/D14_copykat_prediction.txt")
d21_add_copykat <- read_delim("d21_add/D21_add_copykat_prediction.txt")

d2_seu <- readRDS("../../../multiome/single_d2/simple_d2.RDS")
d14_seu <- readRDS("../../../multiome/single/simple_d14.RDS")
d21_add_seu <- readRDS("../../d21/simple_d21.RDS")

# filtering copykat results to only cells we kept
d2_cells_status <- d2_copykat[which(d2_copykat$cell.names %in% Cells(d2_seu)),]
d2_status <- d2_cells_status$copykat.pred
names(d2_status) <- d2_cells_status$cell.names
d2_seu$copykat <- d2_cells_status$copykat.pred

DimPlot_scCustom(d2_seu, group.by = "copykat", pt.size = 1)

d14_cells_status <- d14_copykat[which(d14_copykat$cell.names %in% Cells(d14_seu)),]
d14_seu$copykat <- d14_cells_status$copykat.pred
DimPlot_scCustom(d14_seu, group.by = "copykat", pt.size = 1)

d21_cells_status <- d21_add_copykat[which(d21_add_copykat$cell.names %in% Cells(d21_add_seu)),]
d21_add_seu$copykat <- d21_cells_status$copykat.pred
DimPlot_scCustom(d21_add_seu, group.by = "copykat", pt.size = 1, colors_use = c(NavyAndOrange(), polychrome_pal[1]))

all_cells_barcodes <- c(paste0("d2-", d2_cells_status$cell.names), paste0("d14-", d14_cells_status$cell.names),
                      paste0("d21-", d21_cells_status$cell.names))
all_cells_status <- c(d2_cells_status$copykat.pred, d14_cells_status$copykat.pred, d21_cells_status$copykat.pred)
names(all_cells_status) <- all_cells_barcodes
integrated$copykat <- all_cells_status
DimPlot_scCustom(integrated, group.by = "copykat", colors_use = c(NavyAndOrange(), polychrome_pal[1]), pt.size = 1)

# copykat cna results
d2_cna <- read_delim("d2/D2_copykat_CNA_results.txt")
colnames(d2_cna) <- colnames(d2_cna) %>% str_replace("\\.", "-")
d2_cna <- d2_cna[, c(1:3, which(colnames(d2_cna) %in% d2_cells_status$cell.names))]
colnames(d2_cna)[4:ncol(d2_cna)] <- paste0("d2-", colnames(d2_cna)[4:ncol(d2_cna)])

d14_cna <- read_delim("d14/D14_copykat_CNA_results.txt")
colnames(d14_cna) <- colnames(d14_cna) %>% str_replace("\\.", "-")
d14_cna <- d14_cna[, c(1:3, which(colnames(d14_cna) %in% d14_cells_status$cell.names))]
colnames(d14_cna)[4:ncol(d14_cna)] <- paste0("d14-", colnames(d14_cna)[4:ncol(d14_cna)])

d21_cna <- read_delim("d21_add/D21_add_copykat_CNA_results.txt")
colnames(d21_cna) <- colnames(d21_cna) %>% str_replace("\\.", "-")
d21_cna <- d21_cna[, c(1:3, which(colnames(d21_cna) %in% d21_cells_status$cell.names))]
colnames(d21_cna)[4:ncol(d21_cna)] <- paste0("d21-", colnames(d21_cna)[4:ncol(d21_cna)])

c(colnames(d2_cna)[4:ncol(d2_cna)], colnames(d14_cna)[4:ncol(d14_cna)],colnames(d21_add_cna)[4:ncol(d21_add_cna)]) %>% length()

# combining cna matrix
combined_cna_temp <- inner_join(d2_cna, d14_cna, by = c("chrom", "chrompos", "abspos"))
colnames(combined_cna_temp)[str_detect(colnames(combined_cna_temp), pattern = "[xy]+")]

combined_cna <- inner_join(combined_cna_temp, d21_cna, by = c("chrom", "chrompos", "abspos"))
colnames(combined_cna)[str_detect(colnames(combined_cna), pattern = "[xy]+")]

write_delim(combined_cna, "combined_cna.txt")
combined_cna_barcode_only <- str_extract(colnames(combined_cna)[4:ncol(combined_cna)], "(?<=-)[:graph:]+")

combined_pred <-rbind(d2_cells_status, d14_cells_status, d21_cells_status)
combined_pred_status <- combined_pred[combined_pred$cell.names %in% combined_cna_barcode_only, 2]
write_delim(combined_pred_status, "combined_pred.txt")

# clones ------
tree_2 <- readRDS("heatmap/tree2.RDS")
names(tree_2) <- names(tree_2) %>% str_replace("-", "_")
names(tree_2) %in% Cells(integrated) %>% table()
integrated$tree_2 <- tree_2
DimPlot_scCustom(integrated, group.by = "tree_2", colors_use = color_clones, pt.size = 1)


