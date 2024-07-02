library(tidyverse)
library(epiAneufinder)
library(cowplot)
library(gridExtra)
library(patchwork)
library(scCustomize)

d2_res_table <- read.table("../d2/epiAneufinder_results/results_table.tsv") # 2154
d14_res_table <- read.table("../d14/epiAneufinder_results/results_table.tsv") # 2495 
d21_res_table <- read.table("../d21/epiAneufinder_results/results_table.tsv") # 4265

# ------
colnames(d2_res_table) <- colnames(d2_res_table) %>% str_remove("cell.") %>% str_replace_all("\\.", "-")
colnames(d14_res_table) <- colnames(d14_res_table) %>% str_remove("cell.") %>% str_replace_all("\\.", "-")
colnames(d21_res_table) <- colnames(d21_res_table) %>% str_remove("cell.") %>% str_replace_all("\\.", "-")

endogenous_integrated <- readRDS("../../../../multiome2/atac_integration/atac_integrated_new.RDS")
combined_cells <- Cells(endogenous_integrated) #7720
d2_cells <- combined_cells[str_detect(combined_cells, "d2(?!1)")] %>% str_remove("d2_")
d14_cells <- combined_cells[str_starts(combined_cells, "d14")] %>% str_remove("d14_")
d21_cells <- combined_cells[str_starts(combined_cells, "d21")] %>% str_remove("d21_")

d2_res_table_pass <- d2_res_table[, c(1:3, which(colnames(d2_res_table) %in% d2_cells))] #2113
d14_res_table_pass <- d14_res_table[, c(1:3, which(colnames(d14_res_table) %in% d14_cells))] # 1742
d21_res_table_pass <- d21_res_table[, c(1:3, which(colnames(d21_res_table) %in% d21_cells))] # 3835

colnames(d2_res_table_pass)[4:dim(d2_res_table_pass)[2]] <- paste0(colnames(d2_res_table_pass)[4:dim(d2_res_table_pass)[2]], "-d2")
colnames(d14_res_table_pass)[4:dim(d14_res_table_pass)[2]] <- paste0(colnames(d14_res_table_pass)[4:dim(d14_res_table_pass)[2]], "-d14")
colnames(d21_res_table_pass)[4:dim(d21_res_table_pass)[2]] <- paste0(colnames(d21_res_table_pass)[4:dim(d21_res_table_pass)[2]], "-d21")

d2_res_table_pass$coor <-  paste0(d2_res_table_pass$seq, "-", d2_res_table_pass$start, ":", d2_res_table_pass$end)
d14_res_table_pass$coor <-  paste0(d14_res_table_pass$seq, "-", d14_res_table_pass$start, ":", d14_res_table_pass$end)
d21_res_table_pass$coor <-  paste0(d21_res_table_pass$seq, "-", d21_res_table_pass$start, ":", d21_res_table_pass$end)

combined_res_table_1 <- inner_join(d2_res_table_pass[,4:dim(d2_res_table_pass)[2]], d14_res_table_pass[, 4:dim(d14_res_table_pass)[2]], by = "coor") # 3850
combined_res_table <- inner_join(combined_res_table_1, d21_res_table_pass, by = "coor") # 7685
combined_res_table <- combined_res_table %>% relocate(c("seq", "start", "end", "coor"))
combined_res_table <- combined_res_table %>% select(!coor)
colnames(combined_res_table)[4:length(combined_res_table)] <- paste0("cell-", colnames(combined_res_table[4:length(combined_res_table)]))  

saveRDS(combined_res_table, "combined_res_table.RDS")


