library(epiAneufinder)
library(tidyverse)
library(data.table)
library(cowplot)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
print(as.numeric(args[1]))
sink("color/R output.txt")
combined_res_table <- readRDS("combined_res_table.RDS")
source("../../plot_karyo_annotated_manual_color.R")

find_plot_subclones <- function(res_table, depth = 2) {
  
  print("in function")
  subclones <- split_subclones(res_table, tree_depth=depth,
                               plot_tree=TRUE,
                               plot_path=paste0("color/subclones_", depth, ".pdf"),
                               plot_width=4,
                               plot_height=3)
  
  saveRDS(subclones, paste0("color/subclones_", depth, ".RDS"))
subclones <- readRDS("subclones_2.RDS")

  subclones$day <- subclones$cell %>% str_extract("(?<=-)d\\d+")
  
  annot_dt <-data.frame(cell=subclones$cell,
                        day = subclones$day,
                        annot= paste0("Clone ",subclones$subclone),
                        day_color = factor(subclones$day, levels = c("d2", "d14", "d21")))
  colnames(res_table) <- gsub("\\.","-",colnames(res_table))
  
  
  karyo <- plot_karyo_manual(res_table=res_table,
                             plot_path=paste0("color/karyo_annotated_", depth, ".png"),
                             annot_dt=annot_dt)
  res_list <- list(subclones, annot_dt, karyo)
  return(res_list)
  print("done karyo")
}

depth <- find_plot_subclones(res_table = combined_res_table, depth = as.numeric(args[1]))
saveRDS(depth, paste0("color/depth", args[1], ".RDS"))
sink()
