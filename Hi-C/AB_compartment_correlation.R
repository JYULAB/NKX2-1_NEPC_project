library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# read in files ------

lucap_100kb_list_files <- list.files("../results_100kb", "bedGraph", full.names = TRUE)
timecourse_100kb_list_files <- list.files("/projects/p20023/Viriya/analysis/foxa2/hi-c/cooltools/ab", "r.bedGraph",full.names = TRUE)

all_100kb_files <- c(lucap_100kb_list_files, timecourse_100kb_list_files)
eigenvector_100kb_list <- lapply(all_100kb_files, read_delim, col_names = c("chrom", "start", "end", "score"))

samples <- c("LuCaP105CR", "LuCaP145.2", "LuCaP145.1", "LuCaP147",
             "LuCaP147CR", "LuCaP167CR", "LuCaP173.1", "LuCaP35CR", 
                 "LuCaP70CR", "LuCaP77CR", "LuCaP86.2CR", "LuCaP93", "NCI-H660",
                 "D0", "D14", "D21", "D28")
names(eigenvector_100kb_list) <- samples

list_eigenvectors_100kb <- lapply(seq_along(eigenvector_100kb_list), function(x) {
  colnames(eigenvector_100kb_list[[x]]) <- c("chrom", "start", "end", names(eigenvector_100kb_list)[[x]])
  eigenvector_100kb_list[[x]]
})

# Fig 2C, all PDX and FOXA2OE LNCaP
df_eigenvectors <- purrr::reduce(list_eigenvectors_100kb, full_join, by = c("chrom", "start", "end")) %>% drop_na()
# Fig S1A, 7PDXs, NCI-H660 and LNCaP
df_eigenvectors_9 <- df_eigenvectors[, c(1:3, 5:7, 11:13, 15:17)]
# Only NEPC PDXs and NCI-H660
df_eigenvectors_pdx_nepc <- df_eigenvectors[, c(1:3, 5:6, 10, 15, 16)]


# correlations calculations ---
calc_corr_mat <- function(df) {
  variances <- apply(df[,4:ncol(df)], 1, var)
  
  # Calculate the threshold for the top 10%
  top_10_threshold <- quantile(variances, probs = 0.9)
  top_25_threshold <- quantile(variances, probs = 0.75)
  
  # Filter rows that have variances above the threshold
  top_10_percent <- df[variances >= top_10_threshold, ]
  top_10_mat <- top_10_percent[,4:ncol(top_10_percent)] %>% as.matrix()
  
  top_25_percent <- df[variances >= top_25_threshold, ]
  top_25_mat <- top_25_percent[,4:ncol(top_25_percent)] %>% as.matrix()
  
  top_10_corr_mat <- cor(top_10_mat)
  top_25_corr_mat <- cor(top_25_mat)
  full_corr_mat <- cor(as.matrix(df[,4:ncol(df)]))
  
  l <- list(t10 = top_10_corr_mat, t25 = top_25_corr_mat, full = full_corr_mat)
}

corr <- calc_corr_mat(df_eigenvectors)
corr_pdx_only <- calc_corr_mat(df_eigenvectors_pdx_only)
corr_pdx_nepc <- calc_corr_mat(df_eigenvectors_pdx_nepc)
corr_no_nci <- calc_corr_mat(df_eigenvectors_no_nci)
corr_9 <- calc_corr_mat(df_eigenvectors_9)

# heatmap colors -----
color_fun <- colorRamp2(c(-0.5,0,1), c("blue", "white", "red"))
type_17 <- c("CRPC", "NEPC", "NEPC", "CRPC", 
             "CRPC", "CRPC", "NEPC", "CRPC", 
             "CRPC", "CRPC", "CRPC", "NEPC", "NEPC", rep("LnCaP", 4))
type_colors_17 <- case_when(type_17 == "NEPC" ~ "red",
                            type_17 == "CRPC" ~ "blue",
                            type_17 == "LnCaP" ~ "purple")

hm_anno_17 <- HeatmapAnnotation(samples = anno_text(colnames(corr[[1]]), 
                                                    rot = 30, which = "column", 
                                                    gp = gpar(col = type_colors_17)))

type_colors_16 <- type_colors_17[-13]

hm_anno_16 <- HeatmapAnnotation(samples = anno_text(colnames(corr_no_nci[[1]]), 
                                                    rot = 30, which = "column", 
                                                    gp = gpar(col = type_colors_16)))

type_colors_13 <- case_when(type_17 == "NEPC" ~ "red",
                            type_17 == "CRPC" ~ "blue")
hm_anno_13 <- HeatmapAnnotation(samples = anno_text(colnames(corr_pdx_only[[1]]), 
                                                    rot = 30, which = "column", 
                                                    gp = gpar(col = type_colors_13)))

type_9 <- c("NEPC", "NEPC", "CRPC", "CRPC", "CRPC", "CRPC", "NEPC", "NEPC", rep("LnCaP", 1))
type_colors_9 <- case_when(type_9 == "NEPC" ~ "red",
                           type_9 == "CRPC" ~ "blue",
                           type_9 == "LnCaP" ~ "purple")

hm_anno_9 <- HeatmapAnnotation(samples = anno_text(colnames(corr_9[[1]]), 
                                                   rot = 30, which = "column", 
                                                   gp = gpar(col = type_colors_9)))

type_colors_5 <- c("red", "red", "brown", "brown", "red")
text_5 <- paste0(colnames(corr_pdx_nepc[[1]]), c("\nF+/N+", "\nF+/N+", "\nF-/N+", "\nF-/N+","\nF+/N+"))
hm_anno_5 <- HeatmapAnnotation(samples = anno_text(text_5, 
                                                    rot = 30, which = "column", 
                                                    gp = gpar(col = type_colors_5)))
hm_anno_5_row <- rowAnnotation(samples = anno_text(text_5, 
                                               rot = 30, which = "column", 
                                               gp = gpar(col = type_colors_5)))

print("okay")

# All samples heatmap example -----

corr_hm_25 <- Heatmap(corr[[2]],
                          name = "PC1 correlation",
                          col = color_fun,
                          cluster_rows = TRUE,
                          row_names_gp = gpar(col = type_colors_17),
                          show_row_dend = TRUE,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(sprintf("%.1f", corr[[2]][i, j]), x, y, gp = gpar(fontsize = 10))
                          },
                          
                          
                          column_title = "Top 25% Most Variable Eigenvector Bins at 100kb",
                          column_title_side = "bottom",
                          cluster_columns = TRUE,
                          show_column_dend = FALSE,
                          show_column_names = FALSE,
                          column_names_gp = gpar(col = type_colors_17),
                          bottom_annotation = hm_anno_17,
                          
                          heatmap_legend_param = list(direction = "horizontal"))
draw(corr_hm_25, heatmap_legend_side = "top")

# automated heatmaps -----

make_heatmap <- function(corr_list, which = NULL, col = type_colors_17, anno = hm_anno_17, row_anno = NULL) {
  
  if(names(corr_list)[which] == "full") {
    title <- paste0("All Eigenvector Bins at 100kb")
  } else {
    percent <- names(corr_list)[which] %>% str_remove("t")
    title <-  paste0("Top ", percent ,"% Most Variable Eigenvector Bins at 100kb")
  }
  
  if(is.null(row_anno) == TRUE) {
    row_cond <- TRUE
  } else {
    row_cond <- FALSE
  }
  Heatmap(corr_list[[which]],
          name = "PC1 correlation",
          col = color_fun,
          cluster_rows = TRUE,
          row_names_gp = gpar(col = col),
          show_row_names = row_cond,
          right_annotation = row_anno,
          show_row_dend = TRUE,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.1f", corr_list[[which]][i, j]), x, y, gp = gpar(fontsize = 10))
          },


          column_title = title,
          column_title_side = "bottom",
          cluster_columns = TRUE,
          show_column_dend = FALSE,
          show_column_names = FALSE,
          column_names_gp = gpar(col = col),
          bottom_annotation = anno,

          heatmap_legend_param = list(direction = "horizontal"))
}

hm_corr_25 <- make_heatmap(corr, which = 2)
draw(hm_corr_25)

hm_pdx_nepc_25 <- make_heatmap(corr_pdx_nepc, which = 2, type_colors_5, hm_anno_5,hm_anno_5_row)
draw(hm_pdx_nepc_25, heatmap_legend_side = "top")

hm_9 <- make_heatmap(corr_9, which = 2, type_colors_9, hm_anno_9)
draw(hm_9, heatmap_legend_side = "top")

