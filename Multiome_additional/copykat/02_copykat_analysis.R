library(tidyverse)

source("../../../../multiome2/atac_integration/epi_loose/meta_highlight_manual.R")
source("../../../../multiome2/atac_integration/epi_loose/Meta_Highlight_Plot_manual.R")

d2_cna <- read_delim("d2/D2_copykat_CNA_results.txt")
colnames(d2_cna) <- c(colnames(d2_cna)[1:3], paste0("d2-", colnames(d2_cna)[4:ncol(d2_cna)], "-1"))

d14_cna <- read_delim("d14/D14_copykat_CNA_results.txt")
colnames(d14_cna) <- c(colnames(d14_cna)[1:3], paste0("d14-", colnames(d14_cna)[4:ncol(d14_cna)], "-1"))

d21_cna <- read_delim("d21/D21_add_copykat_CNA_results.txt")
colnames(d21_cna) <- c(colnames(d21_cna)[1:3], paste0("d21-", colnames(d21_cna)[4:ncol(d21_cna)], "-1"))

d2_cna_mat <- d2_cna[4:ncol(d2_cna)]
d14_cna_mat <- d14_cna[4:ncol(d14_cna)]
d21_cna_mat <- d21_cna[4:ncol(d21_cna)]

# filtering on heatmaps ----

# Define the filtering threshold and range
threshold <- 0.8  # 80%
lower_bound <- -0.1
upper_bound <- 0.1

# Calculate the proportion of values within the range for each row
d21_proportion_in_range <- apply(d21_cna_mat, 2, function(row) {
  mean(row >= lower_bound & row <= upper_bound)
})

cna_based_filtered_d21 <- d21_cna_mat[, d21_proportion_in_range < threshold]
cna_based_filtered_d21_fin <- cbind(d21_cna[,1:3], cna_based_filtered_d21)

saveRDS(cna_based_filtered_d21_fin, "filtered_matrix_fin.RDS")

d2_proportion_in_range <- apply(d2_cna_mat, 2, function(row) {
  mean(row >= lower_bound & row <= upper_bound)
})

cna_based_filtered_d2 <- d2_cna_mat[, d2_proportion_in_range < threshold]
cna_based_filtered_d2_fin <- cbind(d2_cna[,1:3], cna_based_filtered_d2)
saveRDS(cna_based_filtered_d2_fin, "cna_based_filtered_d2_fin.RDS")

d14_proportion_in_range <- apply(d14_cna_mat, 2, function(row) {
  mean(row >= lower_bound & row <= upper_bound)
})

cna_based_filtered_d14 <- d14_cna_mat[, d14_proportion_in_range < threshold]
cna_based_filtered_d14_fin <- cbind(d14_cna[,1:3], cna_based_filtered_d14)
saveRDS(cna_based_filtered_d14_fin, "cna_based_filtered_d14_fin.RDS")


combined_cna_temp_1_heatmap <- inner_join(cna_based_filtered_d2_fin, cna_based_filtered_d14_fin, by = c("chrom", "chrompos", "abspos"))
combined_cna_1_heatmap <- inner_join(combined_cna_temp_1_heatmap, cna_based_filtered_d21_fin, by = c("chrom", "chrompos", "abspos"))
write_delim(combined_cna_1_heatmap, "combined_cna_1_heatmap.txt")

