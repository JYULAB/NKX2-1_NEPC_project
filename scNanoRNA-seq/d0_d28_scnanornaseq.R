library(Seurat)
library(tidyverse)
library(scCustomize)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
set.seed(1234)

polychrome_pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
polychrome_pal <- polychrome_pal[c(1, 3, 6:36)]
stepped <- DiscretePalette_scCustomize(num_colors = 24, palette = "stepped")
stepped <- stepped[c(7,6,5, 16,15,14,  12,11,10,9, 24, 24)]
source("../../multiome2/atac_integration/epi_loose/meta_highlight_manual.R")
source("../../multiome2/atac_integration/epi_loose/Meta_Highlight_Plot_manual.R")

# loading in scNanoGPX data ---------
mtx_gene_d28 <- read_tsv("../LnCap_1month/matrix_gene.tsv", comment = "#")
mtx_gene_unique_genes_d28 <- distinct(mtx_gene_d28, gene_name, .keep_all = TRUE)
gene_matrix_d28 <- mtx_gene_unique_genes_d28[,8:dim(mtx_gene_unique_genes_d28)[2]]

# D0
mtx_gene_d0 <- read_tsv("../LnCap_parental/matrix_gene.tsv", comment = "#")
mtx_gene_unique_genes_d0 <- distinct(mtx_gene_d0, gene_name, .keep_all = TRUE)
gene_matrix_d0 <- mtx_gene_unique_genes_d0[,8:dim(mtx_gene_unique_genes_d0)[2]]

# make seurat -----
d28_matrix_gene <- as.sparse(gene_matrix_d28, row.names = mtx_gene_unique_genes_d28$gene_name)
d0_matrix_gene <- as.sparse(gene_matrix_d0, row.names = mtx_gene_unique_genes_d0$gene_name)

# D28
d28 <- CreateSeuratObject(d28_matrix_gene, project = "D28")
d28[["percent.mt"]] <- PercentageFeatureSet(d28, pattern = "^MT-")

d28_v <- VlnPlot(d28, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 ncol = 3, pt.size = 0, combine = FALSE)
d28_cutoff <- c(2300, 6000, 10)
d28_v1 <- lapply(seq_along(d28_v), function(x) d28_v[[x]] + geom_hline(yintercept = d28_cutoff[[x]]) +
                   theme(legend.position="none"))
d28_v1[[1]] + d28_v1[[2]] + d28_v1[[3]] + plot_annotation("D28")

d28 <- subset(d28, subset = percent.mt < 10 & nFeature_RNA > 2300 & nCount_RNA)


# D0
d0 <- CreateSeuratObject(d0_matrix_gene, project = "D0")
d0[["percent.mt"]] <- PercentageFeatureSet(d0, pattern = "^MT-")

d0_v <- VlnPlot(d0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, pt.size = 0, combine = FALSE)
d0_cutoff <- c(3000, 10000, 10)
d0_v1 <- lapply(seq_along(d0_v), function(x) d0_v[[x]] + geom_hline(yintercept = d0_cutoff[[x]]) +
                  theme(legend.position="none")) 
d0_v1[[1]] + d0_v1[[2]] + d0_v1[[3]] + plot_annotation("D0")

d0 <- subset(d0, subset = percent.mt < 10 & nFeature_RNA > 2300 & nCount_RNA)
all_genes <- rownames(d28)

process <- function(obj) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj, features = all_genes)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj))
  obj <- RunUMAP(obj, dims = 1:30)
}

d28 <- process(d28)
d14 <- process(d14)
d0 <- process(d0)
saveRDS(d0, "seurat_d0.RDS")
saveRDS(d14, "seurat_d14.RDS")
saveRDS(d28, "seurat_d28.RDS")

# d0 <- readRDS("../analysis_clean/seurat_d0.RDS")
# d14 <- readRDS("../analysis_clean/seurat_d14.RDS")
# d28 <- readRDS("../analysis_clean/seurat_d28.RDS")

# mutation annotation ------
anno_d28 <- read_tsv("../LnCap_1month/annovar.hg38_multianno.tsv")
anno_d0 <- read_tsv("../LnCap_parental/annovar.hg38_multianno.tsv")

anno_d28 <- anno_d28 %>% mutate(chrompos = paste0(CHROM, "_", POS))
anno_d0 <- anno_d0 %>% mutate(chrompos = paste0(CHROM, "_", POS))

# selecting for high confindece mutations ------
mtx_mutations_d28 <- read_tsv("../LnCap_1month/matrix_SNV_dp.tsv.gz")
mtx_mutations_d0 <- read_tsv("../LnCap_parental/matrix_SNV_dp.tsv")

is_mutated <- function(mutation) {
  rd <- lapply(mutation, function(x) str_split_1(x, "/")) # get read depth per cell in mutation
  mut <- sapply(rd, function(allele) {
    ifelse(any(as.numeric(allele[2:length(allele)]) >= 2), TRUE, FALSE) # mutated if any alternate allele has 2 supporting reads
  })
  return(mut)
}

# only the matrix portion, excluing metadata, only including high quality cells 
mat_d0 <- as.matrix(mtx_mutations_d0[,colnames(d0)])
mat_d28 <- as.matrix(mtx_mutations_d28[,colnames(d28)])

mutation_status_d0 <- apply(mat_d0, MARGIN = 1, function(mutation) is_mutated(mutation)) %>% t()
mutation_status_d28 <- apply(mat_d28, MARGIN = 1, function(mutation) is_mutated(mutation)) %>% t()

saveRDS(mutation_status_d0, "mut_d0.RDS")
saveRDS(mutation_status_d28, "mut_d28.RDS")

mutation_status_d0 <- readRDS("../analysis_clean/mut_d0.RDS")
mutation_status_d28 <- readRDS("../analysis_clean/mut_d28.RDS")

# keep mutation if it is mutated in more than 5% of cells in each condition
keep_d0 <- ifelse(rowSums(mutation_status_d0) >= dim(mutation_status_d0)[2] *0.05, TRUE, FALSE)
keep_d28 <- ifelse(rowSums(mutation_status_d28) >= dim(mutation_status_d28)[2] *0.05, TRUE, FALSE)


filtered_mutations <- data_frame(Condition = c(rep("D0", sum(keep_d0)),
                                               rep("D28", sum(keep_d28))),
                                 Type = c(anno_d0$Func_type[keep_d0],
                                          anno_d28$Func_type[keep_d28]))

ggplot(filtered_mutations, aes(x = Condition)) +
  geom_bar(aes(fill = Type)) +
  scale_fill_manual(values = polychrome_pal[1:14])


filtered_anno_d0 <- anno_d0[keep_d0,] %>% mutate(chrompos = paste0(CHROM, "_", POS))
filtered_anno_d28 <- anno_d28[keep_d28,] %>% mutate(chrompos = paste0(CHROM, "_", POS))

library(Vennerable)
venn_filtered <- Venn(Sets = list(as.character(filtered_anno_d0$chrompos), as.character(filtered_anno_d28$chrompos)),
                      SetNames = c("D0", "D28"))
plot(venn_filtered)


# binary mutation status ---------
filtered_mat_d0 <- mat_d0[keep_d0,]
filtered_mat_d28 <- mat_d28[keep_d28,]

# heatmap param ------
col_anno_day_df <- data.frame(timepoint = c(rep("D0", 1858), rep("D28", 1729)))

col_anno <- HeatmapAnnotation(df = col_anno_day_df,
                              col = list(timepoint = c("D0" = "#FFA500", "D28" = "#009E73")))
col_anno_polychrome <- HeatmapAnnotation(df = col_anno_day_df,
                                         col = list(timepoint = c("D0" = "#5A5156FF", "D28" = "#3283FEFF")))
col_fun_1 <- circlize::colorRamp2(c(0, 1), c("white", "red"))
# 1093 shared
mut_d28_shared <- filtered_anno_d28$chrompos %in% filtered_anno_d0$chrompos
mut_shared <- filtered_anno_d28[mut_d28_shared,]
dim(mut_shared)

# VAF calculation ----
calc_vaf <- function(mutation) {
  rd <- lapply(mutation, function(x) str_split_1(x, "/")) # get read depth per cell in mutation
  
  # mutation coverage / all coverage
  cov_freq <- sapply(rd, function(allele) {
    cov_freq <- sum(as.numeric(allele[2:length(allele)])) / sum(as.numeric(allele))
  })
  
  
}

vaf_d0 <- apply(filtered_mat_d0, MARGIN = 1, function(mutation) calc_vaf(mutation)) %>% t()
vaf_d28 <- apply(filtered_mat_d28, MARGIN = 1, function(mutation) calc_vaf(mutation)) %>% t()

vaf_d0_df <- vaf_d0 %>% as.data.frame()
colnames(vaf_d0_df) <- paste0("D0-", colnames(vaf_d0_df))
vaf_d0_df$chrompos <- filtered_anno_d0$chrompos


vaf_d28_df <- vaf_d28 %>% as.data.frame()
colnames(vaf_d28_df) <- paste0("D28-", colnames(vaf_d28_df))
vaf_d28_df$chrompos <- filtered_anno_d28$chrompos

vaf_all <- full_join(vaf_d0_df, vaf_d28_df, by = "chrompos")
vaf_all <- vaf_all %>% relocate("chrompos")

# heatmap ----
rownames(vaf_all) <- vaf_all$chrompos
vaf_shared <- vaf_all[mut_shared$chrompos,2:dim(vaf_all)[2]] %>% as.matrix()
dim(vaf_shared)
hist(vaf_shared)

d28_loc <- colnames(vaf_shared) %in% paste0("D28-", colnames(d28))
splitting_rules <-  case_when(rowMeans(vaf_shared, na.rm = TRUE) >= 0.9 ~ "Mut",
                              rowMeans(vaf_shared, na.rm = TRUE) <= 0.3 ~ "WT",
                              rowMeans(vaf_shared[,d28_loc], na.rm = TRUE) >= 0.9 ~ "D28 Mut",
                              rowMeans(vaf_shared[,!d28_loc], na.rm = TRUE) >= 0.9 ~ "D28 WT",
                              .default = "Mixed")

hm_vaf_shared_s_row <- Heatmap(vaf_shared,
                               name = "VAF",
                               
                               cluster_rows = TRUE,
                               show_row_dend = FALSE,
                               show_row_names = FALSE,
                               row_split = splitting_rules,
                               row_title_rot = 0,
                               
                               cluster_columns = FALSE,
                               show_column_dend = FALSE,
                               show_column_names = FALSE,
                               
                               
                               show_heatmap_legend = TRUE,
                               top_annotation = col_anno_polychome
)

hmap_vaf_shared_s_row <- draw(hm_vaf_shared_s_row, merge_legend = TRUE)
draw(hmap_vaf_shared_s_row, merge_legend = TRUE)

# small number imputation ------
vaf_shared_coverage_imputed <- apply(vaf_shared, MARGIN = 2, function(x) ifelse(is.na(x), -0.1, x))
vaf_shared_coverage_imputed_kmeans <- kmeans(t(vaf_shared_coverage_imputed), 2)
clones_df <- data.frame(clone = paste0("Clone ", vaf_shared_coverage_imputed_kmeans$cluster),
                                 timepoint = str_extract(names(vaf_shared_coverage_imputed_kmeans$cluster), "D\\d+"))
table(clones_df$clone, clones_df$timepoint)
rownames(clones_df) <- names(vaf_shared_coverage_imputed_kmeans$cluster)
col_anno_kmeans <- HeatmapAnnotation(df = clones_df,
                                    col = list(timepoint = c("D0" = "#5A5156FF", "D28" = "#3283FEFF"),
                                               clone = c("Clone 1" = "orange", "Clone 2" = "navy")))

d28_loc_2 <- colnames(vaf_shared_962) %in% paste0("D28-", colnames(d28))
splitting_rules_2 <-  case_when(rowMeans(vaf_shared_962, na.rm = TRUE) >= 0.9 ~ "Mut",
                              rowMeans(vaf_shared_962, na.rm = TRUE) <= 0.3 ~ "WT",
                              rowMeans(vaf_shared_962[,d28_loc_2], na.rm = TRUE) >= 0.9 ~ "D28 Mut",
                              rowMeans(vaf_shared_962[,!d28_loc_2], na.rm = TRUE) >= 0.9 ~ "D28 WT",
                              .default = "Mixed")


hm_vaf_shared_s_row_c <- Heatmap(vaf_shared,
                               name = "VAF",
                               
                               # col = col_fun,
                               cluster_rows = FALSE,
                               show_row_dend = FALSE,
                               show_row_names = FALSE,
                               row_split = splitting_rules,
                               row_title_rot = 0,
                               
                               cluster_columns = FALSE,
                               column_split = clones_df$clone,
                               cluster_column_slices = FALSE,
                               column_gap = unit(0, "mm"),
                               column_title = "1093 mutations",
                               show_column_dend = FALSE,
                               show_column_names = FALSE,
                               
                               show_heatmap_legend = TRUE,
                               top_annotation = col_anno_kmeans
)


hmap_vaf_shared_s_row_c <- draw(hm_vaf_shared_s_row_c, merge_legend = TRUE)

hm_vaf_shared_s_row_imputed <- Heatmap(vaf_imputed_shared,
                                       name = "VAF",
                                       
                                       cluster_rows = TRUE,
                                       show_row_dend = FALSE,
                                       show_row_names = FALSE,
                                       row_split = splitting_rules,
                                       row_title_rot = 0,
                                       
                                       cluster_columns = FALSE,
                                       column_split = clones_df$clone,
                                       show_column_dend = FALSE,
                                       show_column_names = FALSE,
                                       
                                       show_heatmap_legend = TRUE,
                                       top_annotation = col_anno_kmeans
)
hmap_vaf_shared_s_row_imputed <- draw(hm_vaf_shared_s_row_imputed, merge_legend = TRUE)

# germline ------
# if a mutation has any alt alelle with 1 read, call it mutated for that cell

is_mutated_germline_fast <- function(mutation) {
  rd <- strsplit(mutation, "/") # get read depth per cell in mutation
  mut <- vapply(rd, function(allele) {
    # Ensure the first allele is 0 and any alternate allele (from the second onward) has ≥1 supporting read
    if (as.numeric(allele[1]) != 0) {
      return(FALSE)
    }
    any(as.numeric(allele[-1]) >= 1)
  }, logical(1)) # vapply is faster and safer than sapply
  return(mut)
}

mutation_status_1_read_d0_1 <- apply(mat_d0, MARGIN = 1, function(mutation) is_mutated_germline_fast(mutation)) %>% t()
mutation_status_1_read_d28_1 <- apply(mat_d28, MARGIN = 1, function(mutation) is_mutated_germline_fast(mutation)) %>% t()

# if a mutation has all zeros, mark as not expressed
is_expressed <- function(mutation) {
  ifelse(mutation %in% c("0/0", "0/0/0", "0/0/0/0"), FALSE, TRUE) # if 0/0, mark as not expressed
}

exp_mut_d0_1 <- apply(mat_d0, MARGIN = 1, function(mutation) is_expressed(mutation)) %>% t()
exp_mut_d28_1 <- apply(mat_d28, MARGIN = 1, function(mutation) is_expressed(mutation)) %>% t()

exp_mut_d0_hq <- exp_mut_d0[,colnames(d0)]
exp_mut_d28_hq <- exp_mut_d28[,colnames(d28)]

keep_d0_germline <- rowSums(mutation_status_1_read_d0_1) >= rowSums(exp_mut_d0_1) *0.95 & 
  rowSums(mat_d0 == "0/0") / ncol(mat_d0) <= 0.5
keep_d28_germline <- rowSums(mutation_status_1_read_d28_1) >= rowSums(exp_mut_d28_1) *0.95 & 
                              rowSums(mat_d28 == "0/0") / ncol(mat_d28) <= 0.5

filtered_anno_d0_germline <- anno_d0[keep_d0_germline,] %>% mutate(chrompos = paste0(CHROM, "_", POS))
filtered_anno_d28_germline <- anno_d28[keep_d28_germline,] %>% mutate(chrompos = paste0(CHROM, "_", POS))

venn_filtered_germline <- Venn(Sets = list(as.character(filtered_anno_d0_germline$chrompos),
                                           as.character(filtered_anno_d28_germline$chrompos)),
                      SetNames = c("D0","D28"))
plot(venn_filtered_germline)

mut_d28_shared_germline <- filtered_anno_d28_germline$chrompos %in% filtered_anno_d0_germline$chrompos
mut_shared_germline <- filtered_anno_d28_germline[mut_d28_shared_germline,]
# only the 139 shared germlines included
vaf_germline <- vaf_all[mut_shared_germline$chrompos,2:dim(vaf_all)[2]] %>% as.matrix()

hm_vaf_germline_s_row_manual_kmeans <- Heatmap(vaf_germline,
                                               name = "VAF",
                                               
                                               cluster_rows = FALSE,
                                               show_row_dend = FALSE,
                                               show_row_names = FALSE,
                                               row_title_rot = 0,
                                               
                                               cluster_columns = FALSE,
                                               show_column_dend = FALSE,
                                               show_column_names = FALSE,
                                               column_title = "139 shared germline mutations",
                                               column_split = clones_df$clone,
                                               column_gap = unit(0, "mm"),
                                               
                                               show_heatmap_legend = TRUE,
                                               top_annotation = col_anno_kmeans
)
hmap_vaf_germline_s_row_vaf_manual_kmeans <- draw(hm_vaf_germline_s_row_manual_kmeans, merge_legend = TRUE)


hm_vaf_germline_s_row_manual_kmeans_no_clones <- Heatmap(vaf_germline,
                                               name = "VAF",
                                               col = col_fun_1,
                                               
                                               cluster_rows = FALSE,
                                               show_row_dend = FALSE,
                                               show_row_names = FALSE,
                                               row_title_rot = 0,
                                               
                                               cluster_columns = FALSE,
                                               show_column_dend = FALSE,
                                               show_column_names = FALSE,
                                               column_title = "139 shared germline mutations",
                                               
                                               show_heatmap_legend = TRUE,
                                               top_annotation = col_anno_polychrome
)
hmap_vaf_germline_s_row_vaf_manual_kmeans_no_clones <- draw(hm_vaf_germline_s_row_manual_kmeans_no_clones, merge_legend = TRUE)

# including all individual day germline mutations
vaf_germline_all <- vaf_all[wt_status_all$chrompos,2:dim(vaf_all)[2]] %>% as.matrix()

row_split_def <- data.frame(shared = ifelse(wt_status_all$chrompos %in% mut_shared_germline$chrompos, "shared", "not shared"))
rownames(row_split_def) <- wt_status_all$chrompos
row_anno <- rowAnnotation(df = row_split_def,
                              col = list(shared = c("shared" = "green", "not shared" = "black")))
row_split_order <- c(venn_filtered_germline@IntersectionSets$`01`,
                     venn_filtered_germline@IntersectionSets$`10`,
                     venn_filtered_germline@IntersectionSets$`11`)

# row_order(hmap_binary_c) taken from heatmap below, which clusters rows by d28 only, d0 only and shared binary mutations status
hm_vaf_all_germline_s_row_manual_kmeans_no_clones <- Heatmap(vaf_germline_all,
                                                         name = "VAF",
                                                         # col = col_fun_1,
                                                         
                                                         cluster_rows = FALSE,
                                                         show_row_dend = FALSE,
                                                         show_row_names = FALSE,
                                                         row_title_rot = 0,
                                                         row_order = row_order(hmap_binary_c),
                                                         
                                                         cluster_columns = FALSE,
                                                         show_column_dend = FALSE,
                                                         show_column_names = FALSE,
                                                         column_title = "All found germline mutations",
                                                         
                                                         show_heatmap_legend = TRUE,
                                                         top_annotation = col_anno_polychrome
)
hmap_vaf_germline_s_row_vaf_manual_kmeans_no_clones <- draw(
  hm_vaf_all_germline_s_row_manual_kmeans_no_clones, merge_legend = TRUE)


# binary ------

count_wt_mut_cells <- function(mutation) {
  # Split the strings into components
  rd <- str_split(mutation, "/", simplify = TRUE)
  
  # Convert to numeric for faster operations
  rd <- apply(rd, 2, as.numeric)
  # Preallocate result vector
  mut <- character(nrow(rd))
  
  # WT condition: Reference allele is non-zero, others are zero
  is_wt <- rd[, 1] != 0 & rowSums(rd[, -1, drop = FALSE] != 0) == 0
  mut[is_wt] <- "WT"
  
  # MUT condition: At least one non-reference allele is non-zero
  is_mut <- rowSums(rd != 0) > 0 & !is_wt
  mut[is_mut] <- "MUT"
  
  # NA condition: All alleles are zero
  mut[!is_wt & !is_mut] <- NA
  
  return(mut)
}

wt_status_d0 <- apply(mat_d0, MARGIN = 1, function(mutation) count_wt_mut_cells(mutation)) %>% t()
wt_status_d28 <- apply(mat_d28, MARGIN = 1, function(mutation) count_wt_mut_cells(mutation)) %>% t()

wt_status_d0_df <- wt_status_d0 %>% as.data.frame()
wt_status_d0_df$chrompos <- anno_d0$chrompos

wt_status_d28_df <- wt_status_d28 %>% as.data.frame()
wt_status_d28_df$chrompos <- anno_d28$chrompos

# this version did not include the values for positions that's not shared by both timelines
wt_status_shared <- full_join(wt_status_d0_df[wt_status_d0_df$chrompos %in% filtered_anno_d0_germline$chrompos,], 
                           wt_status_d28_df[wt_status_d28_df$chrompos %in% filtered_anno_d28_germline$chrompos,],
                           by = "chrompos")
rownames(wt_status_shared) <- wt_status_shared$chrompos
wt_status_shared <- wt_status_shared %>% relocate("chrompos")
all(rownames(wt_status_shared$chrompos) == rownames(row_split_def))

wt_status_all <- full_join(wt_status_d0_df[wt_status_d0_df$chrompos %in% 
                                               c(filtered_anno_d0_germline$chrompos, filtered_anno_d28_germline$chrompos),], 
                           wt_status_d28_df[wt_status_d28_df$chrompos %in% 
                                              c(filtered_anno_d0_germline$chrompos, filtered_anno_d28_germline$chrompos),],
                           by = "chrompos")

rownames(wt_status_all) <- wt_status_all$chrompos
wt_status_all <- wt_status_all %>% relocate("chrompos")
all(rownames(wt_status_all$chrompos) == rownames(row_split_def))

wt_status_shared_mat <- as.matrix(wt_status_shared[, c(2:ncol(wt_status_shared))])
wt_status_all_mat <- as.matrix(wt_status_all[, c(2:ncol(wt_status_all))])

wt_status_shared_mat_encoded <- apply(wt_status_shared_mat, 2, function(x) case_when(
  x == "WT" ~ 0,
  x == "MUT" ~ 1,
  is.na(x) ~ 2
))

wt_status_all_mat_encoded <- apply(wt_status_all_mat, 2, function(x) case_when(
  x == "WT" ~ 0,
  x == "MUT" ~ 1,
  is.na(x) ~ 2
))
rownames(wt_status_all_mat_encoded) <- wt_status_all$chrompos
wt_status_all_mat_encoded_reordered <- wt_status_all_mat_encoded[wt_status_shared$chrompos,]

# heatmaps  -----
hm_binary_c <- Heatmap(wt_status_shared_mat_encoded,
                       name = "mutated",
                       
                       col = list("1" = "brown3", "0" = "blue", "2" = "#FFFFFF"),
                       heatmap_legend_param = list(border = "black",
                                                   labels = c("Mut", "WT", "NA")),
                       cluster_rows = TRUE,
                       show_row_dend = FALSE,
                       show_row_names = FALSE,
                       right_annotation = row_anno,
                       
                       show_column_names = FALSE,
                       show_column_dend = FALSE,
                       cluster_columns = FALSE,
                       show_heatmap_legend = TRUE,
                       top_annotation = col_anno_polychrome
)

hmap_binary_c <- draw(hm_binary_c, merge_legend = TRUE)

hm_binary <- Heatmap(wt_status_shared_mat_encoded,
                     name = "mutated",
                     
                     col = list("1" = "brown3", "0" = "blue", "2" = "#FFFFFF"),
                     heatmap_legend_param = list(border = "black",
                                                 labels = c("Mut", "WT", "NA")),
                     cluster_rows = FALSE,
                     show_row_dend = FALSE,
                     show_row_names = FALSE,
                     # row_order = row_order(hmap_binary_c),
                     
                     show_column_names = FALSE,
                     show_column_dend = FALSE,
                     cluster_columns = FALSE,
                     show_heatmap_legend = TRUE,
                     top_annotation = col_anno_polychrome
)

hmap_binary <- draw(hm_binary, merge_legend = TRUE)


hm_binary_clone <- Heatmap(wt_status_shared_mat_encoded,
                       name = "mutated",
                       
                       col = list("1" = "brown3", "0" = "blue", "2" = "#FFFFFF"),
                       heatmap_legend_param = list(border = "black",
                                                   labels = c("Mut", "WT", "NA")),
                       cluster_rows = TRUE,
                       show_row_dend = FALSE,
                       show_row_names = FALSE,
                       
                       show_column_names = FALSE,
                       show_column_dend = FALSE,
                       cluster_columns = FALSE,
                       show_heatmap_legend = TRUE,
                       column_split = clones_df$clone,
                       column_gap = unit(0, "mm"),
                       top_annotation = col_anno_kmeans
)
hmap_binary_clone <- draw(hm_binary_clone, merge_legend = TRUE)

pdf("mutation_shared_heatmap.pdf")
draw(hm_binary, merge_legend = TRUE)
draw(hm_binary_c, merge_legend = TRUE)
draw(hm_binary_clone, merge_legend = TRUE)
dev.off()

# all values ----
row_order_df <- case_when(rownames(wt_status_all_mat_encoded_reordered) %in% venn_filtered_germline@IntersectionSets$`01` == TRUE ~ "D28",
                          rownames(wt_status_all_mat_encoded_reordered) %in% venn_filtered_germline@IntersectionSets$`10` == TRUE ~ "D0",
                          rownames(wt_status_all_mat_encoded_reordered) %in% venn_filtered_germline@IntersectionSets$`11` == TRUE ~ "shared")

hm_binary_c_1 <- Heatmap(wt_status_all_mat_encoded_reordered,
                       name = "mutated",
                       
                       col = list("1" = "brown3", "0" = "blue", "2" = "#FFFFFF"),
                       heatmap_legend_param = list(border = "black",
                                                   labels = c("Mut", "WT", "NA")),
                       cluster_rows = FALSE,
                       show_row_dend = FALSE,
                       show_row_names = FALSE,
                       row_order = row_order(hmap_binary_c),
                       
                       show_column_names = FALSE,
                       show_column_dend = FALSE,
                       cluster_columns = FALSE,
                       show_heatmap_legend = TRUE,
                       top_annotation = col_anno_polychrome
)

hmap_binary_c_1 <- draw(hm_binary_c_1, merge_legend = TRUE)

hm_binary_c_2 <- Heatmap(wt_status_mat_encoded_1_reordered,
                         name = "mutated",
                         
                         col = list("1" = "brown3", "0" = "blue", "2" = "#FFFFFF"),
                         heatmap_legend_param = list(border = "black",
                                                     labels = c("Mut", "WT", "NA")),
                         cluster_rows = TRUE,
                         show_row_dend = FALSE,
                         show_row_names = FALSE,
                         # row_order = row_order(hmap_binary_c),
                         row_split = row_order_df,
                         row_gap = unit(0, "mm"),
                         cluster_row_slices = FALSE,
                         
                         show_column_names = FALSE,
                         show_column_dend = FALSE,
                         cluster_columns = FALSE,
                         show_heatmap_legend = TRUE,
                         top_annotation = col_anno_polychrome
)

hmap_binary_c_2 <- draw(hm_binary_c_2, merge_legend = TRUE)

pdf("germline_v2.pdf")
hmap_binary_c_1 <- draw(hm_binary_c_1, merge_legend = TRUE)
draw(hm_binary_c_2, merge_legend = TRUE)
hmap_vaf_germline_s_row_vaf_manual_kmeans_no_clones <- draw(
  hm_vaf_all_germline_s_row_manual_kmeans_no_clones, merge_legend = TRUE)

dev.off()

# integration and umap ------
merge_simple <- merge(d0, d28, add.cell.ids = c("d0", "d28"), merge.data = TRUE)
merge_list <- SplitObject(merge_simple, split.by = "orig.ident")
merge_list <- lapply(X = merge_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = merge_list, nfeatures = 3000)
merge_list <- lapply(X = merge_list, FUN = ScaleData, features = features)
merge_list <- lapply(X = merge_list, FUN = RunPCA, features = features)
anchors_rpca <- FindIntegrationAnchors(object.list = merge_list, 
                                       anchor.features = features, k.anchor = 5, reduction = "rpca")

rna_integrated_rpca <- IntegrateData(anchorset = anchors_rpca)
rna_integrated_rpca <- ScaleData(rna_integrated_rpca)
rna_integrated_rpca <- RunPCA(rna_integrated_rpca)
rna_integrated_rpca <- RunUMAP(rna_integrated_rpca, dims = 1:30)
DimPlot_scCustom(rna_integrated_rpca, group.by = "orig.ident", colors_use = polychrome_pal, pt.size = 1)

DefaultAssay(rna_integrated_rpca) <- "RNA"
our_AR <- c("AR", "KLK2", "KLK3", "TMPRSS2", "NKX3-1", "FKBP5", "PLPP1", "PMEPA1", "PART1", 
            "ALDH1A3", "STEAP4")
our_NE <- c("SYP", "NCAM1", "ENO2", "INSM1", "SOX2", "NKX2-1")

rna_integrated_rpca <- AddModuleScore(rna_integrated_rpca, features = list(our_AR), name = "AR_signature", search = TRUE)
rna_integrated_rpca <- AddModuleScore(rna_integrated_rpca, features = list(our_NE), name = "NE_signature", search = TRUE)


rna_integrated_rpca$clones <- "clone"
Idents(rna_integrated_rpca) <- "clones"
rna_integrated_rpca$clones <- clones_df$clone
rna_integrated_rpca$category <- factor(paste0(rna_integrated_rpca$orig.ident, " ", rna_integrated_rpca$clones))

# plotting ------
pdf("scrnaseq_d0_d28_1.pdf")

hmap_vaf_shared_s_row_c <- draw(hm_vaf_shared_s_row_c, merge_legend = TRUE)
hmap_vaf_germline_s_row_vaf_manual_kmeans <- draw(hm_vaf_germline_s_row_manual_kmeans, merge_legend = TRUE)
hmap_vaf_germline_s_row_vaf_manual_kmeans_no_clones <- draw(hm_vaf_germline_s_row_manual_kmeans_no_clones, merge_legend = TRUE)

DimPlot_scCustom(rna_integrated_rpca, group.by = "orig.ident", colors_use = c("#5A5156FF", "#3283FEFF"), pt.size = 0.7)

FeaturePlot_scCustom(rna_integrated_rpca, order = T, pt.size = 0.7, features = "AR_signature1") + 
  labs(x = "UMAP 1", y = "UMAP 2", title = "AR Signature")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 1), na.value = "lightgrey")
FeaturePlot(rna_integrated_rpca, order = T, pt.size = 0.7, features = "NE_signature1") +
  labs(x = "UMAP 1", y = "UMAP 2", title = "NE Signature")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 0.5), na.value = "lightgrey")
VlnPlot_scCustom(rna_integrated_rpca, "AR_signature1", group.by = "category", colors_use = stepped[c(1,2,5,6)], pt.size = 0)
VlnPlot_scCustom(rna_integrated_rpca, "NE_signature1", group.by = "category", colors_use = stepped[c(1,2,5,6)], pt.size = 0)

Meta_Highlight_Plot_man(rna_integrated_rpca, meta_data_column = "category", pt.size = 0.7,
                    meta_data_highlight = c("D0 Clone 1", "D0 Clone 2"), 
                    order = c("D0 Clone 2", "D0 Clone 1"), NavyAndOrange(flip_order = TRUE))
Meta_Highlight_Plot(rna_integrated_rpca, meta_data_column = "category", pt.size = 3, meta_data_highlight = c("D0 Clone 2"))
Meta_Highlight_Plot_man(rna_integrated_rpca, meta_data_column = "category", pt.size = 0.7,
                    meta_data_highlight = c("D28 Clone 1", "D28 Clone 2"), 
                    order = c("D28 Clone 1", "D28 Clone 2"), highlight_color = c("navy", "orange"))
Meta_Highlight_Plot(rna_integrated_rpca, meta_data_column = "category", pt.size = 3, meta_data_highlight = c("D28 Clone 1"))
dev.off()
