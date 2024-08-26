library(tidyverse)
library(UpSetR)
library(ComplexHeatmap)

samples <- c("93", "1451", "1452")

# this will combine every comparison except those unique to just one sample
create_consensus_bedpe <- function(sample) {
  venn <- read.table(paste0("../NCI", "/output_",sample, "_diffloop1.bedpe"), header = TRUE, sep = "\t", skip = 6)
  colnames(venn) <- stringr::str_remove(colnames(venn), ".diffloop1.bedpe")
  colnames(venn) <- stringr::str_remove(colnames(venn), "X")
  
  venn$alpha <- apply(venn[, -1], 1, function(x) {
    sets <- letters[1:16]
    selected_sets <- sets[which(x == "X")]
    paste(selected_sets, collapse = "")
  })
  
  chosen_loops <-  venn[nchar(venn$alpha) >= 2,]$Name
  print(length(chosen_loops))
  chosen_loops <- str_replace_all(chosen_loops, "\\|", "_")
  
  chosen_loops_files <- paste0("../NCI", "/", sample, "_gain_", chosen_loops)
  chosen_loops_list <- lapply(chosen_loops_files, read_delim, skip = 1, col_names = FALSE)
  chosen_loops_df <- do.call(rbind, chosen_loops_list)
  print(dim(chosen_loops_df))
  write_delim(chosen_loops_df, paste0(sample, "_consensus_loops.bedpe"), delim = "\t", col_names = FALSE)
}

lapply(samples, create_consensus_bedpe)
create_consensus_bedpe("nci")

# after homer of consensus peaks
library(Vennerable)
library(tidyverse)

venn <- read.table("results/output_nepc.bedpe", header = TRUE, sep = "\t", skip = 6)
venn$alpha <- apply(venn[, -1], 1, function(x) {
  sets <- c("A", "B", "C", "D")
  selected_sets <- sets[which(x == "X")]
  paste(selected_sets, collapse = "")
})

order_vector <- c("0", "A", "B", "AB", "C", "AC", "BC", "ABC", "D", "AD", "BD", "ABD", "CD", "ACD", "BCD", "ABCD")

venn$alpha <- factor(venn$alpha, levels = order_vector)

# Sort the data frame based on the factor levels
venn <- venn[order(venn$alpha), ]
venn_plot <- Venn(SetNames = c("93", "145.1", "145.2", "nci"), 
                  Weight = c(0, venn$Features))
plot(venn_plot, type = "ellipses")

chosen_loops <- venn[!venn$alpha %in% c("A", "B", "C", "D"),]$Name
chosen_loops <- str_replace_all(chosen_loops, "\\|", "_")
chosen_loops_files <- paste0("results/nepc_gain_", chosen_loops)

chosen_loops_list <- lapply(chosen_loops_files, read_delim, skip = 1, col_names = FALSE)
chosen_loops_df <- do.call(rbind, chosen_loops_list)
write_delim(chosen_loops_df, "results/NEPC_loop.bedpe", delim = "\t", col_names = FALSE)


library(LoopRig)
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

get_seq2genes_from_anchors_symbol <- function(path, custom_cols = 0) {
  loops <- LoopsToRanges(path, custom_cols = custom_cols)
  loop_an <- paste0("Anchors_", seq(1,length(loops[[1]]$`Anchor 1`)))
  print(paste0("# anchors: ", length(loops[[1]]$`Anchor 1`)))
  loops[[1]]$`Anchor 1`$id <- loop_an
  loops[[1]]$`Anchor 2`$id <- loop_an
  anchor_1_anno <- seq2gene(loops[[1]]$`Anchor 1`, tssRegion = c(-3000,3000), flankDistance = 3000, TxDb = txdb)
  anchor_2_anno <- seq2gene(loops[[1]]$`Anchor 2`, tssRegion = c(-3000,3000), flankDistance = 3000, TxDb = txdb)
  t <- union(anchor_1_anno, anchor_2_anno) %>% unique()
  t_sym <- AnnotationDbi::select(org.Hs.eg.db, keys = t, columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")
  sym <- t_sym$SYMBOL[!is.na(t_sym$SYMBOL)]
  print(length(sym))
  go <- enrichGO(gene = sym,
                 keyType = "SYMBOL", 
                 ont = "BP",
                 OrgDb = org.Hs.eg.db,
                 qvalueCutoff = 0.05)
  go_sim <- simplify(go)
  output_list <- list(anchor_1_anno, anchor_2_anno, go, go_sim, t_sym)
}

nepc_symbols <- get_seq2genes_from_anchors_symbol("results/NEPC_loop.bedpe", custom_cols = 2)
dotplot(nepc_symbols[[4]]) + scale_y_discrete(labels = function(x) str_wrap(x, width = 30))

# annotate loop
path = "results/NEPC_loop.bedpe"
bedpe_with_gene <- function(path, custom_cols = 0) {
  loops <- LoopsToRanges(path, custom_cols = 2)
  loop_an <- paste0("Anchors_", seq(1,length(loops[[1]]$`Anchor 1`)))
  print(paste0("# anchors: ", length(loops[[1]]$`Anchor 1`)))
  loops[[1]]$`Anchor 1`$id <- loop_an
  loops[[1]]$`Anchor 2`$id <- loop_an
  anchor_1_anno <- annotatePeak(loops[[1]]$`Anchor 1`, tssRegion = c(-3000,3000), TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
  anchor_2_anno <- annotatePeak(loops[[1]]$`Anchor 2`, tssRegion = c(-3000,3000), TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
  write_delim(as.data.frame(anchor_1_anno@anno), "NEPC_anchor1_anno_121.bed", delim = "\t")
  write_delim(as.data.frame(anchor_2_anno@anno), "NEPC_anchor2_anno_121.bed", delim = "\t")
  
  anchor_1_seq_list <- lapply(seq(1,2), function(i) {
    seq2gene(loops[[1]][[1]][i], tssRegion = c(-3000,3000), flankDistance = 3000, TxDb = txdb)
  })
  
  
  anchor_1_anno_seq <- sapply(loops[[1]][[1]], seq2gene, tssRegion = c(-3000,3000), flankDistance = 3000, TxDb = txdb)
  anchor_1_anno_seq <- seq2gene(loops[[1]][[1]][3], tssRegion = c(-3000,3000), flankDistance = 3000, TxDb = txdb)
  anchor_2_anno <- seq2gene(loops[[1]]$`Anchor 2`, tssRegion = c(-3000,3000), flankDistance = 3000, TxDb = txdb)
  # t <- union(anchor_1_anno, anchor_2_anno) %>% unique()
  t_1  <- AnnotationDbi::select(org.Hs.eg.db, keys = , columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")
  sym <- t_sym$SYMBOL[!is.na(t_sym$SYMBOL)]
  print(length(sym))
 
  output_list <- list(anchor_1_anno, anchor_2_anno, go, go_sym, t_sym)
}

# gene expression -------
fpkm_df <- read_csv("../../../../cooltools/AB/gene_expression/samples_8_fpkm.csv")
fpkm_df <- column_to_rownames(fpkm_df, "...1")

nepc_recurrent_genes <- nepc_symbols[[5]]$SYMBOL[!is.na(nepc_symbols[[5]]$SYMBOL)]     
nepc_recurrent_genes_expression <- fpkm_df[nepc_recurrent_genes,] %>% drop_na()

# Create a function to perform paired t-test for a gene
t_test_LFC <- function(gene_expression) {
  g1 <- gene_expression[c(1,2,10:15)] # NEPC
  g2 <- gene_expression[c(3:9,16,17)] # nepc
  t_test_result <- t.test(g1, g2)
  LFC <- log2(sum(g1)/sum(g2))
  return(c(t_statistic = t_test_result$statistic, p_value = t_test_result$p.value, lfc = LFC))
}

# Apply the paired_t_test function to each row (gene) of the dataframe
result_matrix <- t(apply(nepc_recurrent_genes_expression, 1, t_test_LFC))

# Create a new dataframe to store the results
results_df <- data.frame(Gene = rownames(nepc_recurrent_genes_expression), t_statistic = result_matrix[, 1], p_value = result_matrix[, 2], lfc = result_matrix[, 3])
results_df <- results_df %>% filter(is.finite(lfc))
results_df$status <- ifelse(results_df$lfc > 1 & results_df$p_value < 0.05, "Up", 
                            ifelse(results_df$lfc < -1 & results_df$p_value < 0.05, "Down", "NS"))
write_csv(results_df, "NEPC_recurrent_genes.csv")

# Fig 1C ----
pdf("NEPC_volcano_plot.pdf", width = 5, height = 3)
ggplot(results_df, aes(x = lfc, y = -log10(p_value))) +
  geom_point(aes(color = status), show.legend = FALSE) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_classic() +
  labs(
    x = "Log2 Fold Change",
    y = "-log10(p-value)",
    subtitle = "Total 3015 genes overlapping with NEPC recurrent loops,
    Upregulated 677, Downregulated 203")
dev.off()
dim(results_df)

up <- results_df %>% filter(status == "Up")
down <- results_df %>% filter(status == "Down")

go_up <- enrichGO(gene = up$Gene,
                  keyType = "SYMBOL", 
                  ont = "BP",
                  OrgDb = org.Hs.eg.db,
                  qvalueCutoff = 0.05)
dotplot(go_up)

go_up_sym <- simplify(go_up)

# Fig 1E ----
pdf("NEPC_GO.pdf", width = 7, height = 10, useDingbats = FALSE)
dotplot(go_up_sym)  + scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
  ggtitle("677 upregulated genes in 
          NEPC recurrent loops")
dev.off()
cygo_up_sym@result

go_down <- enrichGO(gene = down$Gene,
                    keyType = "SYMBOL", 
                    ont = "BP",
                    OrgDb = org.Hs.eg.db,
                    qvalueCutoff = 0.05)
dotplot(go_down)
go_down_sym <- simplify(go_down)
dotplot(go_down_sym)

