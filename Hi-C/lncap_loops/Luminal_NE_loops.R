library(tidyverse)
library(LoopRig)
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


# GO of just the genes overlapping -------
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

w0_genes <- get_seq2genes_from_anchors_symbol("w0_vs_w4.diffloop1.bedpe", 2)
dotplot(w0_genes[[4]]) + scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
  ggtitle("D0 specific loops 2809 genes")

w4_genes <- get_seq2genes_from_anchors_symbol("w0_vs_w4.diffloop2.bedpe", 2)
dotplot(w4_genes[[4]]) + scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
  ggtitle("D28 specific loops 2050 genes")



get_union_of_genes_from_anchors <- function(path, custom_cols = 0, distanceFilter = 5000) {
  # distance filter limits distance to TSS because annotatePeak will assign a gene no matter how far
  loops <- LoopsToRanges(path, custom_cols = custom_cols)
  loop_an <- paste0("Anchors_", seq(1,length(loops[[1]]$`Anchor 1`)))
  print(paste0("# anchors: ", length(loops[[1]]$`Anchor 1`)))
  loops[[1]]$`Anchor 1`$id <- loop_an
  loops[[1]]$`Anchor 2`$id <- loop_an
  anchor_1_anno <- annotatePeak(loops[[1]]$`Anchor 1`, tssRegion = c(-3000,3000), TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
  anchor_2_anno <- annotatePeak(loops[[1]]$`Anchor 2`, tssRegion = c(-3000,3000), TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
  
  if (is.null(distanceFilter) == TRUE) {
    print("No distance filter.")
    gene1 <- as.data.frame(anchor_1_anno)$SYMBOL
    print(paste0("# gene 1: ", length(gene1)))
    gene2 <- as.data.frame(anchor_2_anno)$SYMBOL
    print(paste0("# gene 2: ", length(gene2)))
    genes <- base::union(gene1, gene2)
  } else {
    gene1 <- as.data.frame(anchor_1_anno) %>% 
      filter(abs(distanceToTSS) <= distanceFilter) %>% dplyr::select(SYMBOL) %>% deframe()
    print(paste0("# gene 1: ", length(gene1)))
    gene2 <- as.data.frame(anchor_2_anno) %>% 
      filter(abs(distanceToTSS) <= distanceFilter) %>% dplyr::select(SYMBOL) %>% deframe()
    print(paste0("# gene 2 : ", length(gene2)))
    genes <- base::union(gene1, gene2)
    print(paste0("# final genes: ", length(genes)))
  }
  go <- enrichGO(gene = genes,
                 keyType = "SYMBOL", 
                 ont = "BP",
                 OrgDb = org.Hs.eg.db,
                 qvalueCutoff = 0.05)
  go_sym <- simplify(go)
  output_list <- list(genes, anchor_1_anno, anchor_2_anno, go, go_sym)
}
w0_121 <- get_union_of_genes_from_anchors("w0_vs_w4.diffloop1.bedpe", custom_cols = 2)
dotplot(w0_121[[5]]) + scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
  ggtitle("D0 specific loops 1809 genes")

w4_121 <- get_union_of_genes_from_anchors("w0_vs_w4.diffloop2.bedpe", custom_cols = 2)
dotplot(w4_121[[5]]) + scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
  ggtitle("D28 specific loops 1298 genes")

# Gene expression --------
deg <- read_csv("/projects/p20023/Viriya/analysis/foxa2/rna-seq/degenes_lrt.csv")

w0_genes_w_deg <- w0_genes_symbol[w0_genes_symbol %in% deg$...1]
go_w0 <- enrichGO(gene = w0_genes_w_deg,
               keyType = "SYMBOL", 
               ont = "BP",
               OrgDb = org.Hs.eg.db,
               qvalueCutoff = 0.05)
dotplot(go_w0)
go_sym_w0 <- simplify(go_w0)
dotplot(go_sym_w0) + scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
  ggtitle("3094 D0 specific loops 
1562 genes that are differentially expressed") + theme(plot.title=element_text(hjust=1))

# final --------
w0_genes_symbol <- w0_genes[[5]]$SYMBOL[!is.na(w0_genes[[5]]$SYMBOL)]

deg_w0 <- deg %>% subset(log2FoldChange < 0)
w0_genes_w_w0 <- w0_genes_symbol[w0_genes_symbol %in% deg_w0$...1]
go_w0_w0 <- enrichGO(gene = w0_genes_w_w0,
                     keyType = "SYMBOL", 
                     ont = "BP",
                     OrgDb = org.Hs.eg.db,
                     qvalueCutoff = 0.05)
dotplot(go_w0_w0)
go_sym_w0_w0 <- simplify(go_w0_w0)


deg_w0_1 <- deg %>% subset(log2FoldChange < -1)
w0_genes_w_w0_1 <- w0_genes_symbol[w0_genes_symbol %in% deg_w0_1$...1]
go_w0_w0_1 <- enrichGO(gene = w0_genes_w_w0_1,
                     keyType = "SYMBOL", 
                     ont = "BP",
                     OrgDb = org.Hs.eg.db,
                     qvalueCutoff = 0.05)
dotplot(go_w0_w0_1)
go_sym_w0_w0_1 <- simplify(go_w0_w0_1)


w4_genes_symbol <- w4_genes[[5]]$SYMBOL[!is.na(w4_genes[[5]]$SYMBOL)]
deg_w4 <- deg %>% subset(log2FoldChange > 0)
w4_genes_w_w4 <- w4_genes_symbol[w4_genes_symbol %in% deg_w4$...1]
go_w4_w4 <- enrichGO(gene = w4_genes_w_w4,
                     keyType = "SYMBOL", 
                     ont = "BP",
                     OrgDb = org.Hs.eg.db,
                     qvalueCutoff = 0.05)
dotplot(go_w4_w4)
go_sym_w4_w4 <- simplify(go_w4_w4)

deg_w4_1 <- deg %>% subset(log2FoldChange > 1)
w4_genes_w_w4_1 <- w4_genes_symbol[w4_genes_symbol %in% deg_w4_1$...1]
go_w4_w4_1 <- enrichGO(gene = w4_genes_w_w4_1,
                     keyType = "SYMBOL", 
                     ont = "BP",
                     OrgDb = org.Hs.eg.db,
                     qvalueCutoff = 0.05)
dotplot(go_w4_w4_1)
go_sym_w4_w4_1 <- simplify(go_w4_w4_1)


# Fig S2H
pdf("D0vD28_GO.pdf", width = 7, height = 10, useDingbats = FALSE)

dotplot(go_sym_w0_w0) + scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
  ggtitle("3094 D0 specific loops 
737 genes upregulated in D0") + theme(plot.title=element_text(hjust=1))

dotplot(go_sym_w4_w4) + scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
  ggtitle("1884 D28 specific loops 
800 genes upregulated in D28") + theme(plot.title=element_text(hjust=1))

dotplot(go_sym_w0_w0_1) + scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  ggtitle("3094 D0 specific loops 
356 genes upregulated in D0") + theme(plot.title=element_text(hjust=1))

dotplot(go_sym_w4_w4_1) + scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  ggtitle("1884 D28 specific loops 
441 genes upregulated in D28") + theme(plot.title=element_text(hjust=1))

dev.off()


