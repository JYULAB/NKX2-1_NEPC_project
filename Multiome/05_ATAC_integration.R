library(tidyverse)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(future)
library(stringr)
library(patchwork)
library(scCustomize)
library(EnsDb.Hsapiens.v86)
library(data.table)
set.seed(1234)

d21 <- readRDS("../single_d21/simple_d21.RDS")
d14 <- readRDS("../single/simple_d14.RDS")
d2 <- readRDS("../single_d2/simple_d2.RDS")

DefaultAssay(d14) <- "ATAC"
DefaultAssay(d2) <- "ATAC"
DefaultAssay(d21) <- "ATAC"

# using the singly called macs2 peaks
mac_d14 <- readRDS("../single/d14_macs2_counts.RDS")
mac_d2 <- readRDS("../single_d2/macs2_counts")
mac_d21 <- readRDS("../single_d21/macs2_counts")

peaks_d14 <- str_split(mac_d14@Dimnames[[1]], "-", simplify = TRUE) %>% as.data.frame()
colnames(peaks_d14) <- c("chr", "start", "end")

peaks_d2 <- str_split(mac_d2@Dimnames[[1]], "-", simplify = TRUE) %>% as.data.frame()
colnames(peaks_d2) <- c("chr", "start", "end")

peaks_d21 <- str_split(mac_d21@Dimnames[[1]], "-", simplify = TRUE) %>% as.data.frame()
colnames(peaks_d21) <- c("chr", "start", "end")

gr_d14 <- makeGRangesFromDataFrame(peaks_d14)
gr_d2 <- makeGRangesFromDataFrame(peaks_d2)
gr_d21 <- makeGRangesFromDataFrame(peaks_d21)

combined_peaks <- reduce(x = c(gr_d14, gr_d2, gr_d21))
peakwidths <- width(combined_peaks)
range(peakwidths) # 224 to 2102 original arc ranger output, #200 to 4039 for above procedure

d14_counts <- FeatureMatrix(
  fragments = Fragments(d14),
  features = combined_peaks,
  cells = colnames(d14)
)
saveRDS(d14_counts, "d14_counts.RDS")


d2_counts <- FeatureMatrix(
  fragments = Fragments(d2),
  features = combined_peaks,
  cells = colnames(d2)
)
saveRDS(d2_counts, "d2_counts.RDS")

d21_counts <- FeatureMatrix(
  fragments = Fragments(d21),
  features = combined_peaks,
  cells = colnames(d21)
)
saveRDS(d21_counts, "d21_counts.RDS")

d2 <- RenameCells(d2, add.cell.id = "d2")
d21 <- RenameCells(d21, add.cell.id = "d21")
d14 <- RenameCells(d14, add.cell.id = "d14")

d14_assay <- CreateChromatinAssay(d14_counts, fragments = Fragments(d14))
d14_new <- CreateSeuratObject(d14_assay, assay = "ATAC")
d14_new[["RNA"]] <- d14[["RNA"]]
d14_new[["SCT"]] <- d14[["SCT"]]
d14_new$dataset <- "d14"
saveRDS(d14_new, "d14.RDS")

d2_assay <- CreateChromatinAssay(d2_counts, fragments = Fragments(d2))
d2_new <- CreateSeuratObject(d2_assay, assay = "ATAC")
d2_new[["RNA"]] <- d2[["RNA"]]
d2_new[["SCT"]] <- d2[["SCT"]]
d2_new$dataset <- "d2"
saveRDS(d2_new, "d2.RDS")

d21_assay <- CreateChromatinAssay(d21_counts, fragments = Fragments(d21))
d21_new <- CreateSeuratObject(d21_assay, assay = "ATAC")
d21_new[["RNA"]] <- d21[["RNA"]]
d21_new[["SCT"]] <- d21[["SCT"]]
d21_new$dataset <- "d21"
saveRDS(d21_new, "d21.RDS")

d2_new <- RenameCells(d2_new, add.cell.id = "d2")
d21_new <- RenameCells(d21_new, add.cell.id = "d21")
d14_new <- RenameCells(d14_new, add.cell.id = "d14")

d2_new <- FindTopFeatures(d2_new, min.cutoff = 10)
d2_new <- RunTFIDF(d2_new)
d2_new <- RunSVD(d2_new)
saveRDS(d2_new, "d2.RDS")

d21_new <- FindTopFeatures(d21_new, min.cutoff = 10)
d21_new <- RunTFIDF(d21_new)
d21_new <- RunSVD(d21_new)
saveRDS(d21_new, "d21.RDS")

d14_new <- FindTopFeatures(d14_new, min.cutoff = 10)
d14_new <- RunTFIDF(d14_new)
d14_new <- RunSVD(d14_new)
saveRDS(d14_new, "d14.RDS")

d2_new <- readRDS("../../multiome/atac_integration_2/d2.RDS")
d14_new <- readRDS("../../multiome/atac_integration_2/d14.RDS")
d21_new <- readRDS("../../multiome/atac_integration_2/d21.RDS")

# combined object -----
combined <- merge(
  x = d2_new,
  y = list(d14_new, d21_new)
)

combined[["ATAC"]]
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
DimPlot(combined, group.by = "dataset") + ggtitle("ATAC Combined")
saveRDS(combined, "combined.RDS")
combined <- readRDS("../../multiome/atac_integration_2/combined.RDS")

FeaturePlot_scCustom(combined, "chr3-36992113-36993926")

integration.anchors <- FindIntegrationAnchors(
  object.list = list(d2_new, d14_new, d21_new),
  anchor.features = rownames(d2_new),
  reduction = "rlsi",
  dims = 2:50
)

integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:50
)

integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
integrated <- FindNeighbors(integrated, reduction = "integrated_lsi", dims = 2:30)
integrated <- FindClusters(integrated, graph.name = "ATAC_snn")
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated)

saveRDS(integrated, "atac_integrated.RDS")
integrated <- readRDS("atac_integrated.RDS")

integrated <- FindClusters(integrated, graph.name = "ATAC_snn", resolution = 0.4)
integrated <- FindClusters(integrated, graph.name = "ATAC_snn", resolution = 0.6)

DimPlot(integrated, group.by = "dataset") + ggtitle("ATAC")
DimPlot(integrated, group.by = "ATAC_snn_res.0.4", label = TRUE, label.box = TRUE) + ggtitle("ATAC, res = 0.4")
DimPlot(integrated, group.by = "ATAC_snn_res.0.6", label = TRUE, label.box = TRUE) + ggtitle("ATAC, res = 0.6")

# ATAC AR/NE signatures ------

atac_anno <- Annotation(integrated[["ATAC"]])

transcript <- CollapseToLongestTranscript(atac_anno)
transcript <- Extend(x = transcript, upstream = 5000)

our_AR <- c("AR", "KLK2", "KLK3", "TMPRSS2", "NKX3-1", "FKBP5", "PLPP1", "PMEPA1", "PART1", 
            "ALDH1A3", "STEAP4")
our_NE <- c("SYP", "NCAM1", "ENO2", "INSM1", "SOX2", "NKX2-1")

transcript_AR <- transcript[transcript$gene_name %in% our_AR]
transcript_NE <- transcript[transcript$gene_name %in% our_NE]

make_peaks <- function(transcript) {
  peaks_in_promoter_hits <- findOverlaps(integrated[["ATAC"]], transcript)
  peaks_in_promoter <- integrated[["ATAC"]]@ranges[queryHits(peaks_in_promoter_hits),]
  peaks_in_promoter$gene <- transcript[subjectHits(peaks_in_promoter_hits),]$gene_name
  peaks_in_promoter$peak <- as.character(peaks_in_promoter) %>% str_replace(pattern = ":", replacement = "-")
  return(peaks_in_promoter)
 }

promoter_peak_AR <- make_peaks(transcript_AR)
promoter_peak_NE <- make_peaks(transcript_NE)

modules <- list("AR.sig" = promoter_peak_AR$peak,
                "NE.sig" = promoter_peak_NE$peak)

integrated <- AddChromatinModule(integrated, features = modules, genome = BSgenome.Hsapiens.UCSC.hg38)

# RNA AR/NE signatures ----
DefaultAssay(integrated) <- "RNA"
integrated <- AddModuleScore(integrated, features = list(our_AR), name = "AR_signature", search = TRUE)
integrated <- AddModuleScore(integrated, features = list(our_NE), name = "NE_signature", search = TRUE)

# final pdfs ----
Idents(integrated) <- integrated$dataset
Idents(integrated) <- factor(Idents(integrated), levels = c('d2', 'd14', 'd21'))
integrated$dataset <- Idents(integrated)

pdf("scatac.pdf")
DimPlot_scCustom(integrated, group.by = "dataset", colors_use = ColorBlind_Pal(), pt.size = 0.7) + 
  ggtitle("ATAC") &
  labs(x = "atac_UMAP 1", y = "atac_UMAP 2")

FeaturePlot(integrated, order = T, pt.size = 0.7, features = "AR.sig") + 
  labs(x = "atac_UMAP 1", y = "atac_UMAP 2", title = "AR Signature (ATAC)")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 6), na.value = "lightgrey")
FeaturePlot(integrated, order = T, pt.size = 0.7, features = "NE.sig") +
  labs(x = "atac_UMAP 1", y = "atac_UMAP 2", title = "NE Signature (ATAC)")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 6), na.value = "lightgrey")

FeaturePlot(integrated, order = T, pt.size = 0.7, features = "AR_signature1") + 
  labs(x = "atac_UMAP 1", y = "atac_UMAP 2", title = "AR Signature (RNA)")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 1), na.value = "lightgrey")
FeaturePlot(integrated, order = T, pt.size = 0.7, features = "NE_signature1") +
  labs(x = "atac_UMAP 1", y = "atac_UMAP 2", title = "NE Signature (RNA)")  + 
  scale_color_gradientn(colours = c('lightgrey', 'blue'), limits = c(0, 0.5), na.value = "lightgrey")

dev.off()

