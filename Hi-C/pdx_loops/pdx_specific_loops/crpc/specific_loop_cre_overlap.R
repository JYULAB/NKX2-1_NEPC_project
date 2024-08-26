library(tidyverse)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(LoopRig)
library(rtracklayer)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

ad_cre <- readPeakFile("../../NEPC_recurrent_2/pp/AD_hg38.bed")
ne_cre <- readPeakFile("../../NEPC_recurrent_2/pp/NE_hg38.bed")

loops <- LoopsToRanges("../results/CRPC_loop.bedpe", custom_cols = 2)
loop_an <- paste0("Anchors_", seq(1,length(loops[[1]]$`Anchor 1`)))
print(paste0("# anchors: ", length(loops[[1]]$`Anchor 1`)))
loops[[1]]$`Anchor 1`$id <- loop_an
loops[[1]]$`Anchor 2`$id <- loop_an

count_peak_in_anchors <- function(peaks, loops, sample = NULL) {
  hits1 <- findOverlaps(peaks, loops[[1]]$`Anchor 1`)
  hits2 <- findOverlaps(peaks, loops[[1]]$`Anchor 2`)
  query <- calc_hit_stats_of_list_of_hits(list(hits1, hits2))
  l <- list("hits1" = hits1, "hits2" = hits2, "peaks" = query)
}

calc_hit_stats_of_list_of_hits <- function(hit_list) {
  print(paste0("Query hits 1: ", length(unique(queryHits(hit_list[[1]]))), 
               ", Query hits 2: ", length(unique(queryHits(hit_list[[2]])))))
  query_unique <- union(unique(queryHits(hit_list[[1]])), unique(queryHits(hit_list[[2]])))
  print(paste0("Total unique query hits: ", length(query_unique)
  ))
  print(paste0("Subject hits 1: ", length(unique(subjectHits(hit_list[[1]]))), 
               ", Subject hits 2: ", length(unique(subjectHits(hit_list[[2]])))))
  print(paste0("Total unique subject hits: ", length(union(unique(subjectHits(hit_list[[1]])), 
                                                           unique(subjectHits(hit_list[[2]]))))
  ))
  return(query_unique)
}

hits_atac_in_anchors_ne <- count_peak_in_anchors(ne_cre, loops)
granges_atac_in_anchors_ne <- ne_cre[hits_atac_in_anchors_ne$peaks,]
export(granges_atac_in_anchors_ne, paste0("cre_anchor_only/", "NE_cre_at_AD_recurrent.bed"))

loops_with_ne_cre <- union(subjectHits(hits_atac_in_anchors_ne$hits1), subjectHits(hits_atac_in_anchors_ne$hits2))
loops_ne_cre <- cbind(as.data.frame(loops[[1]]$`Anchor 1`[loops_with_ne_cre])[1:3], 
                      as.data.frame(loops[[1]]$`Anchor 2`[loops_with_ne_cre])[1:3])
write_delim(loops_ne_cre, "AD_recurrent_loops_with_NE_cre.bedpe", col_names = FALSE, delim = "\t")

hits_atac_in_anchors_ad <- count_peak_in_anchors(ad_cre, loops)
granges_atac_in_anchors_ad <- ad_cre[hits_atac_in_anchors_ad$peaks,]
export(granges_atac_in_anchors_ad, paste0("cre_anchor_only/", "AD_cre_at_AD_recurrent.bed"))

loops_with_ad_cre <- union(subjectHits(hits_atac_in_anchors_ad$hits1), subjectHits(hits_atac_in_anchors_ad$hits2))
loops_ad_cre <- cbind(as.data.frame(loops[[1]]$`Anchor 1`[loops_with_ad_cre])[1:3], 
                      as.data.frame(loops[[1]]$`Anchor 2`[loops_with_ad_cre])[1:3])
write_delim(loops_ad_cre, "AD_recurrent_loops_with_AD_cre.bedpe", col_names = FALSE, delim = "\t")
