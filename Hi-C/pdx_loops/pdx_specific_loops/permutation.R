library(tidyverse)
library(regioneR)
library(LoopRig)
library(rtracklayer)
library(GenomicRanges)
set.seed(1234)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

library(rtracklayer)
ad_cre <- import.bed("../../round2/specific_loops/NEPC_recurrent_2/pp/AD_hg38.bed")
ne_cre <- import.bed("../../round2/specific_loops/NEPC_recurrent_2/pp/NE_hg38.bed")

nepc_loops <- LoopsToRanges("NEPC_loop.bedpe", custom_cols = 2)[[1]]
nepc_loop_an <- paste0("Anchors_", seq(1,length(nepc_loops$`Anchor 1`)))
print(paste0("# anchors: ", length(nepc_loops$`Anchor 1`)))
nepc_loops$`Anchor 1`$id <- nepc_loop_an
nepc_loops$`Anchor 2`$id <- nepc_loop_an

crpc_loops <- LoopsToRanges("CRPC_loop.bedpe", custom_cols = 2)[[1]]
crpc_loop_an <- paste0("Anchors_", seq(1,length(crpc_loops$`Anchor 1`)))
print(paste0("# anchors: ", length(crpc_loops$`Anchor 1`)))
crpc_loops$`Anchor 1`$id <- crpc_loop_an
crpc_loops$`Anchor 2`$id <- crpc_loop_an

# functions -----
loop_overlap <- function(peaks, loops) {
  hits1 <- findOverlaps(peaks, loops[[1]])
  hits2 <- findOverlaps(peaks, loops[[2]])
  unique_loop_overlap <- length(union(unique(subjectHits(hits1)), 
                                      unique(subjectHits(hits2))))
  return(unique_loop_overlap)
}

loop_focused <- function(cre, loop, num_randomizations = 50) {
  expected_overlaps <- numeric(num_randomizations)
  for (i in 1:num_randomizations) {
    rand_loops <- GRangesList(circularRandomizeRegions(loop$`Anchor 1`, genome = "hg38"),
                              circularRandomizeRegions(loop$`Anchor 2`, genome = "hg38"))
    rand_overlap <- loop_overlap(loops = rand_loops, peaks = cre)
    expected_overlaps[i] <- rand_overlap
  }
  observed_overlap <- cre_overlap(cre, loop)
  mean_expected_overlap <- mean(expected_overlaps)
  std_expected_overlap <- sd(expected_overlaps)
  z_score <- (observed_overlap - mean_expected_overlap) / std_expected_overlap
  p_value <- (sum(expected_overlaps >= observed_overlap) + 1) / (num_randomizations + 1)
  cat("Observed Overlap:", observed_overlap, "\n")
  cat("Mean Expected Overlap:", mean_expected_overlap, "\n")
  cat("Standard Deviation of Expected Overlap:", std_expected_overlap, "\n")
  cat("Z-Score:", z_score, "\n")
  cat("P-Value:", p_value, "\n")
  return(list(z_score, p_value))
}

# manual calculation ----
crpc_loop_ad_cre <- loop_focused(cre = ad_cre, loop = crpc_loops, n = 1000)
crpc_loop_ne_cre <- loop_focused(cre = ne_cre, loop = crpc_loops, n = 1000)
nepc_loop_ne_cre <- loop_focused(cre = ne_cre, loop = nepc_loops, n = 1000)
nepc_loop_ad_cre <- loop_focused(cre = ad_cre, loop = nepc_loops, n = 1000)
