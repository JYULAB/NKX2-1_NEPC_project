{
  library(LoopRig)
  library(tidyverse)
  library(ChIPseeker)
  library(clusterProfiler)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
}

# modifying LoopRig's LinkedLoops function to our purposes
LinkedLoops <- function (loop_ranges, element_ranges_x, element_ranges_y, range_out_x = FALSE, 
                         range_out_y = FALSE, overlap_threshold = 1) 
{
  if (class(loop_ranges) != "LoopRanges") {
    stop("Please enter an object of LoopRanges class for the loop_ranges parameter")
  }
  if (length(loop_ranges) != 1) {
    stop("Please enter a conseus LoopRanges object with only one range for the loop_ranges parameter")
  }
  if ((range_out_x == TRUE) & (range_out_y == TRUE)) {
    stop("Can only output either element_ranges_x or element_ranges_y as ElementRanges object")
  }
  loop_ranges <- loop_ranges[[1]]
  anchor1_overlaps_1 <- as.data.frame(findOverlaps(loop_ranges[[1]], 
                                                   element_ranges_x, minoverlap = overlap_threshold))
  anchor2_overlaps_1 <- as.data.frame(findOverlaps(loop_ranges[[2]], 
                                                   element_ranges_y, minoverlap = overlap_threshold))
  hits_merge_1 <- merge(anchor1_overlaps_1, anchor2_overlaps_1, 
                        by = "queryHits")
  anchor1_overlaps_2 <- as.data.frame(findOverlaps(loop_ranges[[2]], 
                                                   element_ranges_x, minoverlap = overlap_threshold))
  anchor2_overlaps_2 <- as.data.frame(findOverlaps(loop_ranges[[1]], 
                                                   element_ranges_y, minoverlap = overlap_threshold))
  hits_merge_2 <- merge(anchor1_overlaps_2, anchor2_overlaps_2, 
                        by = "queryHits")
  if (range_out_x == TRUE) {
    range_x_hits <- unique(c(hits_merge_1[[2]], hits_merge_2[[2]]))
    range_x_subset <- unique(element_ranges_x[range_x_hits])
    return(structure(list(el_x_linked = range_x_subset), 
                     class = "ElementRanges"))
  }
  if (range_out_y == TRUE) {
    range_y_hits <- unique(c(hits_merge_1[[3]], hits_merge_2[[3]]))
    range_y_subset <- unique(element_ranges_y[range_y_hits])
    return(structure(list(el_y_linked = range_y_subset), 
                     class = "ElementRanges"))
  }
  
  combined_query_hits <- unique(c(hits_merge_1[[1]], hits_merge_2[[1]]))
  an_1 <- loop_ranges[[1]][combined_query_hits]
  an_2 <- loop_ranges[[2]][combined_query_hits]
  
  df1 <- as.data.frame(an_1)[,1:3]
  df2 <- as.data.frame(an_2)[,1:3]
  df <- cbind(df1,df2)
  return(df)
  
}

make_loops <- function(week, fx, nx, path = "./") {
  print(paste0("../loops/", week, "_mustache_10kb.bedpe"))
  loops <- LoopsToRanges(paste0("../loops/", week, "_mustache_10kb.bedpe"), custom_cols = 2)
  
  elements_1 <- ElementsToRanges(paste0("../../hi-c/peaks/", fx, "_peaks.narrowPeak"), element_names = c("FOXA2"), 
                                 custom_cols = 7, custom_mcols = 4)
  elements_2 <- ElementsToRanges(paste0("../../hi-c/peaks/", nx, "_peaks.narrowPeak"), element_names = c("NKX2-1"), 
                                 custom_cols = 7, custom_mcols = 4)
  
  linked_loops <- LinkedLoops(loop_ranges = loops, element_ranges_x = elements_1[[1]], 
                              element_ranges_y = elements_2[[1]])
  write_delim(linked_loops, paste0(path, week ,"_fx_nx_linked.bedpe"), col_names = FALSE, delim = "\t")
  return(linked_loops)
}

w4_loops_linked <- make_loops(week = "w4", fx = "m921", nx = "m929", "./10kb/")
w3_loops_linked <- make_loops(week = "w3", fx = "m919", nx = "m927", "./10kb/")
w2_loops_linked <- make_loops(week = "w2", fx = "m917", nx = "m925", "./10kb/")
lncap_loops_linked <- make_loops(week = "w0", fx = "m914", nx = "m923", "./10kb/")

