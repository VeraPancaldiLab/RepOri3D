RANDOM_range_keeping_distance_to_TSS <- function(range, seed = 1988) {
  
  
  # Chromosome selection
  range <- range[seqnames(range) %in% paste0("chr", c(1:19, "X"))]
  
  
  center_of_origin <- apply(cbind(start(range),end(range)), 1, function(x) round(mean(x)))
  range_center <- range
  start(range_center) <- center_of_origin
  end(range_center)   <- center_of_origin
  
  
  hits <- distanceToNearest(range_center, genes_starts)
  
  # Remove ranges that dissapeared 
  range_center <- range_center[queryHits(hits)]
  
  
  # min_distances<- hits@elementMetadata[[1]]
  
  direction_and_distance <- start(genes_starts[subjectHits(hits)]) - start(range_center[queryHits(hits)]) 
  
  direction_and_distance <- data.frame(dist = direction_and_distance, strand = genes_starts$Strand[subjectHits(hits)])
  
  
  direction_and_distance$dist_corr_by_strand <- direction_and_distance$dist
  direction_and_distance$dist_corr_by_strand[direction_and_distance$strand == -1] <- -direction_and_distance$dist_corr_by_strand[direction_and_distance$strand == -1]
  
  direction_and_distance$efficiency <- range_center$norm_background_corrected

  
  start_center_distance <- start(range) - start(range_center)
  end_center_distance   <- end(range) - end(range_center)
  
  # In each randomization i change to which gene they have such distance
  
  set.seed(seed)
  random_oris0 <- genes_starts[sample(x       = 1:length(genes_starts), 
                                      size    = length(range_center), 
                                      replace = T)]
  
  random_oris <- random_oris0
 
  random_oris <- as.data.frame(random_oris0, row.names = NULL)
  
  # negative direction means upstream, so it needs to be add it to the gene start ,and vice versa
  
  # +1 genes
  random_oris$start[random_oris$Strand == 1] <- random_oris$start[random_oris$Strand == 1] - direction_and_distance$dist_corr_by_strand[random_oris$Strand == 1]

  # -1 genes
  random_oris$start[random_oris$Strand == -1] <- random_oris$start[random_oris$Strand == -1] + direction_and_distance$dist_corr_by_strand[random_oris$Strand == -1]
  
  random_oris$end <- random_oris$start
  
  random_oris <- makeGRangesFromDataFrame(random_oris, seqnames.field = "seqnames", start.field =  "start", end.field = "end", ignore.strand = T)
  
 
  end(random_oris)   <- end(random_oris) + end_center_distance
  start(random_oris) <- start(random_oris) + start_center_distance

  random_oris
  
}










