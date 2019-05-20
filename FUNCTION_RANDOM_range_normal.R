# generate_random_resp_oris_normal <- function(range, seed = 1988, mm ) {
#   
#   range <- range[seqnames(range) %in% paste0("chr", c(1:19, "X"))]
#   
#   chromosomes_RANGES <- get(paste0("mm", mm, "_chromosomes_RANGES"))
#   seqlevels(range) <- seqlevels(chromosomes_RANGES)
#   seqlengths(range) <- seqlengths(chromosomes_RANGES)
#   
#   set.seed(1988+seed)
#   
#   random_int0 <- random.intervals(chromosomes_RANGES,
#                                   n=length(range),
#                                   ms=width(range)-1) ## THIS -1 WAS ADDED ON THE 9TH of may, because I realised that otherwise the ranges were 1bp longer
#   
#   seqlengths(chromosomes_RANGES) - seqlengths(random_int0)
#   
#   
#   seqlengths(random_int0) <- seqlengths(chromosomes_RANGES)
#   
#   
#   random_int <- GRanges( seqnames = as.character(seqnames(random_int0)),
#                          ranges   = IRanges(start  = start(random_int0) + start(chromosomes_RANGES)[1],
#                                             end    = end(random_int0) + start(chromosomes_RANGES)[1]))
#   
#   random_int 
# }


generate_random_resp_oris_normal <- function(range, seed = 1988, mm ) {
  
  range <- range[seqnames(range) %in% paste0("chr", c(1:19, "X"))]
  
  set.seed(1988+seed)
  
  random_int0 <- random.intervals(chromosomes_no_shade_areas,
                                  n=length(range),
                                  ms=width(range)-1) ## THIS -1 WAS ADDED ON THE 9TH of may, because I realised that otherwise the ranges were 1bp longer
  
  seqlengths(random_int0) <- seqlengths(chromosomes_no_shade_areas)
  
  
  random_int <- GRanges( seqnames = chromosomes_no_shade_areas$chr[as.integer(seqnames(random_int0))],
                         ranges   = IRanges(start  = start(random_int0) + start(chromosomes_no_shade_areas)[as.integer(seqnames(random_int0))],
                                            end    = end(random_int0) + start(chromosomes_no_shade_areas)[as.integer(seqnames(random_int0))]))
  
  seqlevelsStyle(random_int) <- "UCSC"
  random_int
}
