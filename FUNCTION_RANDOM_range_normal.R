

generate_random_resp_oris_normal <- function(range, seed = 1988, mm ) {
  
  range <- range[seqnames(range) %in% paste0("chr", c(1:19, "X"))]
  
  set.seed(1988+seed)
  
  random_int0 <- random.intervals(chromosomes_no_shade_areas,
                                  n=length(range),
                                  ms=width(range)-1)
  
  seqlengths(random_int0) <- seqlengths(chromosomes_no_shade_areas)
  
  
  random_int <- GRanges( seqnames = chromosomes_no_shade_areas$chr[as.integer(seqnames(random_int0))],
                         ranges   = IRanges(start  = start(random_int0) + start(chromosomes_no_shade_areas)[as.integer(seqnames(random_int0))],
                                            end    = end(random_int0) + start(chromosomes_no_shade_areas)[as.integer(seqnames(random_int0))]))
  
  seqlevelsStyle(random_int) <- "UCSC"
  random_int
}
