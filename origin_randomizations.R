setwd("~/Jodkowska_Pancaldi/")

library(GenomicRanges)
library( seqbias ) # for random.intervals()
library( rtracklayer) # for import.chain()

output_folder <- "NORMAL_and_TSS_randomizations/"
if(!file.exists(paste0(output_folder, "randomized_mm10")))       dir.create("randomized_mm10")


# Load chromosomes without low-mappability regions
load("mm10chromosomes_no_shade_areas.RData", verbose = T)
chromosomes_no_shade_areas      <- mm10chromosomes_no_shade_areas
chromosomes_no_shade_areas_UCSC <- chromosomes_no_shade_areas
seqlevelsStyle(chromosomes_no_shade_areas_UCSC) <- "UCSC"

# Gene coordinates in mm10
library("biomaRt")
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl",  host = "aug2017.archive.ensembl.org")

g <- getBM(mart = ensembl, attributes = c("ensembl_gene_id",
                                          "chromosome_name",
                                          "start_position", 
                                          "end_position", 
                                          "strand", 
                                          "gene_biotype"))
genes <- g[g$gene_biotype == "protein_coding",]
genes <- genes[genes$chromosome_name %in% c(1:19, "X"),]


genes_RANGES <- GRanges(seqnames = paste0("chr", genes$chromosome_name),
                        ranges   = IRanges(start = genes$start_position,
                                           end = genes$end_position),
                        Strand = genes$strand)


genes_starts <- genes_RANGES
# STRAND 1
end(genes_starts)[genes_starts$Strand == 1] <- start(genes_starts)[genes_starts$Strand == 1]
# strand -1
start(genes_starts)[genes_starts$Strand == -1] <- end(genes_starts)[genes_starts$Strand == -1]



### FUNCTIONS to randomize origins
source("FUNCTION_RANDOM_range_normal.R") # To randomize anywhere in the genome except for low-mappability regions
source("FUNCTION_RANDOM_range_keeping_distance_to_TSS.R") # To randomize preserving distance and direction to gene start

load("example_oris.RData", verbose = T)


for(oris in c("example_oris")) {
  oris_ranges <- get(oris)
  for(seed_num in 1:10) { # For 10 randomizations, initial seed
    
    ####  normal randomization #### 
    rand_int_oris <- generate_random_resp_oris_normal(oris_ranges, seed = seed_num, mm = 10)
    
    
    #### randomization keeping distance to TSS #### 
    ## Initial seed for TSS randomization
    reset_i <- (seed_num-1) * 1000 # For TSS seeds 
    # Note: 
    # Function RANDOM_range_keeping_distance_to_TSS preserves relocates an origin at a distance and direction from a random gene start that
    # is equal to the distance to the closest gene start in the original set. 
    # However, after randomization, it is possible that the origin has another gene start that is closer in the genome. 
    # In such cases we will re-randomize the origin until it finds a position that preserves the minimum distance to the closest gene start. 
    # If after 1,000 randomizations it doesn't randomly find a position at the the exact distance, we increase the range of distances progressively. 
    
   
    # Assign names to each ori
    names(oris_ranges) <- paste0("ori_", 1:length(oris_ranges))
    
    oris_to_randomize_initial <- oris_ranges
    correct_oris      <- GRanges()
    
    oris_to_randomize <- oris_to_randomize_initial
    
    i <- reset_i
    while(length(correct_oris) < length(oris_to_randomize_initial)) {
      # Seed
      i <- i+1
      tss_random_oris <- RANDOM_range_keeping_distance_to_TSS(oris_to_randomize, seed = i)
      names(tss_random_oris) <- names(oris_to_randomize)
      
      # Check if distance is ok
      # Original distance to nearest
      center_of_origin <- apply(cbind(start(oris_to_randomize),
                                      end(oris_to_randomize)), 1, function(x) round(mean(x)))
      range_center <- oris_to_randomize
      start(range_center) <- center_of_origin
      end(range_center)   <- center_of_origin
      hits_original <- distanceToNearest(range_center, genes_starts)
      
      # Random distance to nearest
      center_of_origin <- apply(cbind(start(tss_random_oris),end(tss_random_oris)), 1, function(x) round(mean(x)))
      range_center <- tss_random_oris
      start(range_center) <- center_of_origin
      end(range_center)   <- center_of_origin
      hits_random <- distanceToNearest(range_center, genes_starts)
      
      dist_original <- hits_original@elementMetadata$distance*genes_starts[subjectHits(hits_original)]$Strand
      dist_random   <- hits_random@elementMetadata$distance*genes_starts[subjectHits(hits_random)]$Strand
      
      
      # To be ok, must be at the same distance - direction to TSS and not on low-mappability region
      ok     <- names(oris_to_randomize)[(dist_original == dist_random) & 
                                           overlapsAny(tss_random_oris, chromosomes_no_shade_areas_UCSC, type = "within")]
      not_ok <- names(oris_to_randomize)[!names(oris_to_randomize) %in% ok]
      
      
      correct_oris      <- append(correct_oris, tss_random_oris[ok])
      oris_to_randomize <- oris_to_randomize_initial[not_ok]
      
      percent_correct   <- length(correct_oris)/length(oris_to_randomize_initial)*100
      
      if(i%%500 == 0) {
        message(c(length(correct_oris), " of ", length(oris_to_randomize_initial), 
                  " origins successfully relocated (", round(percent_correct, 1), "%)"))
      }

      if(i%%1000 == 0) break
    }
    
    # Amplify distance after 1000 seeds 
    margin_bp <- 0 # Increase of the distance from a gene start at which the origin can be relocated
    i <- reset_i
    
    while(length(correct_oris) < length(oris_to_randomize_initial)) {
      if(length(correct_oris)/length(oris_to_randomize_initial)*100 > 99.95) { # If only 0.05% origins remain, avoid re-randomization
        tss_random_oris <- RANDOM_range_keeping_distance_to_TSS(oris_to_randomize, seed = i)
        names(tss_random_oris) <- names(oris_to_randomize)
        
        correct_oris      <- append(correct_oris, tss_random_oris)
      } else {
        if(i %% 100 == 0) margin_bp <- margin_bp + 100
        
        # Seed
        i <- i+1
        tss_random_oris <- RANDOM_range_keeping_distance_to_TSS(oris_to_randomize, seed = i)
        names(tss_random_oris) <- names(oris_to_randomize)
        
        # Check if distance is ok
        # Original distance to nearest
        center_of_origin <- apply(cbind(start(oris_to_randomize),
                                        end(oris_to_randomize)), 1, function(x) round(mean(x)))
        range_center <- oris_to_randomize
        start(range_center) <- center_of_origin
        end(range_center)   <- center_of_origin
        hits_original <- distanceToNearest(range_center, genes_starts)
        
        # Random distance to nearest
        center_of_origin <- apply(cbind(start(tss_random_oris),end(tss_random_oris)), 1, function(x) round(mean(x)))
        range_center <- tss_random_oris
        start(range_center) <- center_of_origin
        end(range_center)   <- center_of_origin
        hits_random <- distanceToNearest(range_center, genes_starts)
        
        dist_original <- hits_original@elementMetadata$distance*genes_starts[subjectHits(hits_original)]$Strand
        dist_random   <- hits_random@elementMetadata$distance*genes_starts[subjectHits(hits_random)]$Strand
        
        
        # To be ok, must be at the same distance - direction to TSS and not on low-mappability region
        ok     <- names(oris_to_randomize)[ abs(dist_original - dist_random) < margin_bp & 
                                              overlapsAny(tss_random_oris, chromosomes_no_shade_areas_UCSC, 
                                                          type = "within")]
        not_ok            <- names(oris_to_randomize)[!names(oris_to_randomize) %in% ok]
        
        # Append "ok" origins to final object
        correct_oris      <- append(correct_oris, tss_random_oris[ok])
        
        # Origins that will need to be re-relocated
        oris_to_randomize <- oris_to_randomize_initial[not_ok]

        # Correctly randomized origins, percentage
        percent_correct   <- length(correct_oris)/length(oris_to_randomize_initial)*100
        
        if(i%%500 == 0) {
          message(c(length(correct_oris), " of ", length(oris_to_randomize_initial), 
                    " origins successfully relocated (", round(percent_correct, 1), "%)"))
        }
      }
    
    }  
    
    # Write tables  
    # TSS randomizations
    df <- as.data.frame(correct_oris)[,1:3]
    write.table(df, col.names = F, row.names = F, quote = F,  sep = "\t",
                file = paste0(output_folder, "randomized_mm10", "/random_", oris ,"_TSS_mm10_", 
                              "_seed", seed_num, ".bed"))
    # Normal randomizations 
    df <- as.data.frame(rand_int_oris)[,1:3]
    write.table(df, col.names = F, row.names = F, quote = F,  sep = "\t",
                file = paste0(output_folder, "randomized_mm10", "/random_", oris ,"_NORMAL_mm10_", 
                              "_seed", seed_num, ".bed"))

    
    print(seed_num)
  }
}

