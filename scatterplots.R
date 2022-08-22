
library(MASS)
library(ggplot2)
library(viridis)
library(GenomicRanges)

setwd("RepOris/classic_efficiencies/mm10/normalized_efficiency/")

wt_MGC5 <- read.table("efficiencies__DOWNSAMPLED_MGC5_mm10.bed", header = T)
wt_SNS <- read.table("efficiencies__DOWNSAMPLED_SNS_H1_WT_II_mm10.bed", header = T)

aph_MGC4 <- read.table("efficiencies_MGC4_mm10.bed", header = T)
aph_MGC7 <- read.table("efficiencies__DOWNSAMPLED_MGC7_mm10.bed", header = T)

cdc6_JMR2 <- read.table("efficiencies__DOWNSAMPLED_JMR2_mm10.bed", header = T)
cdc6_MGC3 <- read.table("efficiencies__DOWNSAMPLED_MGC3_mm10.bed", header = T)

mean_wt  <- cbind(wt_MGC5[,1:3],  (wt_MGC5[4:7]+wt_SNS[4:7])/2)
mean_aph <- cbind(aph_MGC4[,1:3], (aph_MGC4[4:7]+aph_MGC7[4:7])/2)
mean_cdc6 <- cbind(cdc6_JMR2[,1:3], (cdc6_JMR2[4:7]+cdc6_MGC3[4:7])/2)

all_norm_corrected <- cbind(wt_MGC5[,1:3], 
                            wt_MGC5 = wt_MGC5$norm_background_corrected,
                            wt_SNS = wt_SNS$norm_background_corrected,
                            aph_MGC4 = aph_MGC4$norm_background_corrected,
                            aph_MGC7 = aph_MGC7$norm_background_corrected,
                            cdc6_JMR2 = cdc6_JMR2$norm_background_corrected,
                            cdc6_MGC3 = cdc6_MGC3$norm_background_corrected,
                            mean_wt = mean_wt$norm_background_corrected,
                            mean_aph = mean_aph$norm_background_corrected,
                            mean_cdc6 = mean_cdc6$norm_background_corrected)

all_oris_ranges <- makeGRangesFromDataFrame(all_norm_corrected)

# Check which one are responsive
aph_resp <- makeGRangesFromDataFrame(read.table("means_mm10/mean_efficiency_aph_responsive_mm10.bed", header = T))
cdc6_resp <- makeGRangesFromDataFrame(read.table("means_mm10/mean_efficiency_cdc6_responsive_mm10.bed", header = T))
comm <- makeGRangesFromDataFrame(read.table("means_mm10/mean_efficiency_constitutive_all_samples_mm10.bed", header = T))

all_norm_corrected$aph_resp <- overlapsAny(all_oris_ranges, aph_resp)
all_norm_corrected$cdc6_resp <- overlapsAny(all_oris_ranges, cdc6_resp)
all_norm_corrected$comm <- overlapsAny(all_oris_ranges, comm)


# Color points by density with ggplot2
# https://slowkow.com/notes/ggplot2-color-by-density/
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

comparisons <- list("WT replicates" = c("wt_MGC5", "wt_SNS"),
                    "APH replicates" = c("aph_MGC4", "aph_MGC7"),
                    "CDC6 replicates" = c("cdc6_MGC3", "cdc6_JMR2"),
                    "mean WT vs mean APH" = c("mean_wt", "mean_aph"),
                    "mean WT vs mean CDC6" = c("mean_wt", "mean_cdc6"),
                    "WT I vs mean APH" = c("wt_SNS", "mean_aph"),
                    "WT II vs mean APH" = c("wt_MGC5", "mean_aph"),
                    "WT I vs mean CDC6" = c("wt_SNS", "mean_cdc6"),
                    "WT II vs mean CDC6" = c("wt_MGC5", "mean_cdc6"),
                    "WT I vs APH I" = c("wt_SNS", "aph_MGC4"),
                    "WT I vs APH II" = c("wt_SNS", "aph_MGC7"),
                    "WT II vs APH I" = c("wt_MGC5", "aph_MGC4"),
                    "WT II vs APH II" = c("wt_MGC5", "aph_MGC7"),
                    "WT I vs CDC6 I" = c("wt_SNS", "cdc6_MGC3"),
                    "WT I vs CDC6 II" = c("wt_SNS", "cdc6_JMR2"),
                    "WT II vs CDC6 I" = c("wt_MGC5", "cdc6_MGC3"),
                    "WT II vs CDC6 II" = c("wt_MGC5", "cdc6_JMR2"))


for(comp in c(names(comparisons))) {
  dat <- cbind(cond1 = all_norm_corrected[, comparisons[[comp]][1]],
               cond2 = all_norm_corrected[, comparisons[[comp]][2]],
               all_norm_corrected[c("aph_resp", "cdc6_resp", "comm")])
  
  # Calculating percentage at each side of diagonal
  message("Higher in ", comparisons[[comp]][2], " than in ", comparisons[[comp]][1], ": ")
  message(round((prop.table(table(dat$cond2 > dat$cond1))*100)["TRUE"], 2), "%")
  message("*******************")
  
  # max_value <- max(c(dat$cond1, dat$cond2), na.rm = T)
  max_value <- 50
  dat$density <- get_density(dat$cond1, dat$cond2, n = 100)
  gg <- ggplot(dat) + geom_point(aes(cond1, cond2, color = density), size = 0.1) +
    theme_classic() + 
    labs(x = comparisons[[comp]][1], y = comparisons[[comp]][2], 
         title = paste("Normalized efficiencies in", comp)) + 
    xlim(c(0, max_value)) + ylim(c(0, max_value)) + coord_fixed()
  
  # Scatterplot with density scale
  pdf(paste0("with_diagonal_density_scatterplot_", gsub(" ", "_", comp),  ".pdf"),
      height = 5, width = 5)
  plot(gg + geom_abline(slope = 1, intercept = 0, col = "black"))
  
  dev.off()
  
  # Scatterplots by group with coloured comm/resp oris
  pdf(paste0("with_diagonal_highlighted_resp_or_comm_scatterplot_", gsub(" ", "_", comp),  ".pdf"),
      height = 5, width = 5)
  if(comp == "WT replicates") {
    gg <- ggplot(dat) + geom_point(aes(cond1, cond2, color = comm), size = 0.1) +
      theme_classic() + 
      labs(x = comparisons[[comp]][1], y = comparisons[[comp]][2], 
           title = paste("Normalized efficiencies in", comp)) +
      guides(colour = guide_legend(override.aes = list(size=5))) +
      scale_colour_manual(values = c("gray50", "gold")) + 
      xlim(c(0, max_value)) + ylim(c(0, max_value)) + coord_fixed()
    
    plot(gg + geom_abline(slope = 1, intercept = 0, col = "black"))
    
  }
  if(grepl("APH", comp)) {
    gg <- ggplot(dat) + geom_point(aes(cond1, cond2, color = aph_resp), size = 0.1) +
      theme_classic() + 
      labs(x = comparisons[[comp]][1], y = comparisons[[comp]][2], 
           title = paste("Normalized efficiencies in", comp)) +
      guides(colour = guide_legend(override.aes = list(size=5))) +
      scale_colour_manual(values = c("gray50", "red")) + 
      xlim(c(0, max_value)) + ylim(c(0, max_value)) + coord_fixed()
    plot(gg + geom_abline(slope = 1, intercept = 0, col = "black"))
    
  }
  if(grepl("CDC6", comp)) {
    gg <- ggplot(dat) + geom_point(aes(cond1, cond2, color = cdc6_resp), size = 0.1) +
      theme_classic() + 
      labs(x = comparisons[[comp]][1], y = comparisons[[comp]][2], 
           title = paste("Normalized efficiencies in", comp)) +
      guides(colour = guide_legend(override.aes = list(size=5))) +
      scale_colour_manual(values = c("gray50", "blue3")) + 
      xlim(c(0, max_value)) + ylim(c(0, max_value)) + coord_fixed()
    plot(gg + geom_abline(slope = 1, intercept = 0, col = "black"))
    
  }
  dev.off()

}                  
