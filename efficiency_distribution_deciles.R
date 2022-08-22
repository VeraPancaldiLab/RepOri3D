# Figures 2EF and S4AB

library(GenomicRanges)
library(gtools)

setwd("RepOris/classic_efficiencies/mm10/means_mm10/")

wt = read.table("mean_efficiency_wt_mm10.bed", stringsAsFactors = F, header = T)
aph = read.table("mean_efficiency_aph_mm10.bed", stringsAsFactors = F, header = T)
cdc6 = read.table("mean_efficiency_cdc6_mm10.bed", stringsAsFactors = F, header = T)

wt <- makeGRangesFromDataFrame(wt, keep.extra.columns = T)
aph <- makeGRangesFromDataFrame(aph, keep.extra.columns = T)
cdc6 <- makeGRangesFromDataFrame(cdc6, keep.extra.columns = T)

all_merged <- reduce(append(cdc6, append(wt,aph)))


# Efficiency calculated for all merged oris on each sample using Efficiencies.Rmd code
mm10_files <- list.files("../figure2B/normalized_efficiency")[grepl("mm10", list.files("../figure2B/normalized_efficiency"))]


mean_count_corrected <- c()


# Calculating
for(i in 1:length(mm10_files)) {
  x <- read.table(paste0("../figure2B/normalized_efficiency/", mm10_files[i]),
                  header = T, stringsAsFactors = F)
  if(i == 1) {
    mean_count_corrected <- x
  } else {
    mean_count_corrected[, 4:7] <- mean_count_corrected[, 4:7] + x[, 4:7]
  }
}
mean_count_corrected[, 4:7] <- mean_count_corrected[, 4:7]/6

all_merged <- makeGRangesFromDataFrame(mean_count_corrected, keep.extra.columns = T)


all_merged$mean_count_corrected_decile <- quantcut(all_merged$norm_background_corrected, 10)


COMM <- makeGRangesFromDataFrame(read.table("mean_efficiency_constitutive_mm10.bed", stringsAsFactors = F, header = T), keep.extra.columns = F)
APH_R <- makeGRangesFromDataFrame(read.table("mean_efficiency_aph_responsive_mm10.bed", stringsAsFactors = F, header = T), keep.extra.columns = F)
CDC6_R <- makeGRangesFromDataFrame(read.table("mean_efficiency_cdc6_responsive_mm10.bed", stringsAsFactors = F, header = T), keep.extra.columns = F)
APH_CDC6_R <- makeGRangesFromDataFrame(read.table("mean_efficiency_both_aph_cdc6_mm10.bed", stringsAsFactors = F, header = T), keep.extra.columns = F)


x <- c()
for(decile_num in 1:10) {
  subset_merged <- all_merged[all_merged$mean_count_corrected_decile == rev(levels(all_merged$mean_count_corrected_decile))[decile_num]]
  x <- rbind(x, c("COMM" = prop.table(table(overlapsAny(COMM, subset_merged)))["TRUE"]*100,
                  "APH-R" = prop.table(table(overlapsAny(APH_R, subset_merged)))["TRUE"]*100,
                  "CDC6-R" = prop.table(table(overlapsAny(CDC6_R, subset_merged)))["TRUE"]*100,
                  "APH+CDC6-R" = prop.table(table(overlapsAny(APH_CDC6_R, subset_merged)))["TRUE"]*100))

}
rownames(x) <- 1:10

colnames(x) <- gsub(".TRUE", "", colnames(x))

pdf("Figures_2EF.pdf", height = 6, width = 13)
par(mar = c(4,4,3,4), mfrow = c(1,2))


boxplot(rev(by(all_merged$norm_background_corrected, 
               all_merged$mean_count_corrected_decile, list)), 
        outline = F, 
        ylab = "Average efficiency",
        names = c(1:10), 
        # at = n[2, ]+0.5,
        # xaxt = F,
        # boxwex = 4,
        ylim = c(0, 20),
        xlab = "Origin efficiency quantiles", 
        las = 1, 
        col = NA,
        add = F)


barplot(t(x), beside = T, legend = T, col = c("darkgreen", "goldenrod", "lightblue", "green"), 
        args.legend = list(bty = "n"), 
        # main = "Quantiles calculated without merging origins",
        cex.main = 0.9, 
        las = 1,
        xlab = "Origin efficiency quantiles", ylab = "% of origins in each quantile", 
        ylim = c(0,30))
axis(4, at = seq(0, 30, length.out =7),
     las = 1, 
     labels = seq(0, 30, length.out = 7))

mtext(text = "% of origins in each quantile", side = 4, line = 2.5)




dev.off()


##### Using read density (not normalised by ori size)

all_merged$decile_background_corrected <- quantcut(all_merged$background_corrected, 10)

x <- c()
for(decile_num in 1:10) {
  subset_merged <- all_merged[all_merged$decile_background_corrected == rev(levels(all_merged$decile_background_corrected))[decile_num]]
  x <- rbind(x, c("COMM" = prop.table(table(overlapsAny(COMM, subset_merged)))["TRUE"]*100,
                  "APH-R" = prop.table(table(overlapsAny(APH_R, subset_merged)))["TRUE"]*100,
                  "CDC6-R" = prop.table(table(overlapsAny(CDC6_R, subset_merged)))["TRUE"]*100,
                  "APH+CDC6-R" = prop.table(table(overlapsAny(APH_CDC6_R, subset_merged)))["TRUE"]*100))
  
}
rownames(x) <- 1:10
colnames(x) <- gsub(".TRUE", "", colnames(x))



pdf("Figures_Supp4_AB.pdf", height = 6, width = 13)
par(mar = c(4,4,3,4), mfrow = c(1,2))

boxplot(rev(by(all_merged$background_corrected/1000, 
               all_merged$decile_background_corrected, list)), 
        outline = F, 
        ylab = "Efficiency (per base coverage/1000)",
        names = c(1:10), 
        # at = n[2, ]+0.5,
        # xaxt = F,
        # boxwex = 4,
        xlab = "Origin efficiency quantiles", 
        las = 1, 
        col = NA,
        add = F)


barplot(t(x), beside = T, legend = T, col = c("darkgreen", "goldenrod", "lightblue", "green"), 
        args.legend = list(bty = "n"), 
        # main = "Quantiles calculated without merging origins",
        cex.main = 0.9, 
        las = 1,
        xlab = "Origin efficiency quantiles", ylab = "% of origins in each quantile", 
        ylim = c(0,35))

# mtext(text = "Average efficiency", side = 4, line = 2.5)

dev.off()







