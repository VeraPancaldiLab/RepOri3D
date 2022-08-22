# Provide a table with origin number, origin size, location and replication timing regions of 
# CDC6-R, Aph-R, and both CDC6+Aph-R origins, together with COMM ones as a control. 
# Overlap these different groups of origins with previously mapped mESC origins (Cayrou et al., 2015). 
# Do we detect any bias between origin size and efficiency?
library(liftOver)

# Check which one are responsive
p <- "RepOris/classic_efficiencies/mm10/means_mm10/"

aph_all   <- read.table(paste0(p, "mean_efficiency_aph_mm10.bed"), header = T)
cdc6_all  <- read.table(paste0(p, "mean_efficiency_cdc6_mm10.bed"), header = T)
wt_all    <- read.table(paste0(p, "mean_efficiency_wt_mm10.bed"), header = T)
aph_resp  <- read.table(paste0(p, "mean_efficiency_aph_responsive_mm10.bed"), header = T)
cdc6_resp <- read.table(paste0(p, "mean_efficiency_cdc6_responsive_mm10.bed"), header = T)
# both_resp <- read.table(paste0(p, "mean_efficiency_both_aph_cdc6_mm10.bed"), header = T)
comm      <- read.table(paste0(p, "mean_efficiency_constitutive_all_samples_mm10.bed"), header = T)

# cdc6+aph responsive origins
aph_and_cdc6_resp <- read.table(file = "RepOris/classic_origins/mm10/CDC6andAPHresponsive_classic_efficiency.xls",
                                stringsAsFactors = F)
colnames(aph_and_cdc6_resp)[1:3] <- c("chr", "start", "end")
aph_and_cdc6_resp_ranges <- makeGRangesFromDataFrame(aph_and_cdc6_resp) 

aph_all_ranges <- makeGRangesFromDataFrame(aph_all, keep.extra.columns = T)
cdc6_all_ranges <- makeGRangesFromDataFrame(cdc6_all, keep.extra.columns = T)

aph_shared_resp <- aph_all_ranges[overlapsAny(aph_all_ranges,aph_and_cdc6_resp_ranges )]
cdc6_shared_resp <- cdc6_all_ranges[overlapsAny(cdc6_all_ranges,aph_and_cdc6_resp_ranges )]


# origin number
table_for_NAR <- data.frame(total_oris = c(wt_all = nrow(wt_all),
                                           aph_all = nrow(aph_all),
                                           cdc6_all = nrow(cdc6_all),
                                           comm = nrow(comm),
                                           aph_resp = nrow(aph_resp),
                                           cdc6_resp = nrow(cdc6_resp),
                                           aph_and_cdc6_resp = nrow(aph_and_cdc6_resp)))

# origin size, 
table_for_NAR$median_size <- c(wt_all = median(wt_all$end-wt_all$start+1),
                               aph_all = median(aph_all$end-aph_all$start+1),
                               cdc6_all = median(cdc6_all$end-cdc6_all$start+1),
                               comm = median(comm$end-comm$start+1),
                               aph_resp = median(aph_resp$end-aph_resp$start+1),
                               cdc6_resp = median(cdc6_resp$end-cdc6_resp$start+1),
                               aph_and_cdc6_resp = median(c(width(aph_shared_resp), width(cdc6_shared_resp))))
table_for_NAR$median_size <- round(table_for_NAR$median_size, digits = 2)

table_for_NAR$mean_size <- c(wt_all = mean(wt_all$end-wt_all$start+1),
                             aph_all = mean(aph_all$end-aph_all$start+1),
                             cdc6_all = mean(cdc6_all$end-cdc6_all$start+1),
                             comm = mean(comm$end-comm$start+1),
                             aph_resp = mean(aph_resp$end-aph_resp$start+1),
                             cdc6_resp = mean(cdc6_resp$end-cdc6_resp$start+1),
                             aph_and_cdc6_resp = mean(c(width(aph_shared_resp), width(cdc6_shared_resp))))
table_for_NAR$mean_size <- round(table_for_NAR$mean_size, digits = 2)



# replication timing regions
table_for_NAR$median_RT <- NA
table_for_NAR$mean_RT <- NA

RT <- makeGRangesFromDataFrame(read.table("RepOris/NAR_revision/RT_data/RT_46C_ESC_Int62150809_mm10.bedgraph", 
                                          col.names = c("chr", "start", "end", "RT")),
                               keep.extra.columns = T)

for(group in rownames(table_for_NAR)) {
  rt_per_group <- RT$RT[overlapsAny(RT, makeGRangesFromDataFrame(get(group)))]
  table_for_NAR[group, "median_RT"] <- round(median(rt_per_group), 2)
  table_for_NAR[group, "mean_RT"]   <- round(mean(rt_per_group), 2)
  
}



# Overlap with Cayrou et al., 2015
# Downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68347
cayrou_sites <- makeGRangesFromDataFrame(read.table("cayrou2015/GSE68347_Initiation_Sites.bedGraph", col.names = c("chr", "start", "end", "V4")))
cayrou_zones <- makeGRangesFromDataFrame(read.table("cayrou2015/GSE68347_Initiation_Zones.bed", col.names = c("chr", "start", "end")))

ch = import.chain("RepOris/mm9ToMm10.over.chain")
cayrou_sites <- liftOver(cayrou_sites, ch)
cayrou_zones <- liftOver(cayrou_zones, ch)



table_for_NAR$percent_overlap_cayrou_sites <- NA
table_for_NAR$percent_overlap_cayrou_zones <- NA

for(group in rownames(table_for_NAR)) {
  # percent_overlap_sites <- prop.table(table(overlapsAny(makeGRangesFromDataFrame(get(group)), cayrou_sites))["TRUE"])*100
  # percent_overlap_zones <- prop.table(table(overlapsAny(makeGRangesFromDataFrame(get(group)), cayrou_zones))["TRUE"])*100
  
  table_for_NAR[group, "percent_overlap_cayrou_sites"] <- round((prop.table(table(overlapsAny(makeGRangesFromDataFrame(get(group)), cayrou_sites)))*100)["TRUE"], 2)
  table_for_NAR[group, "percent_overlap_cayrou_zones"] <- round((prop.table(table(overlapsAny(makeGRangesFromDataFrame(get(group)), cayrou_zones)))*100)["TRUE"], 2)
  
}

table_for_NAR <- cbind(ori_group = rownames(table_for_NAR), table_for_NAR)



write.table(table_for_NAR, file = "table_size_RT_overlap_cayrou.tsv", 
            row.names = F, sep = "\t", quote = F)



# Ori size distribution
hist(width(makeGRangesFromDataFrame(wt_all)), breaks = 50, xlab = "WT ori size (bp)", main = "Ori size distribution in WT", col = "gray70", border = "gray70")
abline(v = mean(width(makeGRangesFromDataFrame(wt_all))), col = "blue")
abline(v = median(width(makeGRangesFromDataFrame(wt_all))), col = "red")

par(mfrow = c(3,1))
# Ori RT distribution
RT_comm <- RT$RT[overlapsAny(RT, makeGRangesFromDataFrame(comm))]
hist(RT_comm, breaks = 50, xlab = "Replication time (late --> early)", main = "RT distribution in COMM", col = "gray70", border = "gray70")
abline(v = mean(RT_comm), col = "blue")
abline(v = median(RT_comm), col = "red")

RT_aph <- RT$RT[overlapsAny(RT, makeGRangesFromDataFrame(aph_resp))]
hist(RT_aph, breaks = 50, xlab = "Replication time (late --> early)", main = "RT distribution in APH-resp", col = "gray70", border = "gray70")
abline(v = mean(RT_aph), col = "blue")
abline(v = median(RT_aph), col = "red")

RT_cdc6 <- RT$RT[overlapsAny(RT, makeGRangesFromDataFrame(cdc6_resp))]
hist(RT_cdc6, breaks = 50, xlab = "Replication time (late --> early)", main = "RT distribution in CDC6-resp", col = "gray70", border = "gray70")
abline(v = mean(RT_cdc6), col = "blue")
abline(v = median(RT_cdc6), col = "red")

