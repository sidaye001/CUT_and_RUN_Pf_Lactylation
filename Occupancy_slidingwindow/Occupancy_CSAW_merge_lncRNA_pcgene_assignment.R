library(tidyverse)
library(openxlsx)
library(grid)
library(tidyverse)
library(tidytext)
library(rtracklayer)
setwd("./")
####Input assignment results for pcgenes
in.dir1 <- "./Output/Occupancy_slidewindow_assign_final/pcgenes/"




####Input assignment results for lncRNA
in.dir2 <- "./Output/Occupancy_slidewindow_assign_final/lncRNA/"



in.dir.list <- list(in.dir1,in.dir2)

####Output directory
out.dir <- "./Output/Occupancy_slidewindow_assign_final/pcgenes_lncRNA_merge/"

n <- 24
#cutoff for filtering
#Format in MACS2 Output: By default, MACS2 reports -log10(pValue)/-log10(qValue), -log10(0.05)=1.3
#cutoff <- 2
count.files1 <- list.files(in.dir1)
count.files2 <- list.files(in.dir2)
for (i in 1:n){
  ####make sure the input file is excel(read.xlsx) rather than bed file (read.table)
  pcgene_assign <-try({read.xlsx(paste0(in.dir1, count.files1[i]))}, silent = TRUE)
  lncRNA_assign <-try({read.xlsx(paste0(in.dir2, count.files2[i]))}, silent = TRUE)
  sample_name <- sub("_differentialbinding.*", "", count.files1[i])
  if (inherits(pcgene_assign, "try-error")){next}
  
  lncRNA_assign2 <- lncRNA_assign[,c(12:16, 18:19)]
  colnames(lncRNA_assign2) <- c("binding_ID","lncRNA_ID","lncRNA_overlapped","lncRNA_num","lncRNA_Category","lncRNA_overlapped_percentage",
                                "lncRNA_unique_assignment")
  merged_results <- merge(pcgene_assign, lncRNA_assign2, by = "binding_ID", all = TRUE)
  ######optional:further filtering steps
  #merged_results <- merged_results%>%dplyr::filter(signalValue1>cutoff&signalValue2>cutoff)
  write.xlsx(merged_results, paste0(out.dir,sample_name,"_differentialbinding_assigned_merged_pcgenes_lncRNA.xlsx"))
  cat(paste("processing", sample_name),"\n")
}





