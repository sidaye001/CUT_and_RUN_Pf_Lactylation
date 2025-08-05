library(tidyverse)
library(openxlsx)
library(grid)
library(tidyverse)
library(tidytext)
library(rtracklayer)
getwd()
####Input assignment results for pcgenes
in.dir1 <- "./Output/gene_peak_assignment_narrowPeak2/unfiltered/8mM/"
in.dir2 <- "./Output/gene_peak_assignment_narrowPeak2/unfiltered/UNT/"
in.dir3 <- "./Output/gene_peak_assignment_broadPeak2/unfiltered/8mM/"
in.dir4 <- "./Output/gene_peak_assignment_broadPeak2/unfiltered/UNT/"

####Input assignment results for lncRNA
in.dir5 <- "./Output/lncRNA_occupancy_genePeakalignment/narrow/8mM/"
in.dir6 <- "./Output/lncRNA_occupancy_genePeakalignment/narrow/UNT/"
in.dir7 <- "./Output/lncRNA_occupancy_genePeakalignment/broad/8mM/"
in.dir8 <- "./Output/lncRNA_occupancy_genePeakalignment/broad/UNT/"

in.dir.list <- list(in.dir1,in.dir2,in.dir3,in.dir4,in.dir5,in.dir6,in.dir7,in.dir8)

####Output directory
out.dir1 <- "./Output/Occupancy_merged_pcgene_lncRNA_assignment/narrow/8mM/"
out.dir2 <- "./Output/Occupancy_merged_pcgene_lncRNA_assignment/narrow/UNT/"
out.dir3 <- "./Output/Occupancy_merged_pcgene_lncRNA_assignment/broad/8mM/"
out.dir4 <- "./Output/Occupancy_merged_pcgene_lncRNA_assignment/broad/UNT/"

out.dir.list <- list(out.dir1,out.dir2,out.dir3,out.dir4)


#cutoff for filtering
#Format in MACS2 Output: By default, MACS2 reports -log10(pValue)/-log10(qValue), -log10(0.05)=1.3
#cutoff <- 2

for (j in 1:4){
  list.files(in.dir.list[[j]])
  list.files(in.dir.list[[j+4]])
  count.files1 <- list.files(in.dir.list[[j]])
  count.files2 <- list.files(in.dir.list[[j+4]])
  #####make sure the 8mM_statistics and UNT_statistics folders are not in the input directory
  n <- 12
  
  
  ###No BR for intersect bed files
  sample_table_fun <- function(count.files, n){
    sample_table <- data.frame(Con = rep(unique(unlist(lapply(strsplit(count.files, "_"),"[[", 1))), n),
                               Sample =rep(unique(unlist(lapply(strsplit(count.files, "_"),"[[", 2))), each=3),
                               Stage = rep(rep(unique(unlist(lapply(strsplit(count.files, "_"),"[[", 3)))), 4),
                               ###The input is narrowPeak files
                               Index=grep("",count.files))
    return(sample_table)
  }
  
  sample_table <-sample_table_fun(count.files1, n)
  #sample_table <-sample_table_fun(count.files2, n) 
  
  
  in.dir_a <- in.dir.list[[j]]
  in.dir_b <- in.dir.list[[j+4]]
  ##pcgenes
  count.files_a <- count.files1
  ##lncRNA
  count.files_b <- count.files2
  cat(paste("processing", in.dir_a),"\n")
  cat(paste("processing", in.dir_b),"\n")
  out.dir <- out.dir.list[[j]]
  for (i in 1:nrow(sample_table)){
    sample_table_m <- sample_table[,c(1:3)]
    sample_name <- paste(sample_table_m[i,],collapse = "_")
    ind <- sample_table$Index[i]
    ####make sure the input file is excel(read.xlsx) rather than bed file (read.table)
    pcgene_assign <-try({read.xlsx(paste0(in.dir_a, count.files_a[ind]))}, silent = TRUE)
    lncRNA_assign <-try({read.xlsx(paste0(in.dir_b, count.files_b[ind]))}, silent = TRUE)
    if (inherits(pcgene_assign, "try-error")){next}
    
    if(j<3){
      lncRNA_assign2 <- lncRNA_assign[,c(4,13:15,17:24)]
      colnames(lncRNA_assign2) <- c("peak_id","lncRNA_ID","lncRNA_overlapped","lncRNA_Category","lncRNA_overlapped_percentage",
                                    "lncRNA_unique_assignment","lncRNA_strandness","lncRNA_refgene_ID","lncRNA_refgene_strandness",
                                    "lncRNA_distance_to_refgene","lncRNA_refgene_ProductDescription","lncRNA_refgene_Symbol")
      merged_results <- merge(pcgene_assign, lncRNA_assign2, by = "peak_id", all = TRUE)
      ######optional:further filtering steps
      #merged_results <- merged_results%>%dplyr::filter(signalValue1>cutoff&signalValue2>cutoff)
      write.xlsx(merged_results, paste0(out.dir,sample_name,"_allpeaks_assign_merged.xlsx"))
      cat(paste("processing", sample_name),"\n")
      
    }else{
      lncRNA_assign2 <- lncRNA_assign[,c(4,13:16,18:19)]
      colnames(lncRNA_assign2) <- c("peak_id","lncRNA_ID","lncRNA_overlapped","lncRNA_num","lncRNA_Category","lncRNA_overlapped_percentage",
                                    "lncRNA_unique_assignment")
      merged_results <- merge(pcgene_assign, lncRNA_assign2, by = "peak_id", all = TRUE)
      ######optional:further filtering steps
      #merged_results <- merged_results%>%dplyr::filter(signalValue1>cutoff&signalValue2>cutoff)
      write.xlsx(merged_results, paste0(out.dir,sample_name,"_allpeaks_assign_merged.xlsx"))
      cat(paste("processing", sample_name),"\n")
    }
    
  }
}




