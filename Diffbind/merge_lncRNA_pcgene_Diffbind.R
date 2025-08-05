library(tidyverse)
library(openxlsx)
library(grid)
library(tidyverse)
library(tidytext)
library(rtracklayer)
getwd()
####Input assignment results for pcgenes
in.dir1 <- "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Output/DBA_edgeR_assign/narrow/"
in.dir2 <- "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Output/DBA_edgeR_assign/broad/"


####Input assignment results for lncRNA
in.dir3 <- "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Output/DBA_edgeR_assign_lncRNA/narrow/"
in.dir4 <- "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Output/DBA_edgeR_assign_lncRNA/broad/"


in.dir.list <- list(in.dir1,in.dir2,in.dir3,in.dir4)

####Output directory
out.dir1 <- "./Output/Diffbind_merged_pcgene_lncRNA_assignment/narrow/"
out.dir2 <- "./Output/Diffbind_merged_pcgene_lncRNA_assignment/broad/"

out.dir.list <- list(out.dir1,out.dir2)


#cutoff for filtering
#Format in MACS2 Output: By default, MACS2 reports -log10(pValue)/-log10(qValue), -log10(0.05)=1.3
#cutoff <- 2


for (j in 1:2){
  list.files(in.dir.list[[j]])
  list.files(in.dir.list[[j+2]])
  count.files1 <- list.files(in.dir.list[[j]])
  count.files2 <- list.files(in.dir.list[[j+2]])
  #####make sure the 8mM_statistics and UNT_statistics folders are not in the input directory
  n <- 6
  
  
  ###No BR for intersect bed files
  sample_table_fun <- function(count.files, n){
    sample_table <- data.frame(Con = rep(unique(unlist(lapply(strsplit(count.files, "_"),"[[", 1))), 3),
                               Stage = rep(rep(unique(unlist(lapply(strsplit(count.files, "_"),"[[", 2)))), 2),
                               ###The input is narrowPeak files
                               Index=grep("",count.files))
    return(sample_table)
  }
  
  sample_table <-sample_table_fun(count.files1, n)
  #sample_table <-sample_table_fun(count.files2, n) 
  
  
  in.dir_a <- in.dir.list[[j]]
  in.dir_b <- in.dir.list[[j+2]]
  ##pcgenes
  count.files_a <- count.files1
  ##lncRNA
  count.files_b <- count.files2
  cat(paste("processing", in.dir_a),"\n")
  cat(paste("processing", in.dir_b),"\n")
  out.dir <- out.dir.list[[j]]
  for (i in 1:nrow(sample_table)){
    sample_table_m <- sample_table[,c(1:2)]
    sample_name <- paste(sample_table_m[i,],collapse = "_")
    ind <- sample_table$Index[i]
    ####make sure the input file is excel(read.xlsx) rather than bed file (read.table)
    pcgene_assign <-try({read.xlsx(paste0(in.dir_a, count.files_a[ind]))}, silent = TRUE)
    lncRNA_assign <-try({read.xlsx(paste0(in.dir_b, count.files_b[ind]))}, silent = TRUE)
    if (inherits(pcgene_assign, "try-error")){next}
    
    if(j<2){
      lncRNA_assign2 <- lncRNA_assign[,c(12:18,20:24)]
      colnames(lncRNA_assign2) <- c("binding_ID","lncRNA_ID","lncRNA_strandness","lncRNA_overlapped","lncRNA_Category","lncRNA_overlapped_percentage",
                                    "lncRNA_unique_assignment","lncRNA_refgene_ID","lncRNA_refgene_strandness",
                                    "lncRNA_distance_to_refgene","lncRNA_refgene_ProductDescription","lncRNA_refgene_Symbol")
      merged_results <- merge(pcgene_assign, lncRNA_assign2, by = "binding_ID", all = TRUE)
      ######optional:further filtering steps
      #merged_results <- merged_results%>%dplyr::filter(signalValue1>cutoff&signalValue2>cutoff)
      write.xlsx(merged_results, paste0(out.dir,sample_name,"_allpeaks_assign_merged_narrow.xlsx"))
      cat(paste("processing", sample_name),"\n")
      
    }else{
      lncRNA_assign2 <- lncRNA_assign[,c(12:18)]
      colnames(lncRNA_assign2) <- c("binding_ID","lncRNA_ID","lncRNA_overlapped","lncRNA_num","lncRNA_Category","lncRNA_overlapped_percentage",
                                    "lncRNA_unique_assignment")
      merged_results <- merge(pcgene_assign, lncRNA_assign2, by = "binding_ID", all = TRUE)
      ######optional:further filtering steps
      #merged_results <- merged_results%>%dplyr::filter(signalValue1>cutoff&signalValue2>cutoff)
      write.xlsx(merged_results, paste0(out.dir,sample_name,"_allpeaks_assign_merged_broad.xlsx"))
      cat(paste("processing", sample_name),"\n")
    }
    
  }
}




