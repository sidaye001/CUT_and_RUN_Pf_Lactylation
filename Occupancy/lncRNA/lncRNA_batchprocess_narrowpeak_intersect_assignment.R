library(tidyverse)
library(bedtoolsr)
library(openxlsx)
library(grid)
library(matrixStats)
library(tidyverse)
library(tidytext)
library(RColorBrewer)
library(parallel)
library(ComplexHeatmap)
library(circlize)
library(doParallel)
library(edgeR)
# Load the rtracklayer package
library(rtracklayer)
setwd("/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/")
##This script is only for intersect.bed files gene peak assignment for lncRNA in Pf
##gene peak assignment for lncRNA only consider those peaks can be assigned to lncRNAs, lncRNAs are treated as exons with no introns
##gene peak assignment for lncRNA did not consider the intergenic peaks which can not be assigned

#Tutorial of intalling bedtoolsr :http://phanstiel-lab.med.unc.edu/bedtoolsr-install.html
#step 1: install devtools if not already installed
#install.packages("devtools")

#step 2: 
#devtools::install_github("PhanstielLab/bedtoolsr")

#step3: brew install bedtools

#step4: adding bedtools to the PATH
#1):determe which shell you are using
#echo $SHELL 
#it is either '/bin/bash' or '/bin/zsh'
#2):which bedtools #finding the path to bedtools
#3):emacs ~/.zshrc
#4):adding export PATH="/opt/homebrew/bin/bedtools:$PATH" to the end of script
#5):source ~/.zshrc #apply the change
#6):double-check
#echo $PATH
#7):manually specify the path to use bedtools from in R with:
#options(bedtools.path = "/path/to/bedtools")
#specifically, in my mac: options(bedtools.path = "/opt/homebrew/bin/")

######Broad peak file has only 9 columns while narrow peak has 10 columns(10th column shows the center of the peak)

######every time run bedtoolsr in Rstudio, need to run : options(bedtools.path = "/opt/homebrew/bin/")
options(bedtools.path = "/opt/homebrew/bin/")

#####Notice: please note that there are alternative splicing and multiple UTR for one gene exists#####
#####Notice: please note that first exon include 5UTR#####

####Step1: to get TSS loci and promoter regions for each gene####
transcript_bed <- read.table("./Input/Pf_lncRNA/Pf_lncRNA_merged.bed")
##To get TSS loci
transcript_bed2 <- transcript_bed%>%mutate(TSS=ifelse(V6=="+",V2,V3))

#######To specify the promoter region######
Promoter_range <- 1000
Promoter_bed <- transcript_bed2 %>% transmute(V1=V1,TSS=ifelse(V6=="+",V2,V3),Promoter_end=ifelse(V6=="+",TSS-Promoter_range,TSS+Promoter_range),
                                             V4=V4, V5=V5, V6=V6)
###To rearrange the start and end site
Promoter_bed2 <- Promoter_bed%>%transmute(V1=V1,V2=ifelse(V6=="+",Promoter_end,TSS),V3=ifelse(V6=="+",TSS,Promoter_end),
                                          V4=V4, V5=V5, V6=V6)
###Make those negative starting site as 0s
Promoter_bed2$V2 <- ifelse(Promoter_bed2$V2<0,0,Promoter_bed2$V2)

dim(transcript_bed)
dim(Promoter_bed2)



## using qval = 0.01 and gene size 69350000 bp in macs2 calling 
## using pval = 0.05 and gene size 27000000 bp in macs2 calling, 202308
## narrow peaks called by macs2
#narrow.peaks.file <-  "../cut_and_run_pk/Sample6-MK-Cut-20220903_S6.vs.Sample3-MK-Cut-20220903_S3_peaks.narrowPeak"

narrow.peaks.file1<-  "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Input/intersect_bed2/8mM_narrowpeak_p005/"
narrow.peaks.file2<-  "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Input/intersect_bed2/UNT_narrowpeak_p005/"
list.files(narrow.peaks.file1)
list.files(narrow.peaks.file2)
count.files1 <- list.files(narrow.peaks.file1)
count.files2 <- list.files(narrow.peaks.file2)

n <- 12

###No BR(Bioreplicates) for intersect bed files
sample_table_fun <- function(count.files, n){
  sample_table <- data.frame(Con = rep(unique(unlist(lapply(strsplit(count.files, "_"),"[[", 1))), n),
                             Sample =rep(unique(unlist(lapply(strsplit(count.files, "_"),"[[", 2))), each=3),
                             Stage = rep(rep(unique(unlist(lapply(strsplit(count.files, "_"),"[[", 3)))), 4),
                             ###The input is narrowPeak files
                             Index=grep("",count.files))
  return(sample_table)
}

sample_table1 <-sample_table_fun(count.files1, n)
sample_table2 <-sample_table_fun(count.files2, n) 


in.dir <- narrow.peaks.file1
in.dir <- narrow.peaks.file2
out.dir <- "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Output/lncRNA_occupancy_genePeakalignment/narrow/8mM/"
out.dir <- "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Output/lncRNA_occupancy_genePeakalignment/narrow/UNT/"
sample_table <- sample_table1
sample_table <- sample_table2
count.files_input <- count.files1
count.files_input <- count.files2

# Define a function to find the intersection of two ranges
range_intersection <- function(chr, start1, end1, start2, end2) {
  intersect_start <- pmax(start1, start2)
  intersect_end <- pmin(end1, end2)
  # Ensure the intersection is valid
  valid <- intersect_start <= intersect_end
  intersect_start[!valid] <- NA
  intersect_end[!valid] <- NA
  return(data.frame(chr=chr, start = intersect_start, end = intersect_end))
}

#####Basic logic for assigning the peaks to genes and different locations
#1.There are six categories for assignment: intergenic, exonic, intronic, 5'UTR, 3'UTR and promoter regions(1000bp upstream TSS)
#2.Peaks locating at exonic, intronic, 5'UTR, 3'UTR and promoter will be assigned to a gene, peaks in intergenic will not be assigned to any gene
#3.Each category will be assigned separately, since epigenetic fators do not uniquely function on one gene
#4.exon_bed is for exonic regions, intron_bed is for intronic regions, Promoter_bed is for promoter regions, utr_bed/utr5_bed/utr3_bed is for 5'UTR or 3'UTR regions
#5.TSS loci is in transcript_bed
#6. There are alternative splicing, so multiple UTR exists

statistics_table1 <- sample_table[,c(1:3)]
statistics_table2 <- data.frame(lncRNA=rep(NA,12),
                                Promoter=rep(NA,12))

statistics_table <- cbind(statistics_table1,statistics_table2)
statistics_table <- cbind(statistics_table,statistics_table2)
statistics_table_unique <- statistics_table

exon_bed <- transcript_bed

lncRNA_pcgene_assign <- read.table("./Input/Pf_lncRNA/lncRNA_closest_genes.bed")
###Remove those distance to protein coding genes > 2000bp, V14 is absolute distance
lncRNA_pcgene_assign <- lncRNA_pcgene_assign%>%dplyr::filter(V14<=2000)
lncRNA_pcgene_assign2 <- lncRNA_pcgene_assign[,c(4,6,10,12,14)]
colnames(lncRNA_pcgene_assign2) <- c("lncRNA_ID","lncRNA_strandness","refgene_ID","refgene_strandness","distance_to_refgene")

product <- read.csv("./Input/5720_total_Pf_product_description.csv")
product2 <- product[,c(1,3,4)]
colnames(product2)[1] <- "refgene_ID"
lncRNA_pcgene_assign2 <- left_join(lncRNA_pcgene_assign2, product2,by="refgene_ID")

for (i in 1:nrow(sample_table)){
  sample_table_m <- sample_table[,c(1:3)]
  sample_name <- paste(sample_table_m[i,],collapse = "_")
  ind <- sample_table$Index[i]
  narrow.peaks <-try({read.table(paste0(in.dir, count.files_input[ind]),header = F, sep = '\t', quote = NULL)}, silent = TRUE)
  #V1,V2,V3 come from rep1, V11,V12,V13 come from rep2, V7/V17 is signalValue, V8/V18 is -log10(pValue), V9/V19 is -log10(qValue)
  #V5/V15 is score range from 0-1000 to show the significance of peaks, V10/V20 is the loci of peak summit in narrowpeak(broad peaks has no peak summit)
  #V21 is the intersection length of two reps, TSS loci is stored in transcript_bed TSS column
  if (inherits(narrow.peaks, "try-error")){next}
  #narrow.peaks <- read.table(paste0(in.dir, count.files1[ind]),header = F, sep = '\t', quote = NULL)
  #narrow.peaks <- read.table(paste0(in.dir, count.files2[ind]),header = F, sep = '\t', quote = NULL)
  
  narrow.peaks$peak_id1 <- paste(narrow.peaks$V1 , paste(narrow.peaks$V2, narrow.peaks$V3, sep= "-"), sep = ":")
  narrow.peaks$peak_id2 <- paste(narrow.peaks$V11 , paste(narrow.peaks$V12, narrow.peaks$V13, sep= "-"), sep = ":")
  intersect_bed <-range_intersection(narrow.peaks$V1, narrow.peaks$V2, narrow.peaks$V3, narrow.peaks$V12, narrow.peaks$V13)
  intersect_bed$peak_id <- paste(intersect_bed$chr , paste(intersect_bed$start, intersect_bed$end, sep= "-"), sep = ":")
  
  ###V21 in narrow.pkeas is the range of intersection region of two bio-replicates(width of intersected regions)
  narrow.peaks <- cbind(intersect_bed,narrow.peaks)
  #sort the table by intersected regions 
  peaks.all.sort <- narrow.peaks %>% arrange(V1, V2, V3)
  
  ##Overlap peaks with different categories,wo=T means writing the original A and B entries plus the number of base pairs of overlap between the two features
  ###!!!V34 is the width of overlaps between each category with intersected regions###
  exon.peaks <- bedtoolsr::bt.intersect(a = peaks.all.sort, b = exon_bed, wo = T)
  if (nrow(exon.peaks) == 0) {
    exon.peaks <- data.frame(Category = character())
  } else {
    exon.peaks$Category <- "lncRNA"
  }

  promoter.peaks <- bedtoolsr::bt.intersect(a = peaks.all.sort, b = Promoter_bed2, wo = T)
  if (nrow(promoter.peaks) == 0) {
    promoter.peaks<- data.frame(Category = character())
  } else {
    promoter.peaks$Category <- "lncRNA_Promoter"
  }
  
  ####No need to assign peaks to intergenic regions for lncRNA
  lncRNA.peaks <- rbind(exon.peaks,promoter.peaks)
  ###To extract necessary peaks info, the overlapped regions with intergenic should always be 1=100%
  lncRNA.peaks2 <- lncRNA.peaks %>%
    transmute(Chr = V1, str = V2, end = V3, peak_id = V4, signalValue1=V11, pValue1=V12, qValue1=V13, peak_id1=V26,
              signalValue2=V21, pValue2=V22, qValue2=V23, peak_id2=V27,GeneID = V31, overlapped=V34, Category = Category)
  
  all.peaks <-lncRNA.peaks2
  
  #####To calculate the percentage of overlapped regions over intersected regions#####
  all.peaks$intersected_peak_width <-all.peaks$end-all.peaks$str
  all.peaks$overlapped_percentage <- all.peaks$overlapped/all.peaks$intersected_peak_width
  all.peaks <- all.peaks%>%group_by(peak_id)%>% mutate(unique_assignment = Category[which.max(overlapped_percentage)]) %>%ungroup()
  
  colnames(all.peaks)[grep("GeneID",colnames(all.peaks))] <- "lncRNA_ID"
  ####Add lncRNA to pc genes annotation
  all.peaks <- left_join(all.peaks,lncRNA_pcgene_assign2,by="lncRNA_ID")

  ####To calculate the statistics of peaks assignment separately
  sta_allpeaks <- as.data.frame(table(all.peaks$Category))
  sta_allpeaks$sum <- sum(sta_allpeaks$Freq)
  sta_allpeaks$Prop <- (sta_allpeaks$Freq/sta_allpeaks$sum)*100
  if(nrow(sta_allpeaks)==2){
    
  }else{
    df <- data.frame(Var1=colnames(statistics_table)[4:5])
    df2 <- left_join(df,sta_allpeaks, by="Var1")
    df2 <- df2 %>%
      mutate_all(~replace(., is.na(.), 0))
    sta_allpeaks <- df2
  }
  statistics_table[i,4:5] <- sta_allpeaks$Freq
  statistics_table[i,6:7] <- sta_allpeaks$Prop
  
  ####To calculate the statistics of peaks assignment uniquely
  sta_allpeaks2 <- as.data.frame(table(all.peaks$unique_assignment))
  sta_allpeaks2$sum <- sum(sta_allpeaks2$Freq)
  sta_allpeaks2$Prop <- (sta_allpeaks2$Freq/sta_allpeaks2$sum)*100
  if(nrow(sta_allpeaks2)==2){
    
  }else{
    df <- data.frame(Var1=colnames(statistics_table_unique)[4:5])
    df2 <- left_join(df,sta_allpeaks2, by="Var1")
    df2 <- df2 %>%
      mutate_all(~replace(., is.na(.), 0))
    sta_allpeaks2 <- df2
  }
  statistics_table_unique[i,4:5] <- sta_allpeaks2$Freq
  statistics_table_unique[i,6:7] <- sta_allpeaks2$Prop
  
  write.xlsx(all.peaks, paste0(out.dir,sample_name,"_allpeaks_assign_lncRNA.xlsx"))
  cat(paste("processing", sample_name),"\n")
}

out.dir.statistics <- "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Output/lncRNA_occupancy_genePeakalignment/narrow/UNT/"
write.xlsx(statistics_table, paste0(out.dir.statistics,"UNT_statistics/UNT_allpeaks_assign_statistics_separate_lncRNA.xlsx"))
write.xlsx(statistics_table_unique, paste0(out.dir.statistics,"UNT_statistics/UNT_allpeaks_assign_statistics_unique_lncRNA.xlsx"))

out.dir.statistics <- "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Output/lncRNA_occupancy_genePeakalignment/narrow/8mM/"
write.xlsx(statistics_table, paste0(out.dir.statistics,"8mM_statistics/8mM_allpeaks_assign_statistics_separate_lncRNA.xlsx"))
write.xlsx(statistics_table_unique, paste0(out.dir.statistics,"8mM_statistics/8mM_allpeaks_assign_statistics_unique_lncRNA.xlsx"))
