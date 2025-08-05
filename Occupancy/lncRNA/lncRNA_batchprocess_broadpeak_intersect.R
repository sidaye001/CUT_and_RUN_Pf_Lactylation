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

######Notice: the peaks may overlap with multiple exons and introns

######every time run bedtoolsr in Rstudio, need to run : options(bedtools.path = "/opt/homebrew/bin/")
options(bedtools.path = "/opt/homebrew/bin/")
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

broad.peaks.file1<-  "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Input/intersect_bed2/8mM_broadpeak_p001/"
broad.peaks.file2<-  "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Input/intersect_bed2/UNT_broadpeak_p001/"
list.files(broad.peaks.file1)
list.files(broad.peaks.file2)
count.files1 <- list.files(broad.peaks.file1)
count.files2 <- list.files(broad.peaks.file2)

n <- 12

####This function needs to be changed if it is .narrowpeak files
###No BR for intersect bed files
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


in.dir <- broad.peaks.file1
#in.dir <- broad.peaks.file2
out.dir <- "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Output/lncRNA_occupancy_genePeakalignment/broad/8mM/"
#out.dir <- "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Output/lncRNA_occupancy_genePeakalignment/broad/UNT/"
sample_table <- sample_table1
#sample_table <- sample_table2
count.files_input <- count.files1
#count.files_input <- count.files2

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

exon_bed <- transcript_bed
for (i in 1:nrow(sample_table)){
  sample_table_m <- sample_table[,c(1:3)]
  sample_name <- paste(sample_table_m[i,],collapse = "_")
  ind <- sample_table$Index[i]
  broad.peaks <-try({read.table(paste0(in.dir, count.files_input[ind]),header = F, sep = '\t', quote = NULL)}, silent = TRUE)
  #V1,V2,V3 come from rep1, V10,V11,V12 come from rep2, V7/V16 is signalValue, V8/V17 is -log10(pValue), V9/V18 is -log10(qValue)
  #V5/V14 is score range from 0-1000 to show the significance of peaks, broad peaks has no peak summit
  #V20 is the intersection length of two reps, TSS loci is stored in transcript_bed TSS column
  if (inherits(broad.peaks, "try-error")){next}
  #broad.peaks <- read.table(paste0(in.dir, count.files1[ind]),header = F, sep = '\t', quote = NULL)
  #broad.peaks <- read.table(paste0(in.dir, count.files2[ind]),header = F, sep = '\t', quote = NULL)
  broad.peaks$peak_id1 <- paste(broad.peaks$V1 , paste(broad.peaks$V2, broad.peaks$V3, sep= "-"), sep = ":")
  broad.peaks$peak_id2 <- paste(broad.peaks$V10 , paste(broad.peaks$V11, broad.peaks$V12, sep= "-"), sep = ":")
  intersect_bed <-range_intersection(broad.peaks$V1, broad.peaks$V2, broad.peaks$V3, broad.peaks$V11, broad.peaks$V12)
  intersect_bed$peak_id <- paste(intersect_bed$chr , paste(intersect_bed$start, intersect_bed$end, sep= "-"), sep = ":")
  
  
  broad.peaks <- cbind(intersect_bed,broad.peaks)
  #sort the table by intersected regions 
  peaks.all.sort <- broad.peaks %>% arrange(V1, V2, V3)
  dim(peaks.all.sort)
  ## Get the transcripts
  #gtf.trans <- gtf %>% dplyr::filter(V3 == 'transcript')
  #parse.str <- strsplit(gtf.trans$V9, split = ' ')
  #inds <- unlist(lapply(parse.str , function(x) which(grepl("gene_id", x)) + 1))
  #gtf.trans$V9 <- gsub(";", "", unlist(lapply(1:length(inds), function(i) parse.str[[i]][inds[[i]]])))
  
  #gtf.trans$V9 <- gsub("\"", "", gtf.trans$V9)
  
  
  ## overlap the genes with peaks
  genic.peaks <- bedtoolsr::bt.intersect(a = peaks.all.sort, b = transcript_bed, wo = T)
  
  
  ##To extract the genes info covered by each peak
  #V32 is the TSS loci#
  gene.peaks <- genic.peaks%>%group_by(V4)%>% summarise(
    V4=list(V4),
    covered_genes=list(V29),
    num=n()
  )%>%as.data.frame()
  
  gene.peaks$V4 <- unlist(lapply(gene.peaks$V4, unique))
  
  genic.peaks_m <- left_join(genic.peaks,gene.peaks, by="V4")
  
  ######To separate the peaks into two categories: covered 1 genes only and >=2 genes
  genic.peaks2 <- genic.peaks_m%>%dplyr::filter(num>=2)
  genic.peaks1 <- genic.peaks_m%>%dplyr::filter(num==1)
  
  #####No need to process the peaks covered more than one gene
  genic.peaks2_final <- genic.peaks2%>%transmute(Chr = V1, str = V2, end = V3, peak_id = V4, signalValue1=V11, pValue1=V12, qValue1=V13, peak_id1=V24,
                                                 signalValue2=V20, pValue2=V21, qValue2=V22, peak_id2=V25,GeneID = covered_genes, overlapped=V3-V2,num=num,Category = 'Covered more than one lncRNA')
  
  #####To calculate the percentage of overlapped regions over intersected regions#####
  genic.peaks2_final$intersected_peak_width <-genic.peaks2_final$end-genic.peaks2_final$str
  genic.peaks2_final$overlapped_percentage <- genic.peaks2_final$overlapped/genic.peaks2_final$intersected_peak_width
  genic.peaks2_final <- genic.peaks2_final%>%group_by(peak_id)%>% mutate(unique_assignment = Category[which.max(overlapped_percentage)]) %>%ungroup()
  # Properly convert list elements to character strings
  genic.peaks2_final$GeneID <- sapply(genic.peaks2_final$GeneID, function(x) paste(x, collapse = ","))
  
  #####For those peaks cover only one gene, run the narrow peak pipeline(please note that the columns may be different from narrow peaks results since broad peak results has no peak summit columns)
  #####genic.peaks1 need promoter peaks info 
  if(nrow(genic.peaks1)==0){
    ## If there are no broad peak covering only 1 gene, no need to put in narrow peaks pipeline
    genic.peaks_final <-genic.peaks2_final 
  }else{
    #dim(genic.peaks2)
    #dim(genic.peaks1)
    #####For those peaks only covered 1 gene, use the same pipeline as narrow peak
    #####For those peaks only covered >=2 genes, stop at here would be enough
    
    ##V23 is length of the intersected broad peaks
    #class(genic.peaks1$covered_genes)
    #####Need to turn list into characters
    genic.peaks1$covered_genes <- as.character(genic.peaks1$covered_genes)
    ######Notice: the peaks may overlap with multiple exons and introns
    exon.peaks <- bedtoolsr::bt.intersect(a = genic.peaks1, b = exon_bed, wo = T)
    if (nrow(exon.peaks) == 0) {
      exon.peaks <- data.frame(Category = character())
    } else {
      exon.peaks$Category <- "lncRNA"
    }
    
    promoter.peaks <- bedtoolsr::bt.intersect(a = genic.peaks1, b = Promoter_bed2, wo = T)
    if (nrow(promoter.peaks) == 0) {
      promoter.peaks<- data.frame(Category = character())
    } else {
      promoter.peaks$Category <- "Promoter"
    }
    
    genic.peaks1_final <- rbind(exon.peaks,promoter.peaks)
    genic.peaks1_final2 <- genic.peaks1_final %>%
      transmute(Chr = V1, str = V2, end = V3, peak_id = V4, signalValue1=V11, pValue1=V12, qValue1=V13, peak_id1=V24,
                signalValue2=V20, pValue2=V21, qValue2=V22, peak_id2=V25,GeneID = V29, overlapped=V41,num=1,Category = Category)
    
    genic.peaks1_final2$intersected_peak_width <-genic.peaks1_final2$end-genic.peaks1_final2$str
    genic.peaks1_final2$overlapped_percentage <- genic.peaks1_final2$overlapped/genic.peaks1_final2$intersected_peak_width
    genic.peaks1_final2 <- genic.peaks1_final2%>%group_by(peak_id)%>% mutate(unique_assignment = Category[which.max(overlapped_percentage)]) %>%ungroup()
    
    genic.peaks_final <-rbind(genic.peaks2_final,genic.peaks1_final2)
    
  }
  
  all.peaks <- genic.peaks_final
  all.peaks<- all.peaks[!duplicated(all.peaks), ]
  ####To remove duplicated rows
  write.xlsx(all.peaks, paste0(out.dir,sample_name,"_allpeaks_assign_lncRNA.xlsx"))
  cat(paste("processing", sample_name),"\n")
}



