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
setwd("./")
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

broad.peaks.file1<-  "./Input/Occupancy_slidewindow/peak_results/"
#narrow.peaks.file2<-  "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Input/intersect_bed2/UNT_narrowpeak_p005/"
list.files(broad.peaks.file1)
#list.files(narrow.peaks.file2)
count.files1 <- list.files(broad.peaks.file1)
#count.files2 <- list.files(narrow.peaks.file2)

n <- 24
#n <- 6

###############broad peak assignment###################
###############broad peak assignment###################
###############broad peak assignment###################

in.dir <- broad.peaks.file1
#in.dir <- narrow.peaks.file2
#out.dir <- "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Output/gene_peak_assignment_narrowPeak2/unfiltered/8mM/"
out.dir <- "./Output/Occupancy_slidewindow_assign_final/lncRNA/"
#sample_table <- sample_table1
#sample_table <- sample_table2
count.files_input <- count.files1
#count.files_input <- count.files2

#sample_name <- c("H4K12la_rings_Diffbind_assign_lncRNA","H4K12la_schizont_Diffbind_assign_lncRNA","H4K12la_trophs_Diffbind_assign_lncRNA", "Kla_rings_Diffbind_assign_lncRNA","Kla_schizont_Diffbind_assign_lncRNA","Kla_trophs_Diffbind_assign_lncRNA")


exon_bed <- transcript_bed
for (i in 1:n){
  broad.peaks <-try({read.xlsx(paste0(broad.peaks.file1, count.files1[i]))}, silent = TRUE)
  sample_name <- sub("_Occupancy.*", "", count.files1[i])
  #V1,V2,V3 come from rep1, V10,V11,V12 come from rep2, V7/V16 is signalValue, V8/V17 is -log10(pValue), V9/V18 is -log10(qValue)
  #V5/V14 is score range from 0-1000 to show the significance of peaks, broad peaks has no peak summit
  #V20 is the intersection length of two reps, TSS loci is stored in transcript_bed TSS column
  if (inherits(broad.peaks, "try-error")){next}
  #sort the table by intersected regions 
  peaks.all.sort <- broad.peaks %>% arrange(seqnames, start, end)
  if (nrow(peaks.all.sort)==0){next}
  dim(peaks.all.sort)
  peaks.all.sort$peak_id <- paste(peaks.all.sort$seqnames , paste(peaks.all.sort$start, peaks.all.sort$end, sep= "-"), sep = ":")
  colnames(peaks.all.sort)[1] <- "chr"
  
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
  gene.peaks <- genic.peaks%>%group_by(V12)%>% summarise(
    V12=list(V12),
    covered_genes=list(V16),
    num=n()
  )%>%as.data.frame()
  
  gene.peaks$V12 <- unlist(lapply(gene.peaks$V12, unique))
  
  genic.peaks_m <- left_join(genic.peaks,gene.peaks, by="V12")
  
  ######To separate the peaks into two categories: covered 1 genes only and >=2 genes
  genic.peaks2 <- genic.peaks_m%>%dplyr::filter(num>=2)
  genic.peaks1 <- genic.peaks_m%>%dplyr::filter(num==1)
  
  #####No need to process the peaks covered more than one gene
  genic.peaks2_final <- genic.peaks2%>%transmute(Chr = V1, str = V2, end = V3, num.tests = V4, num.up.logFC=V5, num.down.logFC=V6, PValue=V7, FDR=V8,
                                                 direction=V9, rep.test=V10, rep.logFC=V11, binding_ID=V12,GeneID = covered_genes, overlapped=V3-V2,num=num,Category = 'Covered more than one gene', binding_width=V3-V2)
  
  #####To calculate the percentage of overlapped regions over intersected regions#####
  genic.peaks2_final$overlapped_percentage <- genic.peaks2_final$overlapped/genic.peaks2_final$binding_width
  genic.peaks2_final <- genic.peaks2_final%>%group_by(binding_ID)%>% mutate(unique_assignment = Category[which.max(overlapped_percentage)]) %>%ungroup()
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
      transmute(Chr = V1, str = V2, end = V3, num.tests = V4, num.up.logFC=V5, num.down.logFC=V6, PValue=V7, FDR=V8,
                direction=V9, rep.test=V10, rep.logFC=V11, binding_ID=V12, GeneID = V16, overlapped=V19,num=1,Category = Category, binding_width=V3-V2)
    
    genic.peaks1_final2$overlapped_percentage <- genic.peaks1_final2$overlapped/genic.peaks1_final2$binding_width
    genic.peaks1_final2 <- genic.peaks1_final2%>%group_by(binding_ID)%>% mutate(unique_assignment = Category[which.max(overlapped_percentage)]) %>%ungroup()
    
    genic.peaks_final <-rbind(genic.peaks2_final,genic.peaks1_final2)
    
  }
  
  all.peaks <- genic.peaks_final
  all.peaks<- all.peaks[!duplicated(all.peaks), ]
  ####To remove duplicated rows
  write.xlsx(all.peaks, paste0(out.dir,sample_name,"_Occupancy_slidewindow_assign.xlsx"))
  cat(paste("processing", sample_name),"\n")
}



