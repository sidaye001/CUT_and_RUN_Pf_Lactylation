library(dplyr)
library(tidyverse)
library(openxlsx)
library(bedtoolsr)
library(rtracklayer)

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

#####Diffbind assignment runs the same gene peak assignment pipeline, please notice that narrow peak assignment is different from broadpeak assignment

setwd("./")
##This script is only for intersect.bed files gene peak assignment
######Broad peak file has only 9 columns while narrow peak has 10 columns(10th column shows the center of the peak)

######every time run bedtoolsr in Rstudio, need to run : options(bedtools.path = "/opt/homebrew/bin/")
options(bedtools.path = "/opt/homebrew/bin/")

#####Notice: please note that there are alternative splicing and multiple UTR for one gene exists#####
#####Notice: please note that first exon include 5UTR#####

####Step1: to get TSS loci and promoter regions for each gene####s
gtf.file <- "./Input/Pf_genome/PlasmoDB-63_Pfalciparum3D7.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
###To extract transcript loci and get TSS loci and promoter regions
transcript_gtf <- gtf%>%dplyr::filter(V3=="transcript")

#####In Pf gff file, the first exon includes 5UTR, so if needs to stratify the categories, it is better to use CDS to show exon regions without 5UTR
#exon_gtf <- gtf%>%dplyr::filter(V3=="exon")
exon_gtf <- gtf%>%dplyr::filter(V3=="CDS")

##transform gtf file into bed file
function_gtf_to_bed <- function(x){
  x <- x%>% dplyr::select(V1, V4, V5, V9, V3, V7)
  x <- x %>% dplyr::rename(V1=V1, V2=V4, V3=V5, V4=V9, V5=V3, V6=V7)
  x$V5 <- '0'
  
  return(x)
}

transcript_bed <- function_gtf_to_bed(transcript_gtf)
transcript_names <- unlist(lapply(strsplit(transcript_bed$V4, ' '), '[[',4))
transcript_names <- unlist(lapply(strsplit(transcript_names, ';'),'[[',1))
transcript_names <- unlist(lapply(strsplit(transcript_names, '\"'),'[[',2))
transcript_bed$V4 <- transcript_names
##To get TSS loci
transcript_bed <- transcript_bed%>%mutate(TSS=ifelse(V6=="+",V2,V3))

exon_bed <- function_gtf_to_bed(exon_gtf)
exon_names <- unlist(lapply(strsplit(exon_bed$V4, ' '), '[[',4))
exon_names <- unlist(lapply(strsplit(exon_names, ';'),'[[',1))
exon_names <- unlist(lapply(strsplit(exon_names, '\"'),'[[',2))
exon_bed$V4 <- exon_names

####To get intron loci from exon_bed
# Sort exons by transcript and start position
exons_sorted <- exon_bed %>%
  arrange(V4, V2)

###Group by chr, geneID and strand
intron_bed <- exons_sorted %>%
  group_by(V1, V4, V6) %>%
  arrange(V2) %>%
  mutate(next_start = lead(V2), next_end = lead(V3)) %>%
  filter(!is.na(next_start)) %>%
  transmute(V1 = V1, V2 = V3, V3 = next_start, V4 = V4, score = ".", V6 = V6) %>%
  filter(V3 > V2)


#######To specify the promoter region######
Promoter_range <- 1000
Promoter_bed <- transcript_bed %>% transmute(V1=V1,TSS=ifelse(V6=="+",V2,V3),Promoter_end=ifelse(V6=="+",TSS-Promoter_range,TSS+Promoter_range),
                                             V4=V4, V5=V5, V6=V6)
###To rearrange the start and end site
Promoter_bed2 <- Promoter_bed%>%transmute(V1=V1,V2=ifelse(V6=="+",Promoter_end,TSS),V3=ifelse(V6=="+",TSS,Promoter_end),
                                          V4=V4, V5=V5, V6=V6)
###Make those negative starting site as 0s
Promoter_bed2$V2 <- ifelse(Promoter_bed2$V2<0,0,Promoter_bed2$V2)

dim(transcript_bed)
dim(Promoter_bed)


gff_file <- "./Input/Pf_genome/PlasmoDB-63_Pfalciparum3D7.gff"
gff <- import(gff_file)

function_gff_to_utr_bed <- function(x){
  # Filter the GFF data to extract UTR information
  utr <- subset(gff, type %in% c("five_prime_UTR", "three_prime_UTR"))
  head(utr)
  
  #####Turn utr region from Granges object into bed file format#####
  seqnames <- as.character(seqnames(utr))
  ranges <- as.data.frame(ranges(utr))
  strand <- as.character(strand(utr))
  type <- mcols(utr)$type
  ID <- unlist(mcols(utr)$Parent)
  
  utr_bed <- data.frame(
    seqnames = seqnames,
    start = ranges$start,
    end = ranges$end,
    ID = ID,
    type = type,
    strand = strand,
    width = ranges$width)
  ###Use \\. to represent the literal dot
  utr_bed$ID <- unlist(lapply(strsplit(utr_bed$ID,"\\."),'[[',1))
  return(utr_bed)
}

utr_bed <- function_gff_to_utr_bed(gff)
utr5_bed1 <- utr_bed%>%dplyr::filter(type=="five_prime_UTR")
utr3_bed1 <- utr_bed%>%dplyr::filter(type=="three_prime_UTR")
utr5_bed <- data.frame(V1=utr5_bed1$seqnames, 
                       V2=utr5_bed1$start,
                       V3=utr5_bed1$end,
                       V4=utr5_bed1$ID,
                       V5=utr5_bed1$width,
                       V6=utr5_bed1$strand)

utr3_bed <- data.frame(V1=utr3_bed1$seqnames, 
                       V2=utr3_bed1$start,
                       V3=utr3_bed1$end,
                       V4=utr3_bed1$ID,
                       V5=utr3_bed1$width,
                       V6=utr3_bed1$strand)

#####To check the distribution of length of utr regions######
#hist(utr_bed$width, nclass=1000)
## using qval = 0.01 and gene size 69350000 bp in macs2 calling 
## using pval = 0.05 and gene size 27000000 bp in macs2 calling, 202308
## narrow peaks called by macs2
#narrow.peaks.file <-  "../cut_and_run_pk/Sample6-MK-Cut-20220903_S6.vs.Sample3-MK-Cut-20220903_S3_peaks.narrowPeak"

#####Same as Diffbind output, in total has 11 columns
broad.peaks.file1<-  "./Input/Occupancy_slidewindow/peak_results/"
#narrow.peaks.file2<-  "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Input/intersect_bed2/UNT_narrowpeak_p005/"
list.files(broad.peaks.file1)
#list.files(narrow.peaks.file2)
count.files1 <- list.files(broad.peaks.file1)
#count.files2 <- list.files(narrow.peaks.file2)

n <- 24

###############broad peak assignment###################
###############broad peak assignment###################
###############broad peak assignment###################

in.dir <- broad.peaks.file1
#in.dir <- narrow.peaks.file2
#out.dir <- "/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/Output/gene_peak_assignment_narrowPeak2/unfiltered/8mM/"
out.dir <- "./Output/Occupancy_slidewindow_assign_final/pcgenes/"
#sample_table <- sample_table1
#sample_table <- sample_table2
count.files_input <- count.files1
#count.files_input <- count.files2

#####Basic logic for assigning the peaks to genes and different locations
#1.There are six categories for assignment: intergenic, exonic, intronic, 5'UTR, 3'UTR and promoter regions(1000bp upstream TSS)
#2.Peaks locating at exonic, intronic, 5'UTR, 3'UTR and promoter will be assigned to a gene, peaks in intergenic will not be assigned to any gene
#3.Each category will be assigned separately, since epigenetic fators do not uniquely function on one gene
#4.exon_bed is for exonic regions, intron_bed is for intronic regions, Promoter_bed is for promoter regions, utr_bed/utr5_bed/utr3_bed is for 5'UTR or 3'UTR regions
#5.TSS loci is in transcript_bed
#6. There are alternative splicing, so multiple UTR exists
#sample_name <- c("H4K12la_rings_Diffbind_assign","H4K12la_schizont_Diffbind_assign","H4K12la_trophs_Diffbind_assign", "Kla_rings_Diffbind_assign","Kla_schizont_Diffbind_assign","Kla_trophs_Diffbind_assign")

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
  
  ## Exons including 5UTR and TSS
  gtf.exon <- gtf %>% dplyr::filter(V3 == 'exon')
  gtf.exon.sort <- gtf.exon %>% arrange(V1, V4, V5)
  parse.str <- strsplit(gtf.exon$V9, split = ' ')
  inds <- unlist(lapply(parse.str , function(x) which(grepl("gene_id", x)) + 1))
  gtf.exon$V9 <- gsub(";", "", unlist(lapply(1:length(inds), function(i) parse.str[[i]][inds[[i]]])))
  gtf.exon$V9 <- gsub("\"", "", gtf.exon$V9)
  gtf.exon <- gtf.exon %>% group_by(V9) %>% mutate(exon.ord = ifelse(V7 == '+', 1:n(), seq(n(), 1, by = -1)),
                                                   multiple.exon = ifelse(n() > 1, T, F))
  
  ## Filter for first exon coordinates (exon1 coordinates)
  tmp.neg <- gtf.exon %>% filter(V7 == "-") %>% group_by(V9) %>%  dplyr::slice(which.max(V5))
  tmp.pos <- gtf.exon %>% filter(V7 == "+") %>% group_by(V9) %>%  dplyr::slice(which.min(V5))
  gtf.exon1 <- bind_rows(tmp.pos, tmp.neg)
  gtf.exon1.sort <- gtf.exon1 %>% arrange(V1, V4, V5)
  
  ####To assign the intergenic peak to genes
  ####Different narrow peaks, the intergenic peaks are non-genic peaks and can be classified as three categories:
  #1. non-genic peaks overlapped with promoter regions
  #2. non-genic peaks has no overlaps with promoter regions and <2000bp to nearest gene and can be assigned to its nearest gene
  #3. The remaining non-genic peaks are intergenic peaks can not be assigned
  #V4 is binding ID(peak ID)
  nongenic.peaks <- peaks.all.sort %>% dplyr::filter(!(peak_id %in% (unique(genic.peaks_m$V12))))
  nongenic.peaks <- nongenic.peaks %>% arrange(chr, start, end) 
  
  nongenic.promoter.peaks <- bedtoolsr::bt.intersect(a = nongenic.peaks, b = Promoter_bed2, wo = T)
  if (nrow(nongenic.promoter.peaks) == 0) {
    nongenic.promoter.peaks<- data.frame(Category = character())
    intergenic.peaks <- nongenic.peaks
  } else {
    nongenic.promoter.peaks$Category <- "Promoter"
    intergenic.peaks <- nongenic.peaks %>% dplyr::filter(!(peak_id %in% (unique(nongenic.promoter.peaks$V12))))
  }
  ###V19 is the overlapped length between consensus peak binding range between promoter boundaries.
  nongenic.promoter.peaks_final <- nongenic.promoter.peaks%>%transmute(Chr = V1, str = V2, end = V3, num.tests = V4, num.up.logFC=V5, num.down.logFC=V6, PValue=V7, FDR=V8,
                                                                       direction=V9, rep.test=V10, rep.logFC=V11, binding_ID=V12, GeneID = V16, overlapped=V19,num=1,Category = Category, binding_width=V3-V2)
  
  #####To calculate the percentage of overlapped regions over intersected regions#####
  nongenic.promoter.peaks_final$overlapped_percentage <- nongenic.promoter.peaks_final$overlapped/nongenic.promoter.peaks_final$binding_width
  nongenic.promoter.peaks_final <- nongenic.promoter.peaks_final%>%group_by(binding_ID)%>% mutate(unique_assignment = Category[which.max(overlapped_percentage)]) %>%ungroup()
  ####If the overlapped_percentage=1 & num=1, it means those peaks are non.genic peaks overlap with promoter regions only
  
  intergenic.peaks <- intergenic.peaks %>% arrange(chr, start, end) 
  dim(intergenic.peaks)
  ##D="b" reports distance with respect to B, k means reports the k closest hits
  intergenic.peaks.genes.dist <- bedtoolsr::bt.closest(a = intergenic.peaks, b = gtf.exon1.sort, D = "b", k = 1)
  ####Need to be modified and colnames are different from narrow peaks since broad peaks have no summit columns
  if (nrow(intergenic.peaks.genes.dist)==0){
    non_genic_final <-nongenic.promoter.peaks_final
  }else{
    gene.id <- unlist(lapply(strsplit(intergenic.peaks.genes.dist$V21, split = ''), function(x){
      n = length(x)
      gene = paste(x[1:(n-1)], collapse = '')
      return(gene)
    }))
    
    intergenic.peaks.genes.dist$V28 <- gene.id 
    ## V23 <= 0 means the peak is at upstream of TSS, V23 >0 means peak is downstream of TES
    
    ##V23 shows the distance to nearest gene
    #####Those intergenic peaks are too far away from the genic regions 2000bp will not be assigned a gene
    intergenic.peaks.genes.dist2<- intergenic.peaks.genes.dist %>% dplyr::mutate(V28=ifelse(abs(V23) < 2000,V28,NA))
    ###If there is no gene assigned, then the num become 0
    intergenic.peaks.genes.dist2$num <- ifelse(is.na(intergenic.peaks.genes.dist2$V28),0,1)
    ###To extract necessary peaks info, the overlapped regions with intergenic should always 1=100%
    intergenic.peaks_final <- intergenic.peaks.genes.dist2 %>% 
      transmute(Chr = V1, str = V2, end = V3, num.tests = V4, num.up.logFC=V5, num.down.logFC=V6, PValue=V7, FDR=V8,
                direction=V9, rep.test=V10, rep.logFC=V11, binding_ID=V12,GeneID = V28, overlapped=V3-V2,num=num,Category = 'Intergenic', binding_width=V3-V2)
    #####Intergenic peaks, the overlapped_percentage=1
    intergenic.peaks_final$overlapped_percentage <- intergenic.peaks_final$overlapped/intergenic.peaks_final$binding_width
    intergenic.peaks_final$unique_assignment <- 'Intergenic'
    non_genic_final <- rbind(nongenic.promoter.peaks_final,intergenic.peaks_final)
  }
  
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
    
    utr5.peaks <- bedtoolsr::bt.intersect(a = genic.peaks1, b = utr5_bed, wo = T)
    # Check if the data frame is empty
    if (nrow(utr5.peaks) == 0) {
      # If the data frame is empty, create it with the new column
      utr5.peaks <- data.frame(Category = character())
    } else {
      # If the data frame is not empty, add the new column
      utr5.peaks$Category <- "UTR5"
    }
    utr3.peaks <- bedtoolsr::bt.intersect(a = genic.peaks1, b = utr3_bed, wo = T)
    if (nrow(utr3.peaks) == 0) {
      utr3.peaks <- data.frame(Category = character())
    } else {
      utr3.peaks$Category <- "UTR3"
    }
    
    ######Notice: the peaks may overlap with multiple exons and introns
    exon.peaks <- bedtoolsr::bt.intersect(a = genic.peaks1, b = exon_bed, wo = T)
    if (nrow(exon.peaks) == 0) {
      exon.peaks <- data.frame(Category = character())
    } else {
      exon.peaks$Category <- "Exon"
    }
    intron.peaks <- bedtoolsr::bt.intersect(a = genic.peaks1, b = intron_bed, wo = T)
    if (nrow(intron.peaks) == 0) {
      intron.peaks <- data.frame(Category = character())
    } else {
      intron.peaks$Category <- "Intron"
    }
    
    promoter.peaks <- bedtoolsr::bt.intersect(a = genic.peaks1, b = Promoter_bed2, wo = T)
    if (nrow(promoter.peaks) == 0) {
      promoter.peaks<- data.frame(Category = character())
    } else {
      promoter.peaks$Category <- "Promoter"
    }
    genic.peaks1_final <- rbind(utr5.peaks,utr3.peaks,exon.peaks,intron.peaks,promoter.peaks)
    genic.peaks1_final2 <- genic.peaks1_final %>%
      transmute(Chr = V1, str = V2, end = V3, num.tests = V4, num.up.logFC=V5, num.down.logFC=V6, PValue=V7, FDR=V8,
                direction=V9, rep.test=V10, rep.logFC=V11, binding_ID=V12, GeneID = V26, overlapped=V29,num=1,Category = Category,binding_width=V3-V2)
    
    genic.peaks1_final2$overlapped_percentage <- genic.peaks1_final2$overlapped/genic.peaks1_final2$binding_width
    genic.peaks1_final2 <- genic.peaks1_final2%>%group_by(binding_ID)%>% mutate(unique_assignment = Category[which.max(overlapped_percentage)]) %>%ungroup()
    
    genic.peaks_final <-rbind(genic.peaks2_final,genic.peaks1_final2)
    
  }
  
  all.peaks <- rbind(genic.peaks_final,non_genic_final)
  all.peaks<- all.peaks[!duplicated(all.peaks), ]
  ####To remove duplicated rows
  write.xlsx(all.peaks, paste0(out.dir,sample_name,"_Occupancy_slidewindow_assign.xlsx"))
  cat(paste("processing", sample_name),"\n")
}


