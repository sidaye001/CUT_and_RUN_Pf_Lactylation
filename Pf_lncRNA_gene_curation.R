library(openxlsx)
require(dplyr)
library(tidyverse)

Pf_lncRNA <- read.xlsx("./Input/Pf_lncRNA/Pf_lncRNA.xlsx")
write.table(Pf_lncRNA, file = "./Input/Pf_lncRNA/Pf_lncRNA.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
colnames(Pf_lncRNA) <- c("V1","V2","V3","V4","V5")

Pf_lncRNA2 <- read.csv("./Input/Pf_lncRNA/12864_2022_9017_MOESM3_ESM.csv")
Pf_lncRNA2 <- Pf_lncRNA2[,c(2,3,4,1,6,5)]
write.table(Pf_lncRNA2, file = "./Input/Pf_lncRNA/Pf_lncRNA2.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

colnames(Pf_lncRNA) <- c("V1","V2","V3","V4","V5")
colnames(Pf_lncRNA2)<- c("V1","V2","V3","V4","V5")
merged_lncRNA <- rbind(Pf_lncRNA,Pf_lncRNA2)
write.table(merged_lncRNA, file = "./Input/Pf_lncRNA/Pf_lncRNA_merged.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

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
write.table(transcript_bed,"./Input/Pf3D7_v63_transcript.bed",quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#sort -k1,1 -k2,2n Pf_lncRNA_merged.bed > Pf_lncRNA_merged_sorted.bed
#sort -k1,1 -k2,2n Pf3D7_v63_transcript.bed > Pf3D7_v63_transcript_sorted.bed
#bedtools closest -a Pf_lncRNA_merged_sorted.bed -b Pf3D7_v63_transcript_sorted.bed -d > lncRNA_closest_genes.bed

