setwd("./")
require(Rsamtools)
require(dplyr)
library(reshape2)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(Rmisc)
library(tidyr)
library(InPAS)
library(remotes)
library(GenomicRanges)
library(Gviz)
library(rtracklayer)
library(GenomeInfoDb)
require(BRGenomics)
library(BSgenome)
library(Biostrings)
library(openxlsx)

#options(ucscChromosomeNames=FALSE)
#setwd("/Users/sidaye/Documents/R/Cut_Run_Manish/Manish_lac/")
######To forge BSgenome package
#forgeBSgenomeDataPkg("./seqs_srcdir/BSGenome.Pf3D7.PlasmoDB.63-seed.txt")


install.packages("./Input/Pf_genome/BSGenome.Pf3D7.PlasmoDB.51_1.0.tar.gz", repos = NULL, type="source")
library(BSGenome.Pf3D7.PlasmoDB.51)
Pf3D7.seqinfo<-data.frame(cbind(Chr=seqlevels(BSGenome.Pf3D7.PlasmoDB.51),chr=c(paste0("Chr",substring(names(seqinfo(BSGenome.Pf3D7.PlasmoDB.51))[1:14],7,8)),"Mito","Apico"), data.frame(seqinfo(BSGenome.Pf3D7.PlasmoDB.51))))
genome=Pf3D7

renameChrs<-function(gr){  
  gr<-data.frame(gr)
  gr$seqnames<-as.character(gr$seqnames)
  if(grepl("Pf3D7_",gr$seqnames[1],fixed=T)){
    gr$seqnames<-Pf3D7.seqinfo$chr[match(gr$seqnames,Pf3D7.seqinfo$Chr)]
    gr<-makeGRangesFromDataFrame(data.frame(gr),keep.extra.columns = T)
    seqlengths(gr)<-Pf3D7.seqinfo$seqlengths[match(names(seqlengths(gr)),Pf3D7.seqinfo$chr)]
    isCircular(gr)<-Pf3D7.seqinfo$isCircular[match(names(seqlengths(gr)),Pf3D7.seqinfo$chr)]
    genome(gr)<-"51"
  }else{
    gr$seqnames<-Pf3D7.seqinfo$Chr[match(gr$seqnames,Pf3D7.seqinfo$chr)]
    gr<-makeGRangesFromDataFrame(data.frame(gr),keep.extra.columns = T)
    seqlengths(gr)<-Pf3D7.seqinfo$seqlengths[match(names(seqlengths(gr)),Pf3D7.seqinfo$Chr)]
    isCircular(gr)<-Pf3D7.seqinfo$isCircular[match(names(seqlengths(gr)),Pf3D7.seqinfo$Chr)]
    genome(gr)<-"51"
  }
  return(sort(gr))
}

GFF<-sort(import.gff3("./Input/Pf_genome/PlasmoDB-51_Pfalciparum3D7.gff"))
seqlengths(GFF)[1:16]<-seqlengths(Pf3D7)[c(1:14,16,15)]
isCircular(GFF)<-c(rep(F,14),T,T)
genome(GFF)<-genome(Pf3D7)
seqinfo(GFF)

# modify GFF for compatibility with later functions
gff<-renameChrs(GFF)
gff$feature<-gff$type
names(gff)<-gff$ID

gff$gene<-substr(gff$ID,regexpr("PF3D7_",gff$ID),regexpr("PF3D7_",gff$ID)+12)
gff[gff@seqnames %in% c("Apico","Mito")]$gene<-substr(gff[gff@seqnames %in% c("Apico","Mito")]$ID,regexpr("PF3D7_",gff[gff@seqnames %in% c("Apico","Mito")]$ID),regexpr("PF3D7_",gff[gff@seqnames %in% c("Apico","Mito")]$ID)+13)
gff$transcript<-gff$ID
gff[gff$type %in% c("protein_coding_gene","pseudogene","ncRNA_gene")]$transcript<-NA
gff[gff$type %in% c("exon","CDS","five_prime_UTR","three_prime_UTR")]$transcript<-gff[gff$type %in% c("exon","CDS","five_prime_UTR","three_prime_UTR")]$Parent
gff$exon<-NA; gff$exon[gff$feature=="exon"]<-gff$ID[gff$feature=="exon"]
gff$source<-NULL; gff$phase<-NULL; gff$Note<-NULL; gff$type<-NULL; gff$protein_source_id<-NULL; gff$score<-NULL


########batch processing#########
narrow.peaks.file1<-  "./Input/intersect_bed2/8mM_narrowpeak_p005/"
narrow.peaks.file2<-  "./Input/intersect_bed2/UNT_narrowpeak_p005/"
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
##################################

#read in FE bedgraphs
####Iput Kla
track11 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_Kla_Rings_BR1_S8_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track12 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_Kla_Rings_BR2_S18_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track21 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_Kla_Trophs_BR1_S8_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track22 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_Kla_Trophs_BR2_S18_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track31 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_Kla_Schizont_BR1_S8_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track32 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_Kla_Schizont_BR2_S18_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track41 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_Kla_Rings_BR1_S3_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track42 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_Kla_Rings_BR2_S13_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track51 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_Kla_Trophs_BR1_S3_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track52 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_Kla_Trophs_BR2_S13_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track61 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_Kla_Schizont_BR1_S3_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track62 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_Kla_Schizont_BR2_S13_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

####Iput H4K12la
track71 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H4K12la_Rings_BR1_S9_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track72 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H4K12la_Rings_BR2_S19_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track81 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H4K12la_Trophs_BR1_S9_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track82 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H4K12la_Trophs_BR2_S19_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track91 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H4K12la_Schizont_BR1_S9_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track92 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H4K12la_Schizont_BR2_S19_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track101 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H4K12la_Rings_BR1_S4_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track102 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H4K12la_Rings_BR2_S14_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track111 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H4K12la_Trophs_BR1_S4_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track112 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H4K12la_Trophs_BR2_S14_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track121 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H4K12la_Schizont_BR1_S4_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track122 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H4K12la_Schizont_BR2_S14_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

####Iput H3k4me3
track131 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H3K4me3_Rings_BR1_S7_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track132 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H3K4me3_Rings_BR2_S17_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track141 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H3K4me3_Trophs_BR1_S7_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track142 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H3K4me3_Trophs_BR2_S17_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track151 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H3K4me3_Schizont_BR1_S7_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track152 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H3K4me3_Schizont_BR2_S17_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track161 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H3K4me3_Rings_BR1_S2_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track162 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H3K4me3_Rings_BR2_S12_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track171 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H3K4me3_Trophs_BR1_S2_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track172 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H3K4me3_Trophs_BR2_S12_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track181 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H3K4me3_Schizont_BR1_S2_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track182 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H3K4me3_Schizont_BR2_S12_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

####Iput H4k12ac
track191 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H4K12ac_Rings_BR1_S10_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track192 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H4K12ac_Rings_BR2_S20_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track201 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H4K12ac_Trophs_BR1_S10_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track202 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H4K12ac_Trophs_BR2_S20_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track211 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H4K12ac_Schizont_BR1_S10_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track212 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/8mM_narrowpeak_p005/SPMR_FE/8mM_H4K12ac_Schizont_BR2_S20_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track221 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H4K12ac_Rings_BR1_S5_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track222 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H4K12ac_Rings_BR2_S15_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track231 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H4K12ac_Trophs_BR1_S5_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track232 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H4K12ac_Trophs_BR2_S15_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

track241 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H4K12ac_Schizont_BR1_S5_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
track242 <-renameChrs(keepSeqlevels(import.bedGraph("./Input/bdg_FE/UNT_narrowpeak_p005/SPMR_FE/UNT_H4K12ac_Schizont_BR2_S15_L001.SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))


plotme<-data.frame(sort(subsetByOverlaps(gff[gff$feature%in% c("CDS","ncRNA"),],track11)))
plotme$STRAND<-"plus"
plotme$STRAND[plotme$strand=="-"]<-"minus"
class(plotme$seqnames)

#############Be careful in normalization
track11m <- track11
track11m$score <- log2(track11m$score)
track11m$score[track11m$score<0] <- 0
track11m$score <- track11m$score/max(track11m$score)

track12m <- track12
track12m$score <- log2(track12m$score)
track12m$score[track12m$score<0] <- 0
track12m$score <- track12m$score/max(track12m$score)

track21m <- track21
track21m$score <- log2(track21m$score)
track21m$score[track21m$score<0] <- 0
track21m$score <- track21m$score/max(track21m$score)

track22m <- track22
track22m$score <- log2(track22m$score)
track22m$score[track22m$score<0] <- 0
track22m$score <- track22m$score/max(track22m$score)

track31m <- track31
track31m$score <- log2(track31m$score)
track31m$score[track31m$score<0] <- 0
track31m$score <- track31m$score/max(track31m$score)

track32m <- track32
track32m$score <- log2(track32m$score)
track32m$score[track32m$score<0] <- 0
track32m$score <- track32m$score/max(track32m$score)

track41m <- track41
track41m$score <- log2(track41m$score)
track41m$score[track41m$score<0] <- 0
track41m$score <- track41m$score/max(track41m$score)

track42m <- track42
track42m$score <- log2(track42m$score)
track42m$score[track42m$score<0] <- 0
track42m$score <- track42m$score/max(track42m$score)

track51m <- track51
track51m$score <- log2(track51m$score)
track51m$score[track51m$score<0] <- 0
track51m$score <- track51m$score/max(track51m$score)

track52m <- track52
track52m$score <- log2(track52m$score)
track52m$score[track52m$score<0] <- 0
track52m$score <- track52m$score/max(track52m$score)

track61m <- track61
track61m$score <- log2(track61m$score)
track61m$score[track61m$score<0] <- 0
track61m$score <- track61m$score/max(track61m$score)

track62m <- track62
track62m$score <- log2(track62m$score)
track62m$score[track62m$score<0] <- 0
track62m$score <- track62m$score/max(track62m$score)

track71m <- track71
track71m$score <- log2(track71m$score)
track71m$score[track71m$score<0] <- 0
track71m$score <- track71m$score/max(track71m$score)

track72m <- track72
track72m$score <- log2(track72m$score)
track72m$score[track72m$score<0] <- 0
track72m$score <- track72m$score/max(track72m$score)

track81m <- track81
track81m$score <- log2(track81m$score)
track81m$score[track81m$score<0] <- 0
track81m$score <- track81m$score/max(track81m$score)

track82m <- track82
track82m$score <- log2(track82m$score)
track82m$score[track82m$score<0] <- 0
track82m$score <- track82m$score/max(track82m$score)

track91m <- track91
track91m$score <- log2(track91m$score)
track91m$score[track91m$score<0] <- 0
track91m$score <- track91m$score/max(track91m$score)

track92m <- track92
track92m$score <- log2(track92m$score)
track92m$score[track92m$score<0] <- 0
track92m$score <- track92m$score/max(track92m$score)

track101m <- track101
track101m$score <- log2(track101m$score)
track101m$score[track101m$score<0] <- 0
track101m$score <- track101m$score/max(track101m$score)

track102m <- track102
track102m$score <- log2(track102m$score)
track102m$score[track102m$score<0] <- 0
track102m$score <- track102m$score/max(track102m$score)

track111m <- track111
track111m$score <- log2(track111m$score)
track111m$score[track111m$score<0] <- 0
track111m$score <- track111m$score/max(track111m$score)

track112m <- track112
track112m$score <- log2(track112m$score)
track112m$score[track112m$score<0] <- 0
track112m$score <- track112m$score/max(track112m$score)

track121m <- track121
track121m$score <- log2(track121m$score)
track121m$score[track121m$score<0] <- 0
track121m$score <- track121m$score/max(track121m$score)

track122m <- track122
track122m$score <- log2(track122m$score)
track122m$score[track122m$score<0] <- 0
track122m$score <- track122m$score/max(track122m$score)

track131m <- track131
track131m$score <- log2(track131m$score)
track131m$score[track131m$score<0] <- 0
track131m$score <- track131m$score/max(track131m$score)

track132m <- track132
track132m$score <- log2(track132m$score)
track132m$score[track132m$score<0] <- 0
track132m$score <- track132m$score/max(track132m$score)

track141m <- track141
track141m$score <- log2(track141m$score)
track141m$score[track141m$score<0] <- 0
track141m$score <- track141m$score/max(track141m$score)

track142m <- track142
track142m$score <- log2(track142m$score)
track142m$score[track142m$score<0] <- 0
track142m$score <- track142m$score/max(track142m$score)

track151m <- track151
track151m$score <- log2(track151m$score)
track151m$score[track151m$score<0] <- 0
track151m$score <- track151m$score/max(track151m$score)

track152m <- track152
track152m$score <- log2(track152m$score)
track152m$score[track152m$score<0] <- 0
track152m$score <- track152m$score/max(track152m$score)

track161m <- track161
track161m$score <- log2(track161m$score)
track161m$score[track161m$score<0] <- 0
track161m$score <- track161m$score/max(track161m$score)

track162m <- track162
track162m$score <- log2(track162m$score)
track162m$score[track162m$score<0] <- 0
track162m$score <- track162m$score/max(track162m$score)

track171m <- track171
track171m$score <- log2(track171m$score)
track171m$score[track171m$score<0] <- 0
track171m$score <- track171m$score/max(track171m$score)

track172m <- track172
track172m$score <- log2(track172m$score)
track172m$score[track172m$score<0] <- 0
track172m$score <- track172m$score/max(track172m$score)

track181m <- track181
track181m$score <- log2(track181m$score)
track181m$score[track181m$score<0] <- 0
track181m$score <- track181m$score/max(track181m$score)

track182m <- track182
track182m$score <- log2(track182m$score)
track182m$score[track182m$score<0] <- 0
track182m$score <- track182m$score/max(track182m$score)

track191m <- track191
track191m$score <- log2(track191m$score)
track191m$score[track191m$score<0] <- 0
track191m$score <- track191m$score/max(track191m$score)

track192m <- track192
track192m$score <- log2(track192m$score)
track192m$score[track192m$score<0] <- 0
track192m$score <- track192m$score/max(track192m$score)

track201m <- track201
track201m$score <- log2(track201m$score)
track201m$score[track201m$score<0] <- 0
track201m$score <- track201m$score/max(track201m$score)

track202m <- track202
track202m$score <- log2(track202m$score)
track202m$score[track202m$score<0] <- 0
track202m$score <- track202m$score/max(track202m$score)

track211m <- track211
track211m$score <- log2(track211m$score)
track211m$score[track211m$score<0] <- 0
track211m$score <- track211m$score/max(track211m$score)

track212m <- track212
track212m$score <- log2(track212m$score)
track212m$score[track212m$score<0] <- 0
track212m$score <- track212m$score/max(track212m$score)

track221m <- track221
track221m$score <- log2(track221m$score)
track221m$score[track221m$score<0] <- 0
track221m$score <- track221m$score/max(track221m$score)

track222m <- track222
track222m$score <- log2(track222m$score)
track222m$score[track222m$score<0] <- 0
track222m$score <- track222m$score/max(track222m$score)

track231m <- track231
track231m$score <- log2(track231m$score)
track231m$score[track231m$score<0] <- 0
track231m$score <- track231m$score/max(track231m$score)

track232m <- track232
track232m$score <- log2(track232m$score)
track232m$score[track232m$score<0] <- 0
track232m$score <- track232m$score/max(track232m$score)

track241m <- track241
track241m$score <- log2(track241m$score)
track241m$score[track241m$score<0] <- 0
track241m$score <- track241m$score/max(track241m$score)

track242m <- track242
track242m$score <- log2(track242m$score)
track242m$score[track242m$score<0] <- 0
track242m$score <- track242m$score/max(track242m$score)

####
#Pf3D7_04_v3:1,104,184-1,142,000
#Pf3D7_01_v3:488,082-514,518

##########Option1:Choose the specific loci for each gene##########
###gdv1
Chrname <- "Chr09"
###var genes
Chrname <- "Chr04"
###METAP1a
Chrname <- "Chr05"
###METAP1a
Chrname <- "Chr07"


#Pf3D7_04_v3:1,104,184-1,142,000
#Pf3D7_01_v3:488,082-514,518
Chrname <- "Chr04"
Chrname <- "Chr01"

##########Option2:Input the list of regions need to be plotted##########
region_list <- read.xlsx("./Input/Cut&Run_IGV_Schematics.xlsx")
region_list$Chr <- unlist(lapply(strsplit(region_list$binding_ID, ":"),"[[",1))
region_list$Loci <- unlist(lapply(strsplit(region_list$binding_ID, ":"),"[[",2))

region_list$start <- unlist(lapply(strsplit(region_list$Loci, "-"),"[[",1))
region_list$start <- as.numeric(gsub(",", "", region_list$start))

region_list$end <- unlist(lapply(strsplit(region_list$Loci, "-"),"[[",2))
region_list$end <- as.numeric(gsub(",", "", region_list$end))

extension <- 100
region_list$start2 <- region_list$start-extension
region_list$end2 <- region_list$end+extension

region_list$Chrname <- paste0("Chr",unlist(lapply(strsplit(region_list$Chr, "_"),"[[",2)))

for (i in 1:nrow(region_list)){
  Chrname <- region_list$Chrname[i]
  start_loci <- region_list$start2[i]
  end_loci <- region_list$end2[i]
  geneID <- region_list$Genes[i]
  track11a <- DataTrack(
    track11m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_Kla_Rings_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track12a <- DataTrack(
    track12m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_Kla_Rings_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track21a <- DataTrack(
    track21m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_Kla_Trophs_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track22a <- DataTrack(
    track22m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_Kla_Trophs_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track31a <- DataTrack(
    track31m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_Kla_Schizont_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track32a <- DataTrack(
    track32m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_Kla_Schizont_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track41a <- DataTrack(
    track41m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_Kla_Rings_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track42a <- DataTrack(
    track42m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_Kla_Rings_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track51a <- DataTrack(
    track51m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_Kla_Trophs_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track52a <- DataTrack(
    track52m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_Kla_Trophs_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track61a <- DataTrack(
    track61m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_Kla_Schizont_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track62a <- DataTrack(
    track62m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_Kla_Schizont_BR2", ylim=c(-.1,1.1), 
  ) 
  
  
  track71a <- DataTrack(
    track71m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_H4K12la_Rings_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track72a <- DataTrack(
    track72m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_H4K12la_Rings_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track81a <- DataTrack(
    track81m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_H4K12la_Trophs_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track82a <- DataTrack(
    track82m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_H4K12la_Trophs_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track91a <- DataTrack(
    track91m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_H4K12la_Schizont_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track92a <- DataTrack(
    track92m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_H4K12la_Schizont_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track101a <- DataTrack(
    track101m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_H4K12la_Rings_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track102a <- DataTrack(
    track102m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_H4K12la_Rings_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track111a <- DataTrack(
    track111m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_H4K12la_Trophs_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track112a <- DataTrack(
    track112m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_H4K12la_Trophs_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track121a <- DataTrack(
    track121m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_H4K12la_Schizont_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track122a <- DataTrack(
    track122m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_H4K12la_Schizont_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track131a <- DataTrack(
    track131m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 1000,
    name="8mM_H3K4me3_Rings_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track132a <- DataTrack(
    track132m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 1000,
    name="8mM_H3K4me3_Rings_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track141a <- DataTrack(
    track141m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 1000,
    name="8mM_H3K4me3_Trophs_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track142a <- DataTrack(
    track142m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 1000,
    name="8mM_H3K4me3_Trophs_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track151a <- DataTrack(
    track151m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 1000,
    name="8mM_H3K4me3Schizont_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track152a <- DataTrack(
    track152m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 1000,
    name="8mM_H3K4me3_Schizont_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track161a <- DataTrack(
    track161m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 1000,
    name="UNT_H3K4me3_Rings_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track162a <- DataTrack(
    track162m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 1000,
    name="UNT_H3K4me3_Rings_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track171a <- DataTrack(
    track171m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 1000,
    name="UNT_H3K4me3_Trophs_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track172a <- DataTrack(
    track172m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 1000,
    name="UNT_H3K4me3_Trophs_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track181a <- DataTrack(
    track181m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 1000,
    name="UNT_H3K4me3_Schizont_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track182a <- DataTrack(
    track182m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 1000,
    name="UNT_H3K4me3_Schizont_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track191a <- DataTrack(
    track191m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_H4K12ac_Rings_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track192a <- DataTrack(
    track192m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_H4K12ac_Rings_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track201a <- DataTrack(
    track201m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_H4K12ac_Trophs_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track202a <- DataTrack(
    track202m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_H4K12ac_Trophs_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track211a <- DataTrack(
    track211m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_H4K12ac_Schizont_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track212a <- DataTrack(
    track212m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#E53528",window=-1,windowSize = 200,
    name="8mM_H4K12ac_Schizont_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track221a <- DataTrack(
    track221m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_H4K12ac_Rings_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track222a <- DataTrack(
    track222m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_H4K12ac_Rings_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track231a <- DataTrack(
    track231m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_H4K12ac_Trophs_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track232a <- DataTrack(
    track232m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_H4K12ac_Trophs_BR2", ylim=c(-.1,1.1), 
  ) 
  
  track241a <- DataTrack(
    track241m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_H4K12ac_Schizont_BR1", ylim=c(-.1,1.1), 
  ) 
  
  track242a <- DataTrack(
    track242m,ucscChromosomeNames=FALSE,chromosome = Chrname,
    type="histogram",col.histogram="#55B7E6",window=-1,windowSize = 200,
    name="UNT_H4K12ac_Schizont_BR2", ylim=c(-.1,1.1), 
  ) 
  
  #####To show more than one bigwig files in one track
  overlayTrack1 <- OverlayTrack(trackList = list(track11a, track41a), name = "Combined Track")
  overlayTrack2 <- OverlayTrack(trackList = list(track12a, track42a), name = "Combined Track")
  overlayTrack3 <- OverlayTrack(trackList = list(track21a, track51a), name = "Combined Track")
  overlayTrack4 <- OverlayTrack(trackList = list(track22a, track52a), name = "Combined Track")
  overlayTrack5 <- OverlayTrack(trackList = list(track31a, track61a), name = "Combined Track")
  overlayTrack6 <- OverlayTrack(trackList = list(track32a, track62a), name = "Combined Track")
  
  overlayTrack7 <- OverlayTrack(trackList = list(track71a, track101a), name = "Combined Track")
  overlayTrack8 <- OverlayTrack(trackList = list(track72a, track102a), name = "Combined Track")
  overlayTrack9 <- OverlayTrack(trackList = list(track81a, track111a), name = "Combined Track")
  overlayTrack10 <- OverlayTrack(trackList = list(track82a, track112a), name = "Combined Track")
  overlayTrack11 <- OverlayTrack(trackList = list(track91a, track121a), name = "Combined Track")
  overlayTrack12 <- OverlayTrack(trackList = list(track92a, track122a), name = "Combined Track")
  
  overlayTrack13 <- OverlayTrack(trackList = list(track131a, track161a), name = "Combined Track")
  overlayTrack14 <- OverlayTrack(trackList = list(track132a, track162a), name = "Combined Track")
  overlayTrack15 <- OverlayTrack(trackList = list(track141a, track171a), name = "Combined Track")
  overlayTrack16 <- OverlayTrack(trackList = list(track142a, track172a), name = "Combined Track")
  overlayTrack17 <- OverlayTrack(trackList = list(track151a, track181a), name = "Combined Track")
  overlayTrack18 <- OverlayTrack(trackList = list(track152a, track182a), name = "Combined Track")
  
  overlayTrack19 <- OverlayTrack(trackList = list(track191a, track221a), name = "Combined Track")
  overlayTrack20 <- OverlayTrack(trackList = list(track192a, track222a), name = "Combined Track")
  overlayTrack21 <- OverlayTrack(trackList = list(track201a, track231a), name = "Combined Track")
  overlayTrack22 <- OverlayTrack(trackList = list(track202a, track232a), name = "Combined Track")
  overlayTrack23 <- OverlayTrack(trackList = list(track211a, track241a), name = "Combined Track")
  overlayTrack24 <- OverlayTrack(trackList = list(track212a, track242a), name = "Combined Track")
  
  plotme2 <- plotme%>%dplyr::filter(seqnames==Chrname)
  CDS.track<-GeneRegionTrack(
    plotme2,
    chr=Chrname,
    feature=plotme2$STRAND,
    col="black",
    name="genes",
    transcriptAnnotation = "gene",
    plus="#ff9200",
    minus="#008d00",
    cex.group=0.75,
    cex.title=1,
    just.group="right",
    fontsize.group=10,
    fontcolor.group="black"
  )
  
  ######double check CDS.track
  #head(CDS.track@range)
  genomeAxis.track <- GenomeAxisTrack(
    name=Chrname,littleTicks=F,col="black",showTitle=T,rotation.title=0,
    fill.range="transparent",fontcolor="black",labelPos="below",cex.title=0.70,col.border.title="white",font.size=14,cex.title=1
  )
  
  #pdf(paste0("./Output/figs/IGV_schematic/", geneID,".pdf"), width = 8, height = 18)
  #####All tracks
  
  #plotTracks(c(overlayTrack13,overlayTrack14,overlayTrack15,overlayTrack16,overlayTrack17,overlayTrack18,
  #             overlayTrack19,overlayTrack20,overlayTrack21,overlayTrack22,overlayTrack23,overlayTrack24,
  #             overlayTrack1,overlayTrack2,overlayTrack3,overlayTrack4,overlayTrack5,overlayTrack6,
  #             overlayTrack7,overlayTrack8,overlayTrack9,overlayTrack10,overlayTrack11,overlayTrack12,CDS.track,genomeAxis.track),
  #           from=start_loci,to=end_loci,
  #           background.panel = "white", background.title = "white",col.title = "black",col.axis = "black",fontsize.legend=18, sizes=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,
  #                                                                                                                                     1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,2,0.5)
  #)
  #####Only H4K12 and Pan-Kla tracks
  pdf(paste0("./Output/figs/IGV_schematic2/", geneID,".pdf"), width = 6, height = 12)
  plotTracks(c(overlayTrack1,overlayTrack2,overlayTrack3,overlayTrack4,overlayTrack5,overlayTrack6,
               overlayTrack7,overlayTrack8,overlayTrack9,overlayTrack10,overlayTrack11,overlayTrack12,CDS.track,genomeAxis.track),
             from=start_loci,to=end_loci,
             background.panel = "white", background.title = "white",col.title = "black",col.axis = "black",fontsize.legend=18, sizes=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,2,0.5)
  )
  dev.off()
}


#plotTracks(c(overlayTrack1,overlayTrack2,CDS.track,genomeAxis.track),
#           from=1365000,to=1399000,
#           background.panel = "white", background.title = "white",col.title = "black",col.axis = "black",fontsize.legend=18, sizes=c(1.5,1.5,2,0.5)
#)

#plotTracks(c(overlayTrack1,overlayTrack2,overlayTrack3,overlayTrack4,overlayTrack5,overlayTrack6,CDS.track,genomeAxis.track),
#           from=1365000,to=1399000,
#           background.panel = "white", background.title = "white",col.title = "black",col.axis = "black",fontsize.legend=18, sizes=c(1.5,1.5,1.5,1.5,1.5,1.5,2,0.5)
#)

####Please notice that each zoom in, ranges are different
#####GDV1 from=1365000,to=1399000

#####Zoom in
#####Pf3D7_04_v3:1,104,184-1,142,000
#####Pf3D7_01_v3:488,082-514,518
plotTracks(c(overlayTrack1,overlayTrack2,overlayTrack3,overlayTrack4,overlayTrack5,overlayTrack6,
             overlayTrack7,overlayTrack8,overlayTrack9,overlayTrack10,overlayTrack11,overlayTrack12,CDS.track,genomeAxis.track),
           from=488082,to=514518,
           background.panel = "white", background.title = "white",col.title = "black",col.axis = "black",fontsize.legend=18, sizes=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,2,0.5)
)
#12X6
#####Zoom in
plotTracks(c(overlayTrack13,overlayTrack14,overlayTrack15,overlayTrack16,overlayTrack17,overlayTrack18,
             overlayTrack19,overlayTrack20,overlayTrack21,overlayTrack22,overlayTrack23,overlayTrack24,
             overlayTrack1,overlayTrack2,overlayTrack3,overlayTrack4,overlayTrack5,overlayTrack6,
             overlayTrack7,overlayTrack8,overlayTrack9,overlayTrack10,overlayTrack11,overlayTrack12,CDS.track,genomeAxis.track),
           from=1104184,to=1142000,
           background.panel = "white", background.title = "white",col.title = "black",col.axis = "black",fontsize.legend=18, sizes=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,
                                                                                                                                     1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,2,0.5)
)
#18X8

#####Var genes
plotTracks(c(overlayTrack1,overlayTrack2,overlayTrack3,overlayTrack4,overlayTrack5,overlayTrack6,
             overlayTrack7,overlayTrack8,overlayTrack9,overlayTrack10,overlayTrack11,overlayTrack12,CDS.track,genomeAxis.track),
           from=934500,to=990000,
           background.panel = "white", background.title = "white",col.title = "black",col.axis = "black",fontsize.legend=18, sizes=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,2,0.5)
)
#10X8

plotTracks(c(overlayTrack13,overlayTrack14,overlayTrack15,overlayTrack16,overlayTrack17,overlayTrack18,
             overlayTrack19,overlayTrack20,overlayTrack21,overlayTrack22,overlayTrack23,overlayTrack24,
             overlayTrack1,overlayTrack2,overlayTrack3,overlayTrack4,overlayTrack5,overlayTrack6,
             overlayTrack7,overlayTrack8,overlayTrack9,overlayTrack10,overlayTrack11,overlayTrack12,CDS.track,genomeAxis.track),
           from=934500,to=990000,
           background.panel = "white", background.title = "white",col.title = "black",col.axis = "black",fontsize.legend=18, sizes=c(1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,
                                                                                                                                     1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,2.5,0.5)
)

#16X8

#####METAP1a
plotTracks(c(overlayTrack1,overlayTrack2,overlayTrack3,overlayTrack4,overlayTrack5,overlayTrack6,
             overlayTrack7,overlayTrack8,overlayTrack9,overlayTrack10,overlayTrack11,overlayTrack12,CDS.track,genomeAxis.track),
           from=1135000,to=1142000,
           background.panel = "white", background.title = "white",col.title = "black",col.axis = "black",fontsize.legend=18, sizes=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,2,0.5)
)
#12X6

plotTracks(c(overlayTrack13,overlayTrack14,overlayTrack15,overlayTrack16,overlayTrack17,overlayTrack18,
             overlayTrack19,overlayTrack20,overlayTrack21,overlayTrack22,overlayTrack23,overlayTrack24,
             overlayTrack1,overlayTrack2,overlayTrack3,overlayTrack4,overlayTrack5,overlayTrack6,
             overlayTrack7,overlayTrack8,overlayTrack9,overlayTrack10,overlayTrack11,overlayTrack12,CDS.track,genomeAxis.track),
           from=1365000,to=1399000,
           background.panel = "white", background.title = "white",col.title = "black",col.axis = "black",fontsize.legend=18, sizes=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,
                                                                                                                                     1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,2,0.5)
)
#18X8

#####METAP1a/0927100/0705500
plotTracks(c(overlayTrack1,overlayTrack2,overlayTrack3,overlayTrack4,overlayTrack5,overlayTrack6,
             overlayTrack7,overlayTrack8,overlayTrack9,overlayTrack10,overlayTrack11,overlayTrack12,CDS.track,genomeAxis.track),
           from=276000,to=290000,
           background.panel = "white", background.title = "white",col.title = "black",col.axis = "black",fontsize.legend=18, sizes=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,2,0.5)
)
#12X6

plotTracks(c(overlayTrack13,overlayTrack14,overlayTrack15,overlayTrack16,overlayTrack17,overlayTrack18,
             overlayTrack19,overlayTrack20,overlayTrack21,overlayTrack22,overlayTrack23,overlayTrack24,
             overlayTrack1,overlayTrack2,overlayTrack3,overlayTrack4,overlayTrack5,overlayTrack6,
             overlayTrack7,overlayTrack8,overlayTrack9,overlayTrack10,overlayTrack11,overlayTrack12,CDS.track,genomeAxis.track),
           from=1365000,to=1399000,
           background.panel = "white", background.title = "white",col.title = "black",col.axis = "black",fontsize.legend=18, sizes=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,
                                                                                                                                     1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,2,0.5)
)
#18X8


