library(dplyr)
library(DiffBind)
library(tidyverse)
library(openxlsx)
library(GenomicRanges)
####This script is for differential binding sites analysis on HPC, narrowpeaks and broadpeaks have different pipelines.

####This part is for narrowpeaks differential binding sites analysis##########
####Notice: when compare 8mM vs UNT, the consensus peaksets should be the peaksets of 8mM only, since those peaksets has signal in 8mM and 0 signal in UNT will be filtered out by default setting. But those binding sites also shows the induction of lactylation.
####Need to create new consensus peaksets
samples1 <- read.xlsx('/home/sida.ye001/cut_and_run/Input/DBA/narrowpeaks/sample_sheet_DBA_Kla_rings.xlsx')
samples2 <- read.xlsx('/home/sida.ye001/cut_and_run/Input/DBA/narrowpeaks/sample_sheet_DBA_Kla_trophs.xlsx')
samples3 <- read.xlsx('/home/sida.ye001/cut_and_run/Input/DBA/narrowpeaks/sample_sheet_DBA_Kla_schizont.xlsx')

samples4 <- read.xlsx('/home/sida.ye001/cut_and_run/Input/DBA/narrowpeaks/sample_sheet_DBA_H4K12la_rings.xlsx')
samples5 <- read.xlsx('/home/sida.ye001/cut_and_run/Input/DBA/narrowpeaks/sample_sheet_DBA_H4K12la_trophs.xlsx')
samples6 <- read.xlsx('/home/sida.ye001/cut_and_run/Input/DBA/narrowpeaks/sample_sheet_DBA_H4K12la_schizont.xlsx')

intersected_narrow1 <- "/home/sida.ye001/cut_and_run/Input/Occupancy_merged_pcgene_lncRNA_assignment/narrow/8mM/"
intersected_narrow2 <- "/home/sida.ye001/cut_and_run/Input/Occupancy_merged_pcgene_lncRNA_assignment/narrow/UNT/"
intersected_broad1 <- "/home/sida.ye001/cut_and_run/Input/Occupancy_merged_pcgene_lncRNA_assignment/broad/8mM/"
intersected_broad2 <- "/home/sida.ye001/cut_and_run/Input/Occupancy_merged_pcgene_lncRNA_assignment/broad/UNT/"

count.files <- list.files(intersected_narrow1)

predefined_peaks <- try({read.xlsx(paste0(intersected_narrow1, count.files[10]))}, silent = TRUE)
predefined_gr <- GRanges(
  seqnames = predefined_peaks$Chr,
  ranges = IRanges(predefined_peaks$str, predefined_peaks$end)
)
#predefined_ranges <- predefined_intersect[, 2:4]
# Ensure numeric values for start and end
#predefined_ranges$str <- as.numeric(predefined_ranges$str)
#predefined_ranges$end <- as.numeric(predefined_ranges$end)

# Create GRanges object
#predefined_gr <- GRanges(
#  seqnames = predefined_ranges$Chr,
#  ranges = IRanges(start = predefined_ranges$str, end = predefined_ranges$end)
#)

Diffbind_narrow <- function(samples, predefined){
  ####To check if there is peak file does not exist, please ensure that the paths do not have trailing spaces, as these can cause issues when trying to access the files
  samples <- as.data.frame(lapply(samples, trimws))
  
  for (i in 1:nrow(samples)) {
    peak_file <- samples$Peaks[i]  # Adjust this to your column name for peak files
    if (!file.exists(peak_file)) {
      stop(paste("Peak file does not exist:", peak_file))
    }
    if (file.info(peak_file)$size == 0) {
      stop(paste("Peak file is empty:", peak_file))
    }
  }
  
  
  ##For narrow peaks
  sampleDba <- dba(sampleSheet=samples)
  #Count the reads in the peaks, by default summits=200, to check the windows +200/-200 surrounding summits
  ####To get consensus peaksets of 8mM only
  # Restrict to 8 mM condition
  # Create consensus peakset for 8 mM
  ####Option1:Only use the consensus peaks between 8mM BR1 and BR2 intersected peaks
  #sampleDba_8mM <- dba(sampleDba, mask = sampleDba$samples$Condition == "8mM")
  ##minOverlap = 1 means unions, and minOverlap = 2 means intersections
  #consensus_peaks <- dba.peakset(sampleDba_8mM, consensus = DBA_CONDITION, minOverlap = 2, bRetrieve = TRUE)
  
  
  #consensus_peaks[seqnames(consensus_peaks) == "Pf3D7_04_v3"]
  sampleCount <- dba.count(sampleDba, peaks = predefined, summits = FALSE,bSubControl=T,minCount=0,filter = 0)
  
  
  #sampleCount <- dba.count(sampleDba, peaks = consensus_peaks, summits = F, bRemoveDuplicates=T, bSubControl=T,mapQCth=10,filter=20)
  
  #sampleCount <- dba.count(sampleDba)
  sampleNor <- dba.normalize(sampleCount)
  #sampleNor <- dba.normalize(sampleCount, method=DBA_EDGER)
  ##Perform contrast on normalized data
  #sampleNor <- sampleCount
  sampleCon <- dba.contrast(sampleNor, design = ~ Condition, contrast = c("Condition", "8mM", "UNT"))
  
  # Perform the differential analysis
  #sampleAnalyze <- dba.analyze(sampleCon)
  sampleAnalyze <- dba.analyze(sampleCon, method=DBA_EDGER)
  # Generate the differential binding report,Set th parameter to 1.0 or higher to include all results, by default, th=0.05(FDR)
  #diffReport_df <- dba.report(sampleAnalyze, th=1.0) ### this is for DeSeq2 method
  diffReport_df <- dba.report(sampleAnalyze,method=DBA_EDGER, th=1.0)#### this is for edgeR method
  return(diffReport_df)
}

# Convert to a data frame
#diffReport_df <- as.data.frame(diffReport)
# Alternatively, export to an Excel file
diffReport_df1 <- Diffbind_narrow(samples1, predefined=predefined_gr)
diffReport_df2 <- Diffbind_narrow(samples2)
diffReport_df3 <- Diffbind_narrow(samples3)
diffReport_df4 <- Diffbind_narrow(samples4)
diffReport_df5 <- Diffbind_narrow(samples5)
diffReport_df6 <- Diffbind_narrow(samples6)

write.xlsx(diffReport_df1, "/home/sida.ye001/cut_and_run/Output/DBA/Kla_rings_narrow_edgeR.xlsx")
write.xlsx(diffReport_df2, "/home/sida.ye001/cut_and_run/Output/DBA/Kla_trophs_narrow_edgeR.xlsx")
write.xlsx(diffReport_df3, "/home/sida.ye001/cut_and_run/Output/DBA/Kla_schizont_narrow_edgeR.xlsx")
write.xlsx(diffReport_df4, "/home/sida.ye001/cut_and_run/Output/DBA/H4K12la_rings_narrow_edgeR.xlsx")
write.xlsx(diffReport_df5, "/home/sida.ye001/cut_and_run/Output/DBA/H4K12la_trophs_narrow_edgeR.xlsx")
write.xlsx(diffReport_df6, "/home/sida.ye001/cut_and_run/Output/DBA/H4K12la_schizont_narrow_edgeR.xlsx")


##############The code below is for broad peak diffferential binding sites analysis###################
##############The code below is for broad peak diffferential binding sites analysis###################
##############The code below is for broad peak diffferential binding sites analysis###################
samples11 <- read.xlsx('/home/sida.ye001/cut_and_run/Input/DBA/broadpeaks/sample_sheet_DBA_Kla_rings.xlsx')
samples22 <- read.xlsx('/home/sida.ye001/cut_and_run/Input/DBA/broadpeaks/sample_sheet_DBA_Kla_trophs.xlsx')
samples33 <- read.xlsx('/home/sida.ye001/cut_and_run/Input/DBA/broadpeaks/sample_sheet_DBA_Kla_schizont.xlsx')

samples44 <- read.xlsx('/home/sida.ye001/cut_and_run/Input/DBA/broadpeaks/sample_sheet_DBA_H4K12la_rings.xlsx')
samples55 <- read.xlsx('/home/sida.ye001/cut_and_run/Input/DBA/broadpeaks/sample_sheet_DBA_H4K12la_trophs.xlsx')
samples66 <- read.xlsx('/home/sida.ye001/cut_and_run/Input/DBA/broadpeaks/sample_sheet_DBA_H4K12la_schizont.xlsx')


Diffbind_broad <- function(samples){
  ####To check if there is peak file does not exist, please ensure that the paths do not have trailing spaces, as these can cause issues when trying to access the files
  samples <- as.data.frame(lapply(samples, trimws))
  
  for (i in 1:nrow(samples)) {
    peak_file <- samples$Peaks[i]  # Adjust this to your column name for peak files
    if (!file.exists(peak_file)) {
      stop(paste("Peak file does not exist:", peak_file))
    }
    if (file.info(peak_file)$size == 0) {
      stop(paste("Peak file is empty:", peak_file))
    }
  }
  
  print("All peak files exist and are not empty.")
  
  
  ##"broadPeak" is not specifically supported as a peakFormat (see help page for dba.peakset()).In the broadPeak format, the score is in the fifth column. So you could try:
  ##And because almost all broadpeak scores are 0, so we choose to use 7th column(signal value) as Score in Diffbind analysis
  sampleDba <- dba(sampleSheet=samples, scoreCol=7)
  # Replace the modified consensus peakset back into the DBA object
  #Optional: For broad peaks,bUseSummarizeOverlaps() gives a comprehensive count of reads that overlap with the peak, regardless of the precise location within the peak.by simply counting the number of reads that fall into predefined bins or windows within the peak. This method might not accurately capture the full extent of the peak if the peak spans a broad region.
  #with "summits=F", it will not use the summit-centric approach for counting reads. Instead, it will focus on counting reads across the entire broad peak region as specified
  
  ####Option1:Only use the consensus peaks between 8mM BR1 and BR2 intersected peaks
  sampleDba_8mM <- dba(sampleDba, mask = sampleDba$samples$Condition == "8mM")
  consensus_peaks <- dba.peakset(sampleDba_8mM, consensus = DBA_CONDITION, minOverlap = 2, bRetrieve = TRUE)
  class(consensus_peaks) ###needs to be Granges object
  widths <- width(consensus_peaks)
  filtered_consensus_peaks <- consensus_peaks[widths >= 1000]
  
  # specifies that the peaks in the high concentration consensus peakset are used for counting,minOverlap = 1: This ensures that peaks from the consensus peakset are used as the checking window, but not necessarily all four samples need to overlap.
  #consensus_peaks need to be Granges object
  #sampleCount <- dba.count(sampleDba, peaks = consensus_peaks, bUseSummarizeOverlaps=TRUE, summits=F, bRemoveDuplicates=T, bSubControl=T,mapQCth=10,filter=50)
  sampleCount <- dba.count(sampleDba, peaks = filtered_consensus_peaks, bUseSummarizeOverlaps=TRUE, summits=F, bRemoveDuplicates=T, bSubControl=T,mapQCth=10,filter=20)
  
  
  
  
  #Normaliza the count
  sampleNor <- dba.normalize(sampleCount)
  #sampleNor <- dba.normalize(sampleCount, method=DBA=EDGER)
  ##Perform contrast on normalized data
  sampleCon <- dba.contrast(sampleNor, design = ~ Condition, contrast = c("Condition", "8mM", "UNT"))
  #sampleCon <- dba.contrast(sampleCount, design = ~ Condition, contrast = c("Condition", "8mM", "UNT"),method=DBA=EDGER)
  # Perform the differential analysis
  sampleAnalyze <- dba.analyze(sampleCon, method=DBA_EDGER)
  # Generate the differential binding report,Set th parameter to 1.0 or higher to include all results, by default, th=0.05(FDR)
  diffReport_df <- dba.report(sampleAnalyze,method=DBA_EDGER, th=1.0)
  return(diffReport_df)
}


# Convert to a data frame
#diffReport_df <- as.data.frame(diffReport)
# Alternatively, export to an Excel file

diffReport_df11 <- Diffbind_broad(samples11)
diffReport_df22 <- Diffbind_broad(samples22)
diffReport_df33 <- Diffbind_broad(samples33)
diffReport_df44 <- Diffbind_broad(samples44)
diffReport_df55 <- Diffbind_broad(samples55)
diffReport_df66 <- Diffbind_broad(samples66)

write.xlsx(diffReport_df11, "/home/sida.ye001/cut_and_run/Output/DBA/Kla_rings_broad_edgeR.xlsx")
write.xlsx(diffReport_df22, "/home/sida.ye001/cut_and_run/Output/DBA/Kla_trophs_broad_edgeR.xlsx")
write.xlsx(diffReport_df33, "/home/sida.ye001/cut_and_run/Output/DBA/Kla_schizont_broad_edgeR.xlsx")
write.xlsx(diffReport_df44, "/home/sida.ye001/cut_and_run/Output/DBA/H4K12la_rings_broad_edgeR.xlsx")
write.xlsx(diffReport_df55, "/home/sida.ye001/cut_and_run/Output/DBA/H4K12la_trophs_broad_edgeR.xlsx")
write.xlsx(diffReport_df66, "/home/sida.ye001/cut_and_run/Output/DBA/H4K12la_schizont_broad_edgeR.xlsx")



