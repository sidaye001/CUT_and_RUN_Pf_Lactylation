library(dplyr)
library(DiffBind)
library(tidyverse)
library(openxlsx)
library(GenomicRanges)
####This script is for differential binding sites analysis on HPC, narrowpeaks and broadpeaks have different pipelines.

####This part is for narrowpeaks differential binding sites analysis##########
####Notice: when compare 8mM vs UNT, the consensus peaksets should be the peaksets of 8mM only, since those peaksets has signal in 8mM and 0 signal in UNT will be filtered out by default setting. But those binding sites also shows the induction of lactylation.
####Need to create new consensus peaksets
samples <- read.xlsx('/home/sida.ye001/cut_and_run/Input/DBA/narrowpeaks/sample_sheet_DBA_Kla_rings.xlsx')

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
sampleDba_8mM <- dba(sampleDba, mask = sampleDba$samples$Condition == "8mM")
consensus_peaks <- dba.peakset(sampleDba_8mM, consensus = DBA_CONDITION, minOverlap = 1, bRetrieve = TRUE)


#consensus_peaks[seqnames(consensus_peaks) == "Pf3D7_04_v3"]
sampleCount <- dba.count(sampleDba, peaks = consensus_peaks, summits = F, bRemoveDuplicates=T, bSubControl=T,mapQCth=10,filter=20)

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
# Convert to a data frame
#diffReport_df <- as.data.frame(diffReport)
# Alternatively, export to an Excel file
write.xlsx(diffReport_df, "/home/sida.ye001/cut_and_run/Kla_rings_narrow_edgeR.xlsx")


##############The code below is for broad peak diffferential binding sites analysis###################

samples <- read.xlsx('/home/sida.ye001/cut_and_run/Input/DBA/broadpeaks/sample_sheet_DBA_Kla_schizont.xlsx')
samples <- read.xlsx('/home/sida.ye001/cut_and_run/Input/DBA/broadpeaks/sample_sheet_DBA_Kla_rings.xlsx')
samples <- read.xlsx('/home/sida.ye001/cut_and_run/Input/DBA/broadpeaks/sample_sheet_DBA_Kla_trophs.xlsx')

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
sampleDba <- dba(sampleSheet=samples, scoreCol=5)
# Replace the modified consensus peakset back into the DBA object
#Optional: For broad peaks,bUseSummarizeOverlaps() gives a comprehensive count of reads that overlap with the peak, regardless of the precise location within the peak.by simply counting the number of reads that fall into predefined bins or windows within the peak. This method might not accurately capture the full extent of the peak if the peak spans a broad region.
#with "summits=F", it will not use the summit-centric approach for counting reads. Instead, it will focus on counting reads across the entire broad peak region as specified

####Only use the consensus peaks between 8mM broad peaks to perform differential binding analysis
sampleDba_8mM <- dba(sampleDba, mask = sampleDba$samples$Condition == "8mM")
consensus_peaks <- dba.peakset(sampleDba_8mM, consensus = DBA_CONDITION, minOverlap = 1, bRetrieve = TRUE)
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
sampleAnalyze <- dba.analyze(sampleCon)
# Generate the differential binding report,Set th parameter to 1.0 or higher to include all results, by default, th=0.05(FDR)
diffReport_df <- dba.report(sampleAnalyze, th=1.0)
# Convert to a data frame
#diffReport_df <- as.data.frame(diffReport)
# Alternatively, export to an Excel file
write.xlsx(diffReport_df, "/home/sida.ye001/cut_and_run/Kla_rings.xlsx")
write.xlsx(diffReport_df, "/home/sida.ye001/cut_and_run/Kla_trophs.xlsx")
write.xlsx(diffReport_df, "/home/sida.ye001/cut_and_run/Kla_schizont.xlsx")

write.xlsx(diffReport_df, "/home/sida.ye001/cut_and_run/Kla_schizont_broad.xlsx")
write.xlsx(diffReport_df, "/home/sida.ye001/cut_and_run/Kla_rings_broad.xlsx")
write.xlsx(diffReport_df, "/home/sida.ye001/cut_and_run/Kla_trophs_broad.xlsx")

