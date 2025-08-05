#BiocManager::install("csaw")
#install.packages("statmod")
library(csaw)
library(edgeR)
library(GenomicRanges)
library(openxlsx)

pe.param <- readParam(max.frag=400, pe="both")

#Input the pathways for each sample's bam file
bamfiles_table <- read.xlsx("/home/sida.ye001/cut_and_run/Input/DBA/sample_sheet_DBA_narrow.xlsx",startRow = 1)

#vector to specify the combination of comparison, the order should be 
odd_vec <- seq(1, 24, by = 2)

#output <- "/home/sida.ye001/cut_and_run/Output/CSAW_slidewindow/differential_binding/"
#output_bed <- "/home/sida.ye001/cut_and_run/Output/CSAW_slidewindow/bed/"

output <- "/home/sida.ye001/cut_and_run/Output/CSAW_slidewindow_final/differential_binding/"
output_bed <- "/home/sida.ye001/cut_and_run/Output/CSAW_slidewindow_final/bed/"
for(i in odd_vec){
  #vec <- c(i,i+1,i+24,i+25)
  stage <- bamfiles_table$Tissue[i]
  antibody <- bamfiles_table$Factor[i]
  bam_files <- c(trimws(bamfiles_table$bamReads[i]), trimws(bamfiles_table$bamReads[i+1]), trimws(bamfiles_table$bamReads[i+24]), trimws(bamfiles_table$bamReads[i+25]))
  
  
  
  # Count reads in windows, Splits the genome into fixed-size windows
  data <- windowCounts(bam_files, width=100, param=pe.param)
  
  # Experimental design
  group <- factor(c("High", "High", "Low", "Low"))
  design <- model.matrix(~0 + group) #Compared directly
  #design <- model.matrix(~group) #Compared to 8mM as baseline
  colnames(design) <- levels(group)
  contrast <- makeContrasts(High - Low, levels=design)
  
  # Optional:Filter low-abundance windows, this step is data-dependent
  #Bins genome into 1000 bp chunks to estimate average background.
  
  #Filters out windows with log2 enrichment ≤ 1.
  bam.files_chip <- c(bam_files[1],bam_files[2])
  bin.size <- 1000L #specify the integer
  #windowCounts() to count how many reads fall into bins of size 1000 bp; bin=TRUE tells the function to create fixed-size bins (as opposed to sliding windows).
  #param=pe.param likely defines how the BAM files are read (e.g., paired-end settings).
  binned.ip <- windowCounts(bam.files_chip, bin=TRUE, width=bin.size, param=pe.param)
  data.ip=data[,1:2]
  #Applies global filtering to the window data, removing low-signal/noise regions.
  filter.stat <- filterWindowsGlobal(data.ip, background=binned.ip)
  
  #Optional:setting minimum threshold on increase or decrease
  keep <- filter.stat$filter > log2(1) ###Only care about increase induction for at least 0 fold
  data.filt <- data[keep,]
  #data.filt <- data
  
  #To normalise the data for different library sizes you need to calculate normalisation factors based on large bins:
  binned <- windowCounts(bam_files, bin=TRUE, width=10000, param=pe.param)
  data.filt <- normFactors(binned, se.out=data.filt)
  
  data.filt$norm.factors
  
  #Detection of DB (differentially bound) windows (in our case, the occupancy sites, as we test for differences in ChIP vs. input):
  data.filt.calc <- asDGEList(data.filt)
  data.filt.calc <- estimateDisp(data.filt.calc, design)
  fit <- glmQLFit(data.filt.calc, design, robust=TRUE)
  results <- glmQLFTest(fit, contrast=contrast)
  
  
  #head(results$table)
  
  ####Set up merge gap as 500bp or 1000bp
  merged <- mergeWindows(rowRanges(data.filt), tol=500L)
  table.combined <- combineTests(merged$id, results$table)
  
  
  #is.sig.region <- table.combined$FDR <= 0.1
  #table(table.combined$direction[is.sig.region])
  all.results <- data.frame(as.data.frame(merged$region)[,1:3], table.combined)
  
  sample_name <- paste0(antibody,"_",stage)
  
  #####Optional: Only keep significant results
  #all.results <- all.results[all.results$PValue<0.05,]
  ####Assigned peaks to genes
  
  
  write.xlsx(all.results,paste0(output,sample_name,"_CSAW_slidingwindow_differentialbinding_result.xlsx"))
  cat(paste("processing", sample_name),"\n")
  
  sig <- all.results
  #####Optional: Only keep significant results
  #sig <- sig[sig$PValue<0.05,]
  sig.up=sig[sig$direction=="up",]
  
  starts=sig.up[,2]-1
  
  sig.up[,2]=starts
  
  sig_bed=sig.up[,c(1,2,3)]
  sig_bed$start <- as.numeric(sig_bed$start)
  sig_bed$end <- as.numeric(sig_bed$end)
  #filename="/Users/sidaye/Documents/R/Kla_Schizont_narrow_edgeR_CSAW2_merged.bed"
  options(scipen=999)
  write.table(sig_bed,paste0(output_bed, sample_name,"_CSAW_slidingwindow_significant_differentialbinding.bed" ),sep="\t",col.names=FALSE,quote=FALSE,row.names=FALSE)
}


# Load BAM files
#bam_files <- c("/mathspace/data01/mathbio_lab/parasites/processedData/cut_and_run_lactylation_manish/202311_peakcalling/unfiltered_bam/8mM_treated/8mM_Kla_Rings_BR1_S8_L001.sorted.bam", 
#               "/mathspace/data01/mathbio_lab/parasites/processedData/cut_and_run_lactylation_manish/202311_peakcalling/unfiltered_bam/8mM_treated/8mM_Kla_Rings_BR2_S18_L001.sorted.bam", 
#               "/mathspace/data01/mathbio_lab/parasites/processedData/cut_and_run_lactylation_manish/202311_peakcalling/unfiltered_bam/UNT_treated/UNT_Kla_Rings_BR1_S3_L001.sorted.bam", 
#               "/mathspace/data01/mathbio_lab/parasites/processedData/cut_and_run_lactylation_manish/202311_peakcalling/unfiltered_bam/UNT_treated/UNT_Kla_Rings_BR2_S13_L001.sorted.bam")

# Load BAM files
#bam_files <- c("/mathspace/data01/mathbio_lab/parasites/processedData/cut_and_run_lactylation_manish/202311_peakcalling/unfiltered_bam/8mM_treated/8mM_Kla_Schizont_BR1_S8_L001.sorted.bam", 
#               "/mathspace/data01/mathbio_lab/parasites/processedData/cut_and_run_lactylation_manish/202311_peakcalling/unfiltered_bam/8mM_treated/8mM_Kla_Schizont_BR2_S18_L001.sorted.bam", 
#               "/mathspace/data01/mathbio_lab/parasites/processedData/cut_and_run_lactylation_manish/202311_peakcalling/unfiltered_bam/UNT_treated/UNT_Kla_Schizont_BR1_S3_L001.sorted.bam", 
#               "/mathspace/data01/mathbio_lab/parasites/processedData/cut_and_run_lactylation_manish/202311_peakcalling/unfiltered_bam/UNT_treated/UNT_Kla_Schizont_BR2_S13_L001.sorted.bam")

# Load BAM files
#bam_files <- c("/mathspace/data01/mathbio_lab/parasites/processedData/cut_and_run_lactylation_manish/202311_peakcalling/unfiltered_bam/8mM_treated/8mM_Kla_Rings_BR1_S8_L001.sorted.bam", 
#               "/mathspace/data01/mathbio_lab/parasites/processedData/cut_and_run_lactylation_manish/202311_peakcalling/unfiltered_bam/8mM_treated/8mM_Kla_Rings_BR2_S18_L001.sorted.bam", 
#               "/mathspace/data01/mathbio_lab/parasites/processedData/cut_and_run_lactylation_manish/202311_peakcalling/unfiltered_bam/UNT_treated/UNT_Kla_Rings_BR1_S3_L001.sorted.bam", 
#               "/mathspace/data01/mathbio_lab/parasites/processedData/cut_and_run_lactylation_manish/202311_peakcalling/unfiltered_bam/UNT_treated/UNT_Kla_Rings_BR2_S13_L001.sorted.bam")

