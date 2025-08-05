library(dplyr)
library(tidyverse)
library(openxlsx)
library(edgeR)
library(statmod)
library(limma)
library(EnhancedVolcano)

Output.dir <- "./Output/Bulk_RNA/"
# Sample information data frame
sample_info <- data.frame(
  condition = factor(c("UNT", "UNT", "UNT", "8mM", "8mM", "8mM", 
                       "UNT", "UNT", "UNT", "8mM", "8mM", "8mM")),
  stage = factor(c("Rings", "Rings", "Rings", "Rings", "Rings", "Rings", 
                   "Schizonts", "Schizonts", "Schizonts", "Schizonts", "Schizonts", "Schizonts"))
)

design <- model.matrix(~ condition + stage, data = sample_info)

#countmatrix <- read.xlsx("./Input/countmatrix/all_counts.xlsx")
countmatrix <- read.xlsx("./Input/countmatrix/all_counts_including_lncRNA.xlsx")
###To manually change the sample names
#colnames(countmatrix)[2:13] <- c("UNT_Rings_R1","8mM_Schizonts_R1","8mM_Schizonts_R2","8mM_Schizonts_R3",
#                                 "UNT_Rings_R2","UNT_Rings_R3","8mM_Rings_R1","8mM_Rings_R2",
#                                 "8mM_Rings_R3","UNT_Schizonts_R1","UNT_Schizonts_R2","UNT_Schizonts_R3")

colnames(countmatrix)[2:13] <- c("UNT_Rings_R1","8mM_Schizonts_R1","8mM_Schizonts_R2","UNT_Schizonts_R3",
                                 "UNT_Rings_R2","UNT_Rings_R3","8mM_Rings_R1","8mM_Rings_R2",
                                 "8mM_Rings_R3","UNT_Schizonts_R1","UNT_Schizonts_R2","8mM_Schizonts_R3")

#countmatrix2 <- countmatrix[,c(1,2,6,7,8,9,10,11,12,13,3,4,5)]
countmatrix2 <- countmatrix[,c(1,2,6,7,8,9,10,11,12,5,3,4,13)]
colnames(countmatrix2)
write.xlsx(countmatrix2, paste0(Output.dir,"bulkRNA_raw_read_counts_corrected.xlsx"))
write.xlsx(countmatrix2, paste0(Output.dir,"bulkRNA_raw_read_counts_pcgenes_only.xlsx"))
####Perform DEG analysis by edgeR####
####To check the similarity between samples
Expr.c1.c2 <- countmatrix2[,c(2:ncol(countmatrix2))]
## Remove rows with low counts
CPM  <- cpm(Expr.c1.c2)

nor_countmatrix2 <- cbind(countmatrix2[,1],as.data.frame(CPM))
colnames(nor_countmatrix2)[1] <- "GeneID"
#write.xlsx(nor_countmatrix2, paste0(Output.dir,"bulkRNA_normalized_read_counts.xlsx"))
#write.xlsx(nor_countmatrix2, paste0(Output.dir,"bulkRNA_normalized_read_counts_pcgenes_only.xlsx"))
#keep <-  rowSums(CPM > 0) >= 0
keep <-  rowSums(CPM > 2) >= 3
Expr.c1.c2 <- Expr.c1.c2[keep, ]
print(paste('genes kept:', length(which(keep == T))))
gene.id <- rownames(Expr.c1.c2)

dge <- DGEList(counts=Expr.c1.c2, genes = gene.id)
#Perform TMM normalization
dge <- calcNormFactors(dge)
dge$samples

Group  <- factor(c(rep("0", 3), rep("1", 3), rep("2", 3), rep("3", 3)))
#determine which groups are compared to each other based on the design matrix used in the differential expression analysis
design <- model.matrix(~Group)

dge <- DGEList(counts=Expr.c1.c2, group = Group, genes = gene.id)
#Perform TMM normalization
#dge <- calcNormFactors(dge)
#dge$samples
##To check if batch effect exists
#creates a multidimensional scaling (MDS) plot, which is used to visualize similarities or differences between RNA-Seq samples based on their gene expression profiles.
plotMDS(dge) ###Rings and Schizonts are clustered separately, but 8mM and UNT are not separated clearly
dge2 <- dist(t(dge$counts))
plot(hclust(dge2))

dge<- estimateDisp(dge, design, robust=TRUE)
plotBCV(dge)

#logCPM <- cpm(dge, log=TRUE, prior.count=3, normalized.lib.sizes=TRUE)
#plotMDS(logCPM, cex = 0.8)
#The quasi-likelihood F-test is a more robust and accurate method that accounts for the uncertainty in dispersion estimates
dge<- estimateDisp(dge, design, robust=TRUE)
fit <- glmQLFit(dge, design, robust=TRUE)
qlf <- glmQLFTest(fit)
qlf <- topTags(qlf, n = Inf)$table
#tab <- topTags(qlf,n=Inf,adjust.method="BH")$table
#tt <- tt$table[, c('logFC', 'FDR')]

DEG <- function(Sample.c1, Sample.c2){
  ####Sample.c1=WT group
  ####Sample.c2=treated group
  Expr.c1 <- data.frame(Sample.c1[,2:ncol(Sample.c1)])
  Expr.c2 <- data.frame(Sample.c2[,2:ncol(Sample.c2)])
  Expr.c1.c2 <- cbind(Expr.c1, Expr.c2)
  rownames(Expr.c1.c2) <- Sample.c1[,1]
  
  ## Remove rows with low counts
  CPM  <- cpm(Expr.c1.c2)
  #keep <-  rowSums(CPM > 0) >= 0
  keep <-  rowSums(CPM > 2) >= 3
  Expr.c1.c2 <- Expr.c1.c2[keep, ]
  print(paste('genes kept:', length(which(keep == T))))
  gene.id <- rownames(Expr.c1.c2)
  Group  <- factor(c(rep("0", ncol(Expr.c1)), rep("1", ncol(Expr.c2))))
  #determine which groups are compared to each other based on the design matrix used in the differential expression analysis
  design <- model.matrix(~Group)
  
  dge <- DGEList(counts=Expr.c1.c2, group = Group, genes = gene.id)
  #Perform TMM normalization
  dge <- calcNormFactors(dge)
  #dge <- estimateGLMCommonDisp(dge, design)
  #dge <- estimateGLMTrendedDisp(dge, design)
  #dge <- estimateGLMTagwiseDisp(dge, design)
  #fit <- glmFit(dge, design)
  #fit <- glmLRT(fit, coef = 2)
  #The quasi-likelihood F-test is a more robust and accurate method that accounts for the uncertainty in dispersion estimates
  dge<- estimateDisp(dge, design, robust=TRUE)
  fit <- glmQLFit(dge, design, robust=TRUE)
  qlf <- glmQLFTest(fit)
  ####QLF will directly give FDR
  qlf <- topTags(qlf, n = Inf)$table
  #tab <- topTags(qlf,n=Inf,adjust.method="BH")$table
  
  #tt <- tt$table[, c('logFC', 'FDR')]
  #Compute adjusted p-values (e.g., using Benjamini-Hochberg correction)
  
  #tab <- topTags(fit,n=Inf,adjust.method="BH")$table
  #remove missing values
  tab <- na.omit(qlf)
  return(tab)
}

###To compare 8mM vs UNT at Rings and Schizonts stages
##Sample.c1 is UNT, Sample.c2 is 8mM
#Sample.c1 <- countmatrix2[,c(1,2,3,4)]
#Sample.c2 <-countmatrix2[,c(1,5,6,7)]

#Sample.c1 <- countmatrix2[,c(1,8,9,10)]
#Sample.c2 <-countmatrix2[,c(1,11,12,13)]
## Detecting Bias
#Expr.c1 <- data.frame(Sample.c1[,2:ncol(Sample.c1)])
#Expr.c2 <- data.frame(Sample.c2[,2:ncol(Sample.c2)])
#Expr.c1.c2 <- cbind(Expr.c1, Expr.c2)
#rownames(Expr.c1.c2) <- Sample.c1[,1]
#Group  <- factor(c(rep("0", ncol(Expr.c1)), rep("1", ncol(Expr.c2))))

tab_Rings <- DEG(countmatrix2[,c(1,2,3,4)],countmatrix2[,c(1,5,6,7)])
tab_Schizonts <- DEG(countmatrix2[,c(1,8,9,10)],countmatrix2[,c(1,11,12,13)])

hist(tab_Rings$PValue, n=50)
hist(tab_Schizonts$PValue, n=50)
hist(tab_Rings$FDR, n=100)
hist(tab_Schizonts$FDR, n=100)



#write.xlsx(tab_Rings, paste0(Output.dir,"Rings_8mM_vs_UNT.xlsx"))
#write.xlsx(tab_Schizonts, paste0(Output.dir,"Schizonts_8mM_vs_UNT.xlsx"))

#write.xlsx(tab_Rings, paste0(Output.dir,"Rings_8mM_vs_UNT_pcgenes_only.xlsx"))
#write.xlsx(tab_Schizonts, paste0(Output.dir,"Schizonts_8mM_vs_UNT_pcgenes_only.xlsx"))

#write.xlsx(tab_Rings, paste0(Output.dir,"Rings_8mM_vs_UNT_corrected.xlsx"))
#write.xlsx(tab_Schizonts, paste0(Output.dir,"Schizonts_8mM_vs_UNT_corrected.xlsx"))

write.xlsx(tab_Rings, paste0(Output.dir,"Rings_8mM_vs_UNT_pcgenes_only_corrected.xlsx"))
write.xlsx(tab_Schizonts, paste0(Output.dir,"Schizonts_8mM_vs_UNT_pcgenes_only_corrected.xlsx"))

write.xlsx(tab_Rings, paste0(Output.dir,"Rings_8mM_vs_UNT_pcgenes_only_corrected_including_lncRNA.xlsx"))
write.xlsx(tab_Schizonts, paste0(Output.dir,"Schizonts_8mM_vs_UNT_pcgenes_only_corrected_including_lncRNA.xlsx"))

tab_Rings <- read.xlsx(paste0(Output.dir,"Rings_8mM_vs_UNT_pcgenes_only_corrected.xlsx"))
tab_Schizonts <- read.xlsx(paste0(Output.dir,"Schizonts_8mM_vs_UNT_pcgenes_only_corrected.xlsx"))

tab_Rings <- read.xlsx(paste0(Output.dir,"Rings_8mM_vs_UNT_pcgenes_only_corrected_including_lncRNA.xlsx"))
tab_Schizonts <- read.xlsx(paste0(Output.dir,"Schizonts_8mM_vs_UNT_pcgenes_only_corrected_including_lncRNA.xlsx"))

total.product.Pf <- read.csv("./Input/5720_total_Pf_product_description.csv")
total.product.Pf <-total.product.Pf[,c(1,3,4)]
colnames(total.product.Pf)
colnames(total.product.Pf) <- c("genes","Product.Description","Symbol")

tab_Rings2 <- left_join(tab_Rings,total.product.Pf, by="genes")
tab_Schizonts2 <- left_join(tab_Schizonts,total.product.Pf, by="genes")

write.xlsx(tab_Rings2, paste0(Output.dir,"Rings_8mM_vs_UNT_pcgenes_only_corrected2.xlsx"))
write.xlsx(tab_Schizonts2, paste0(Output.dir,"Schizonts_8mM_vs_UNT_pcgenes_only_corrected2.xlsx"))

write.xlsx(tab_Rings2, paste0(Output.dir,"Rings_8mM_vs_UNT_pcgenes_only_corrected_including_lncRNA2.xlsx"))
write.xlsx(tab_Schizonts2, paste0(Output.dir,"Schizonts_8mM_vs_UNT_pcgenes_only_corrected_including_lncRNA2.xlsx"))

tab_Rings3 <- tab_Rings2%>% mutate(labels=ifelse(is.na(Symbol)|Symbol=="N/A",genes,Symbol))

tab_Schizonts3 <- tab_Schizonts2%>% mutate(labels=ifelse(is.na(Symbol)|Symbol=="N/A",genes,Symbol))
EnhancedVolcano(tab_Rings3,
                title = "",
                subtitle = "8mM vs UNT: Rings",
                lab = tab_Rings3$labels,
                pCutoff = 0.05,         # Significance threshold for FDR/PValue
                FCcutoff = 1,           # Threshold for fold change
                x = 'logFC',
                y = 'PValue',
                legendLabels = c('NS', 'Fold Change', 'PValue', 'PValue & Fold Change'),
                drawConnectors = TRUE,    # Draw lines connecting labels to points
                widthConnectors = 0.5,    # Width of the connector lines
                colConnectors = 'black',   # Color of the connector lines
                max.overlaps=60,
                ylab = "-log10(Pvalue)")

EnhancedVolcano(tab_Schizonts3,
                title = "",
                subtitle = "8mM vs UNT: Schizonts",
                lab = tab_Schizonts3$labels,
                pCutoff = 0.05,         # Significance threshold for FDR/PValue
                FCcutoff = 1.5,           # Threshold for fold change
                x = 'logFC',
                y = 'PValue',
                legendLabels = c('NS', 'Fold Change', 'PValue', 'PValue & Fold Change'),
                drawConnectors = TRUE,    # Draw lines connecting labels to points
                widthConnectors = 0.5,    # Width of the connector lines
                colConnectors = 'black',   # Color of the connector lines
                max.overlaps=60,
                ylab = "-log10(Pvalue)")

###early gametocyte markers
genelist <- c("PF3D7_1335000","PF3D7_1473700","PF3D7_0801900","PF3D7_1102500","PF3D7_0936400","PF3D7_1467600","PF3D7_1472200","PF3D7_1473200","PF3D7_0936500","PF3D7_0423700","PF3D7_0935400","PF3D7_0935390")

filter_table <- tab_Rings3%>%filter(genes%in%genelist)
filter_table <- tab_Schizonts3%>%filter(genes%in%genelist)
