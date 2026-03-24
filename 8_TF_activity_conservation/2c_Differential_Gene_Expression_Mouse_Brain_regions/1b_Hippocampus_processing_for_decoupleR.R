setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Species_Comparison/Mouse_Datasets_for_Comparison/BulkRNA/DESeq2')
options(stringsAsFactors = F)

#### Packages
library('DESeq2')         # 
library('sva')            # 
library('limma')          # 

library(ggplot2)          # 
library(scales)           # 
library("bitops")         # 
library(Vennerable)       # 
library(data.table)       #

library(ComplexHeatmap)   #
library(circlize)         #

theme_set(theme_bw())   

# 2025-08-08
# Parse/annotate Stilling hippocampus Brain data for DecoupleR

###############################################################################################
# 1. Read and preprocess data

hippo.cts.1 <- read.csv("../FASTQ/RNA_regions/Stilling_Hippocampus/STAR/2025-07-30_Stilling_Hippocampus_Aging_counts_MM39.txt", skip = 1, sep = "\t", header = T)
hippo.cts   <- hippo.cts.1[,-c(2:6)]

## clean and prep
rownames(hippo.cts) <- hippo.cts$Geneid
colnames(hippo.cts)[-1] <- gsub("_RNAseq_STAR_Aligned.sortedByCoord.out.bam","",colnames(hippo.cts)[-1])

my.keep <- rowSums(hippo.cts[,-1]>0) > 3

# keep expressed genes and remove gene name column
hippo.cts.cl <- hippo.cts[my.keep, -1]
#############################################################################################

##############################################################################################
# 2. Use SVA to clean up batch effects on expression and DEseq2 for DE analysis

# will run SVA to clean up noise

######################################################
#######   DEG analysis   ++++   Stilling     #######
######################################################

my.meta           <- data.frame(matrix(0,6,4))
colnames(my.meta) <- c("SampleID","Age_Group", "Age_month", "Sex")
my.meta$Age_Group <- factor(c(rep("Y",3),rep("O",3)) , levels = c("Y","O"))
my.meta$Age_month <- c(rep(3,3),rep(19,3))
my.meta$Sex       <- rep(NA,6)
my.meta$SampleID  <- colnames(hippo.cts.cl)
rownames(my.meta) <- my.meta$SampleID
my.meta


# get outprefix
my.outprefix <- paste0(Sys.Date(),"_DEseq2_Stilling_Hippocampus_")

###################################
#######       Run SVA      #######

# build design matrix
sva.dataDesign = data.frame( row.names = my.meta$SampleID , 
                             age       = my.meta$Age_month)

# Set null and alternative models (ignore batch)
mod1    = model.matrix(~ age , data = sva.dataDesign)
n.sv.be = num.sv(hippo.cts.cl, mod1, method="be") #1

# apply SVAseq algortihm
my.svseq = svaseq(as.matrix(hippo.cts.cl), mod1, n.sv=n.sv.be, constant = 0.1)

# remove RIN and SV, preserve age and sex
my.clean <- removeBatchEffect(log2(hippo.cts.cl + 0.1), 
                              covariates = cbind(my.svseq$sv),
                              design     = mod1)

# delog and round data for DEseq2 processing
my.filtered.sva <- round(2^my.clean-0.1)


# legend
my.cols  <- rep("",nrow(my.meta))
my.cols[my.meta$Age_Group %in% "Y"] <- "deepskyblue"
my.cols[my.meta$Age_Group %in% "O"] <- "deepskyblue4"

# get matrix using age as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = my.filtered.sva,
                              colData   = my.meta,
                              design    = ~ Age_month)

# run DESeq normalizations and export results
dds.deseq <- DESeq(dds)

# get DESeq2 normalized expression value
vst.cts <- getVarianceStabilizedData(dds.deseq)

# MDS analysis
mds.result <- cmdscale(1-cor(vst.cts,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

pdf(paste0(my.outprefix,"_MDS_plot.pdf"))
plot(x, y,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main= "Stilling Hippocampus (MDS)",
     cex=3, col = my.cols, pch = 16,
     cex.lab = 1.25,
     cex.axis = 1.25, las = 1)
dev.off()

# extract gene significance by DEseq2
res.age <- results(dds.deseq, name = "Age_month") # FC per month

# exclude genes with NA FDR value
res.age <- res.age[!is.na(res.age$padj),]

# store results
deseq.res.list       <- data.frame(res.age)

### get sex dimorphic changes at FDR5
genes.age <- rownames(res.age)[res.age$padj < 0.05]
my.num.age <- length(genes.age) # 511

if (my.num.age > 2) {
  # heatmap drawing - only if there is at least 2 gene
  my.heatmap.out <- paste0(my.outprefix,"_AGING_Heatmap_FDR5_GENES.pdf")
  
  pdf(my.heatmap.out, onefile = F, height = 10, width = 10)
  my.heatmap.title <- paste0("aging significant (FDR<5%), ", my.num.age, " genes")
  pheatmap::pheatmap(vst.cts[genes.age,],
                     cluster_cols = F,
                     cluster_rows = T,
                     colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                     show_rownames = F, scale="row",
                     main = my.heatmap.title,
                     cellwidth = 15,
                     border    = NA,
                     cellheight = 0.15 )
  dev.off()
}

deseq.res.list <- list("Hippocampus" = data.frame(res.age))

# save R object with all DEseq2 results
my.rdata.age <- paste0(Sys.Date(),"_Stilling_Hippocampus_DEseq2_objects.RData")
save(deseq.res.list, file = my.rdata.age)

my.vst.age <- paste0(Sys.Date(),"_Stilling_Hippocampus_VST_data_objects.RData")
save(vst.cts, file = my.vst.age)
##############################################################################################

#######################
sink(file = paste0(Sys.Date(),"_Stilling_Hippocampus_processing_session_Info.txt"))
sessionInfo()
sink()
