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
# Parse/annotate Benayoun 2019 Brain data for DecoupleR
# just brain regions

###############################################################################################
# 1. Read and preprocess data

my.cts.files <- c("../FASTQ/RNA_regions/Benayoun/STAR/2025-08-02_Benayoun_Cerebellum_Aging_counts_MM39.txt",
                  "../FASTQ/RNA_regions/Benayoun/STAR/2025-08-02_Benayoun_OB_Aging_counts_MM39.txt")

my.cl.counts <- vector(mode = "list", length = 2)
names(my.cl.counts) <- c("Cerebellum", "Olfactory_Bulb")

for (i in 1:length(my.cl.counts)) {
  
  my.cts.1 <- read.csv(my.cts.files[[i]], skip = 1, sep = "\t", header = T)
  my.cts   <- my.cts.1[,-c(2:6)]
  
  ## clean and prep
  rownames(my.cts) <- my.cts$Geneid
  colnames(my.cts)[-1] <- unlist(lapply(strsplit(colnames(my.cts)[-1], "_RNA"),'[[',1))
  
  my.keep <- rowSums(my.cts[,-1]>0) > (ncol(my.cts) -1)/2 
  
  # keep expressed genes and remove gene name column
  my.cl.counts[[i]] <- my.cts[my.keep, -1]
}

# lapply(my.cl.counts,nrow)
# $Cerebellum
# [1] 19182
# 
# $Olfactory_Bulb
# [1] 19144
#############################################################################################

##############################################################################################
# 2. Use SVA to clean up batch effects on expression and DEseq2 for DE analysis

# will run SVA to clean up noise

######################################################
#######   DEG analysis   ++++   Benayoun     #######
######################################################

# Create list object to receive SVA normalized counts
sva.cts        <- vector(mode = "list", length = length(my.cl.counts))
names(sva.cts) <- names(my.cl.counts)

# Create list object to receive VST normalized counts
vst.cts        <- vector(mode = "list", length = length(my.cl.counts))
names(vst.cts) <- names(my.cl.counts)

# Create list object to receive DESeq2 results
deseq.res.list        <- vector(mode = "list", length = length(my.cl.counts))
names(deseq.res.list) <- names(my.cl.counts)


for ( i in 1:length(my.cl.counts)) {
  
  # get meta data
  my.meta           <- data.frame(matrix(0,ncol(my.cl.counts[[i]]),4))
  colnames(my.meta) <- c("SampleID","Age_Group", "Age_month", "Sex")
  my.meta$Sex       <- "M"
  my.meta$Age_month <- as.numeric(unlist(lapply(strsplit(colnames(my.cl.counts[[i]]),"_|m",perl = T), '[[',2)))
  my.meta$Age_Group <- ifelse(my.meta$Age_month == 3, "Y", ifelse(my.meta$Age_month == 12, "M", "O"))
  my.meta$Age_Group <- factor(my.meta$Age_Group , levels = c("Y","M","O"))
  my.meta$SampleID  <- colnames(my.cl.counts[[i]])
  rownames(my.meta) <- my.meta$SampleID
  my.meta
  
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_Benayoun_", names(my.cl.counts)[i])
  
  ###################################
  #######       Run SVA      #######
  
  # build design matrix
  sva.dataDesign = data.frame( row.names = my.meta$SampleID , 
                               age       = my.meta$Age_month)
  
  # Set null and alternative models (ignore batch)
  mod1    = model.matrix(~ age , data = sva.dataDesign)
  n.sv.be = num.sv(my.cl.counts[[i]], mod1, method="be") 
  
  if (n.sv.be >0 ) {
    # apply SVAseq algortihm
    my.svseq = svaseq(as.matrix(my.cl.counts[[i]]), mod1, n.sv=n.sv.be, constant = 0.1)
    
    # remove RIN and SV, preserve age and sex
    my.clean <- removeBatchEffect(log2(my.cl.counts[[i]] + 0.1), 
                                  covariates = cbind(my.svseq$sv),
                                  design     = mod1)
    
    # delog and round data for DEseq2 processing
    my.filtered.sva <- round(2^my.clean-0.1)
    
    sva.cts[[i]] <- my.filtered.sva
    
  } else {
    
    sva.cts[[i]] <- my.cl.counts[[i]]
    
  }
  
  # legend
  my.cols  <-  ifelse(my.meta$Age_month == 3, "deepskyblue", ifelse(my.meta$Age_month == 12, "deepskyblue3", "deepskyblue4"))

  # get matrix using age as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = sva.cts[[i]],
                                colData   = my.meta,
                                design    = ~ Age_month)
  
  # run DESeq normalizations and export results
  dds.deseq <- DESeq(dds)
  
  # get DESeq2 normalized expression value
  vst.cts[[i]] <- getVarianceStabilizedData(dds.deseq)
  
  # MDS analysis
  mds.result <- cmdscale(1-cor(vst.cts[[i]],method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  x <- mds.result[, 1]
  y <- mds.result[, 2]
  
  pdf(paste0(my.outprefix,"_MDS_plot.pdf"))
  plot(x, y,
       xlab = "MDS dimension 1", ylab = "MDS dimension 2",
       main= paste0(names(sva.cts)[[i]]," (MDS)"),
       cex=3, col = my.cols, pch = 16,
       cex.lab = 1.25,
       cex.axis = 1.25, las = 1)
  dev.off()
  
  # extract gene significance by DEseq2
  res.age <- results(dds.deseq, name = "Age_month") # FC per month
  
  # exclude genes with NA FDR value
  res.age <- res.age[!is.na(res.age$padj),]
  
  # store results
  deseq.res.list[[i]]       <- data.frame(res.age)
  
  ### get sex dimorphic changes at FDR5
  genes.age <- rownames(res.age)[res.age$padj < 0.05]
  my.num.age <- length(genes.age)
  
  if (my.num.age > 2) {
    # heatmap drawing - only if there is at least 2 gene
    my.heatmap.out <- paste0(my.outprefix,"_AGING_Heatmap_FDR5_GENES.pdf")
    
    pdf(my.heatmap.out, onefile = F, height = 10, width = 10)
    my.heatmap.title <- paste0(names(sva.cts)[[i]], " aging significant (FDR<5%), ", my.num.age, " genes")
    pheatmap::pheatmap(vst.cts[[i]][genes.age,],
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
}

# save R object with all DEseq2 results
# save R object with all DEseq2 results
my.rdata.age <- paste0(Sys.Date(),"_Benayoun_Brain_DEseq2_objects.RData")
save(deseq.res.list, file = my.rdata.age)

my.vst.age <- paste0(Sys.Date(),"_Benayoun_Brain_VST_data_objects.RData")
save(vst.cts, file = my.vst.age)
##############################################################################################

#######################
sink(file = paste0(Sys.Date(),"_Benayoun_Brain_processing_session_Info.txt"))
sessionInfo()
sink()
