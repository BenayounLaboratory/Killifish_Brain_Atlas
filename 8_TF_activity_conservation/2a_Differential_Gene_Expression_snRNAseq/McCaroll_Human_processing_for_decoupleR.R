setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Species_Comparison/Human_Datasets_for_comparison/snRNA/Muscat_DESeq2/McCaroll_dlPFC')
options(stringsAsFactors = F)

#### Packages
library('Seurat')         # 
library(sctransform)      # 
library("singleCellTK")   # 
library("anndata")

library('muscat')         # 
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

# 2025-07-25
# Analyze McCaroll dlPFC data for DecoupleR
# Use pseudobluked data from Steven

################################################################################################################################################################
#### 1. read data and metadata, create clean pseudobulk list

my.raw.folder <- "/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Species_Comparison/Human_Datasets_for_comparison/snRNA/Raw_Data/McCaroll_dlPFC/"

## read meta data table
my.meta.data <- read.table(paste0(my.raw.folder,"SZvillage_donorMetadata.txt"), header = T, sep = "\t") # 191

# remove schizophrenic patients
my.meta.data.flt <- my.meta.data[my.meta.data$Schizophrenia %in% "Unaffected",] # 97

hist(my.meta.data.flt$Age) # 33 to 93 year old

table(my.meta.data.flt$Sex) # 33 females, 64 males

# pseudobulk count files
my.cts.files <- c(paste0(my.raw.folder, "BA46.astrocyte.Protoplasmic.noOutliers.counts.metacells.txt") ,
                  paste0(my.raw.folder, "BA46.gabaergic.All.noOutliers.counts.metacells.txt"         ) ,
                  paste0(my.raw.folder, "BA46.glutamatergic.All.noOutliers.counts.metacells.txt"     ) ,
                  paste0(my.raw.folder, "BA46.microglia.All.noOutliers.counts.metacells.txt"         ) ,
                  paste0(my.raw.folder, "BA46.oligodendrocyte.All.noOutliers.counts.metacells.txt"   ) ,
                  paste0(my.raw.folder, "BA46.polydendrocyte.All.noOutliers.counts.metacells.txt"    ) )

# get list object
cts.list        <- vector(mode = "list", length = length(my.cts.files))
names(cts.list) <- c("Astrocytes", "GABAergic_Neurons", "Glutamatergic_Neurons", "Microglia", "Oligodendrocytes", "Polydendrocytes")
  
for (i in 1:length(my.cts.files)) {
  # read pseudobulk counts
  tmp           <- read.table(my.cts.files[i], header = T, sep = "\t")
  rownames(tmp) <- tmp$GENE
  # colnames is X + study number
  
  # select healthy samples
  cts.list[[i]] <- tmp[,colnames(tmp) %in% paste0("X",my.meta.data.flt$Study_ID)] ## some subjects in meta data not in count matrix
  
}

unlist(lapply(cts.list, ncol))
### 
# Astrocytes     GABAergic_Neurons Glutamatergic_Neurons             Microglia      Oligodendrocytes       Polydendrocytes 
# 93                    93                    93                    93                    93                    93 
#### some individuals are missing

table(unlist(lapply(cts.list, colnames)))
### 6 for all => same donors everwhere
# will filter metadata

my.meta.data.flt.v2          <- my.meta.data.flt[paste0("X",my.meta.data.flt$Study_ID) %in% colnames(cts.list[[1]]), ] ### 93
my.meta.data.flt.v2$Study_ID <- paste0("X",my.meta.data.flt.v2$Study_ID) 

#### save counts
save(cts.list, my.meta.data.flt.v2, file = paste0(Sys.Date(),"_McCarroll_PB_objects_QC_healthy.RData"))
#############################################################################################

##############################################################################################
# 2. Use SVA to clean up batch effects on expression and DEseq2 for DE analysis

# will run SVA to clean up batch/noise

######################################################
#######   DEG analysis   ++++   Lu     #######
######################################################

# filter metadata and order it

my.meta           <- my.meta.data.flt.v2[,c("Study_ID","Age","Sex" )]
rownames(my.meta) <- my.meta$Study_ID
my.meta           <- my.meta[order(my.meta$Age),]
my.meta # 38 to 93 year old

# Create list object to receive SVA normalized counts
sva.cts        <- vector(mode = "list", length = length(cts.list))
names(sva.cts) <- names(cts.list)

# Create list object to receive VST normalized counts
vst.cts        <- vector(mode = "list", length = length(cts.list))
names(vst.cts) <- names(cts.list)

# Create list object to receive DESeq2 results
deseq.res.list        <- vector(mode = "list", length = length(cts.list))
names(deseq.res.list) <- names(cts.list)

# loop over pseudobulk data
for  (i in 1:length(cts.list)) {
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_Pseudobulk_McCarroll_dlPFC_",names(cts.list)[[i]])
  
  # reorder columns to match metadata
  cts.list[[i]] <- cts.list[[i]][,my.meta$Study_ID]
  
  ###################################
  #######       Run SVA      #######
  
  # build design matrix
  sva.dataDesign = data.frame( row.names = my.meta$Study_ID , 
                               age       = my.meta$Age,
                               sex       = as.factor(my.meta$Sex))
  
  # Set null and alternative models (ignore batch)
  mod1    = model.matrix(~ age + sex, data = sva.dataDesign)
  n.sv.be = num.sv(cts.list[[i]], mod1, method="be") # 2 for Astrocytes
  
  if (n.sv.be > 0) {
    # apply SVAseq algortihm
    my.svseq = svaseq(as.matrix(cts.list[[i]]), mod1, n.sv=n.sv.be, constant = 0.1)
    
    # remove SV, preserve age and sex
    my.clean <- removeBatchEffect(log2(cts.list[[i]] + 0.1), 
                                  covariates = cbind(my.svseq$sv),
                                  design     = mod1)
    
    # delog and round data for DEseq2 processing
    my.filtered.sva <- round(2^my.clean-0.1)
  } else {
    my.filtered.sva <- cts.list[[i]]
    
  }
  
  
  # keep only robustly expressed genes
  sva.cts[[i]] <- my.filtered.sva
  
  # get matrix using age as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = sva.cts[[i]],
                                colData   = my.meta,
                                design    = ~ Age + Sex )
  
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
       main= paste0(names(cts.list)[[i]]," (MDS)"),
       cex=3, col = "grey", pch = 16,
       cex.lab = 1.25,
       cex.axis = 1.25, las = 1)
  text(x,y, my.meta$Age)
  dev.off()
  
  # extract gene significance by DEseq2
  res.age <- results(dds.deseq, name = "Age") # FC per year
  
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
    my.heatmap.title <- paste0(names(cts.list)[[i]], " aging significant (FDR<5%), ", my.num.age, " genes")
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
my.rdata.age <- paste0(Sys.Date(),"_PB_McCarroll_dlPFC_DEseq2_objects.RData")
save(deseq.res.list, deseq.res.list, file = my.rdata.age)

my.vst.age <- paste0(Sys.Date(),"_PB_McCarroll_dlPFC_VST_data_objects.RData")
save(vst.cts, file = my.vst.age)
##############################################################################################

#######################
sink(file = paste0(Sys.Date(),"_McCarroll_dlPFC_processing_session_Info.txt"))
sessionInfo()
sink()
