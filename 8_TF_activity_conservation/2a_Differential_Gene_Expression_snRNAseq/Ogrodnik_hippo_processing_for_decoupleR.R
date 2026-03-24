setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Mouse_Datasets_for_Comparison/Muscat_DESeq2/Ogrodnik_GSE161340')
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

# 2024-04-12
# Parse/annotate Ogrodnik Hippocampus data for DecoupleR

################################################################################################################################################################
#### 0. read data and metadata, create new Seurat and SCE objects
load('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Mouse_Datasets_for_Comparison/Processed_Objects/Ogrodnik_GSE161340/2024-04-12_Ogrodnik_Hippocampus_SingleCellNet_Annotated_Seurat_object.RData')
Ogrodnik.singlets
# An object of class Seurat 
# 52162 features across 16059 samples within 2 assays 
# Active assay: RNA (31053 features, 0 variable features)
# 1 other assay present: SCT
# 2 dimensional reductions calculated: pca, umap

table(Ogrodnik.singlets@meta.data$SampleID,Ogrodnik.singlets@meta.data$SingleCellNet_Hajdarovic)
#         Astrocyte Ependymocyte Microglia_Macrophage Neuron Oligodendrocyte  OPC Pericyte_EndothelialCell Tanycyte VascularandLeptomeningealCell
# HIP_O_1       287            5                  388   2298            1163   48                       15       36                            21
# HIP_O_2       551           27                  433   3137            1425   53                       18      102                            20
# HIP_Y_1       322            7                  190   1881             677   23                        1       22                            12
# HIP_Y_2       273            4                  178   1841             521   29                       14       21                            16

# bring RNA as main assay fpr processing
DefaultAssay(Ogrodnik.singlets) <- "RNA"

# convert to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
Ogrodnik.singlets.sce <- as.SingleCellExperiment(Ogrodnik.singlets)
save(Ogrodnik.singlets.sce, file = paste(Sys.Date(),"Ogrodnik_GSE161340_SingleCellExperimnent_object.RData",sep = "_"))
###############################################################################################


###############################################################################################
# 1. Run muscat for pseudobulking and extraction of samples

######################################################
####### Data preparation   ++++   Ogrodnik  PFC   #######
######################################################

Ogrodnik.singlets.sce.cl <- prepSCE(Ogrodnik.singlets.sce, 
                                 kid    = "SingleCellNet_Hajdarovic"  ,  # population assignments
                                 gid    = "Age_Group"                 ,  # group IDs (ctrl/stim)
                                 sid    = "SampleID"                  ,  # sample IDs (ctrl/stim.1234)
                                 drop   = TRUE                        )  # drop all other colData columns

# store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk  <- length(kids <- levels(Ogrodnik.singlets.sce.cl$cluster_id))
ns  <- length(sids <- levels(Ogrodnik.singlets.sce.cl$sample_id))
names(kids) <- kids; names(sids) <- sids

# Aggregation of single-cell to pseudobulk data
pb <- aggregateData(Ogrodnik.singlets.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# one list item per cell type
assayNames(pb)
# [1] "Astrocyte"                     "Ependymocyte"                  "Microglia_Macrophage"          "Neuron"                        "Oligodendrocyte"              
# [6] "OPC"                           "Pericyte_EndothelialCell"      "Tanycyte"                      "VascularandLeptomeningealCell"
# 
# Number of cells in each sample and cell type
cell.per.samp.tab <- t(table(Ogrodnik.singlets.sce.cl$cluster_id, Ogrodnik.singlets.sce.cl$sample_id))
cell.per.samp.tab


# extract pseudobulk information
counts.pb.tmp <- pb@assays@data

# get the genes with no reads in at least half the samples out, they mess up the algorithm
for (i in 1:length(counts.pb.tmp)) {
  my.good <- which(apply(counts.pb.tmp[[i]]>0, 1, sum) >= nrow(cell.per.samp.tab)/2) # see deseq2 vignette, need to remove too low genes
  counts.pb.tmp[[i]] <- counts.pb.tmp[[i]][my.good,]
}

# cell types with at least 10 cells in all samples
celltype.qc <- colnames(cell.per.samp.tab)[colSums(cell.per.samp.tab  >= 10) == nrow(cell.per.samp.tab)]
celltype.qc
# [1] "Astrocyte"                     "Microglia_Macrophage"          "Neuron"                        "Oligodendrocyte"               "OPC"                          
# [6] "Tanycyte"                      "VascularandLeptomeningealCell"

# extract pseudobulk information for samples that pass the cell number cutoff
counts.pb <- counts.pb.tmp[celltype.qc]

#### save counts
save(counts.pb, file = paste0(Sys.Date(),"_Ogrodnik_GSE161340_muscat_PB_objects_QC_Clean.RData"))
#############################################################################################

##############################################################################################
# 2. Use SVA to clean up batch effects on expression and DEseq2 for DE analysis

# clean up memory and reload only muscat PBs
load('2024-04-12_Ogrodnik_GSE161340_muscat_PB_objects_QC_Clean.RData')

# will run SVA to clean up noise

######################################################
#######   DEG analysis   ++++   Hajdarovic     #######
######################################################

# import metadata and order it

colnames(counts.pb[[1]])
# [1] "HIP_O_1" "HIP_O_2" "HIP_Y_1" "HIP_Y_2"

my.meta           <- data.frame(matrix(0,4,4))
colnames(my.meta) <- c("SampleID","Age_Group", "Age_month", "Sex")
my.meta$Age_Group <- factor(c(rep("Y",2),rep("O",2)) , levels = c("Y","O"))
my.meta$Age_month <- c(rep(4,2),rep(25,2))
my.meta$Sex       <- rep("F",4)
my.meta$SampleID  <- c("HIP_Y_1", "HIP_Y_2", "HIP_O_1", "HIP_O_2")
rownames(my.meta) <- my.meta$SampleID
my.meta

# reorder count tables in sensical order
for  (i in 1:length(counts.pb)) {
  counts.pb[[i]] <- counts.pb[[i]][,my.meta$SampleID]
}


counts.pb <- counts.pb[sort(names(counts.pb))]

# Create list object to receive SVA normalized counts
sva.cts        <- vector(mode = "list", length = length(counts.pb))
names(sva.cts) <- names(counts.pb)

# Create list object to receive VST normalized counts
vst.cts        <- vector(mode = "list", length = length(counts.pb))
names(vst.cts) <- names(counts.pb)

# Create list object to receive DESeq2 results
deseq.res.list        <- vector(mode = "list", length = length(counts.pb))
names(deseq.res.list) <- names(counts.pb)

# loop over pseudobulk data
for  (i in 1:length(counts.pb)) {
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_Pseudobulk_Ogrodnik_GSE161340_",names(counts.pb)[[i]])
  
  ###################################
  #######       Run SVA      #######
  
  # build design matrix
  sva.dataDesign = data.frame( row.names = my.meta$SampleID , 
                               age       = my.meta$Age_month)
  
  # Set null and alternative models (ignore batch)
  mod1    = model.matrix(~ age , data = sva.dataDesign)
  n.sv.be = num.sv(counts.pb[[i]], mod1, method="be") 
  
  if (n.sv.be > 0) {
    # apply SVAseq algortihm
    my.svseq = svaseq(as.matrix(counts.pb[[i]]), mod1, n.sv=n.sv.be, constant = 0.1)
    
    # remove RIN and SV, preserve age and sex
    my.clean <- removeBatchEffect(log2(counts.pb[[i]] + 0.1), 
                                  covariates = cbind(my.svseq$sv),
                                  design     = mod1)
    
    # delog and round data for DEseq2 processing
    my.filtered.sva <- round(2^my.clean-0.1)
    
    # keep only robustly expressed genes
    sva.cts[[i]] <- my.filtered.sva
  } else {
    sva.cts[[i]] <- counts.pb[[i]]
  }

  
  # legend
  my.cols  <- rep("",nrow(my.meta))
  my.cols[my.meta$Age_Group %in% "Y"] <- "deeppink"
  my.cols[my.meta$Age_Group %in% "O"] <- "deeppink4"
  
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
       main= paste0(names(counts.pb)[[i]]," (MDS)"),
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
    my.heatmap.title <- paste0(names(counts.pb)[[i]], " aging significant (FDR<5%), ", my.num.age, " genes")
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
my.rdata.age <- paste0(Sys.Date(),"_PB_Ogrodnik_Hippocampus_DEseq2_objects.RData")
save(deseq.res.list, deseq.res.list, file = my.rdata.age)

my.vst.age <- paste0(Sys.Date(),"_PB_Ogrodnik_Hippocampus_VST_data_objects.RData")
save(vst.cts, file = my.vst.age)
##############################################################################################

#######################
sink(file = paste0(Sys.Date(),"_Ogrodnik_GSE161340_processing_session_Info.txt"))
sessionInfo()
sink()
