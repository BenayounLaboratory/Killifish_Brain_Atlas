setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Mouse_Datasets_for_Comparison/Muscat_DESeq2/Hahn_Hippocampus_CaudatePutamen/')
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

# 2024-03-28
# Parse/annotate Allen Brain data for DecoupleR

################################################################################################################################################################
#### 0. read data and metadata, create new Seurat and SCE objects
load('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Mouse_Datasets_for_Comparison/Processed_Objects/Hahn_GSE212576/2024-03-26_hahn_Hip_PC_Brain_SingleCellNet_Annotated_Seurat_object.RData')
hahn.singlets
# An object of class Seurat 
# 58498 features across 75146 samples within 2 assays 
# Active assay: RNA (32285 features, 0 variable features)
# 1 other assay present: SCT
# 2 dimensional reductions calculated: pca, umap

table(hahn.singlets@meta.data$SampleID,hahn.singlets@meta.data$SingleCellNet_Hajdarovic)
#         Astrocyte Ependymocyte Microglia_Macrophage Neuron Oligodendrocyte  OPC Pericyte_EndothelialCell Tanycyte VascularandLeptomeningealCell
# CP_O_1        223            0                   99   6736            3397    1                       43        1                             0
# CP_O_2         64            0                    7   8784            3628    5                        0        0                             0
# CP_Y_1        343            3                  164   6034            2245    1                       46        4                             3
# CP_Y_2         13            0                  155   7824            2365    0                       15        0                             0
# HIP_O_1      1921            0                   53   3910            1264    1                       10       22                             8
# HIP_O_2       680            5                  206   3304            1127    1                       35        9                            13
# HIP_Y_1      1532            0                  234   5305            1097    0                      117       20                            85
# HIP_Y_2       368            0                  348   9623            1424    3                      113       17                            88

# bring RNA as main assay for processing
DefaultAssay(allen.singlets) <- "RNA"

# split by tissues
hahn.list <- SplitObject(hahn.singlets, split.by = "Tissue")

# convert to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
hahn.hippo.sce <- as.SingleCellExperiment(hahn.list$Hippocampus)
hahn.cp.sce    <- as.SingleCellExperiment(hahn.list$Caudate_Putamen)

save(hahn.hippo.sce, hahn.cp.sce, file = paste(Sys.Date(),"Hahn_GSE212576_SingleCellExperimnent_objects.RData",sep = "_"))
###############################################################################################


###############################################################################################
# 1. Run muscat for pseudobulking and extraction of samples

############################################################
####### Data preparation   ++++   Hahn Hippocampus   #######
############################################################

hahn.hippo.sce.cl <- prepSCE(hahn.hippo.sce, 
                             kid    = "SingleCellNet_Hajdarovic"  ,  # population assignments
                             gid    = "Age_Group"                 ,  # group IDs (ctrl/stim)
                             sid    = "SampleID"                  ,  # sample IDs (ctrl/stim.1234)
                             drop   = TRUE                        )  # drop all other colData columns

# store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk  <- length(kids <- levels(hahn.hippo.sce.cl$cluster_id))
ns  <- length(sids <- levels(hahn.hippo.sce.cl$sample_id))
names(kids) <- kids; names(sids) <- sids

# Aggregation of single-cell to pseudobulk data
pb <- aggregateData(hahn.hippo.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# one list item per cell type
assayNames(pb)
# [1] "Astrocyte"                     "Ependymocyte"                  "Microglia_Macrophage"          "Neuron"                       
# [5] "Oligodendrocyte"               "OPC"                           "Pericyte_EndothelialCell"      "Tanycyte"                     
# [9] "VascularandLeptomeningealCell"

# Number of cells in each sample and cell type
cell.per.samp.tab <- t(table(hahn.hippo.sce.cl$cluster_id, hahn.hippo.sce.cl$sample_id))
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
# [1] "Astrocyte"                "Microglia_Macrophage"     "Neuron"                   "Oligodendrocyte"          "Pericyte_EndothelialCell"

# extract pseudobulk information for samples that pass the cell number cutoff
counts.pb.hippo <- counts.pb.tmp[celltype.qc]


#################################################################
####### Data preparation   ++++   Hahn Caudate Putamen   ########
#################################################################

hahn.cp.sce.cl <- prepSCE(hahn.cp.sce, 
                          kid    = "SingleCellNet_Hajdarovic"  ,  # population assignments
                          gid    = "Age_Group"                 ,  # group IDs (ctrl/stim)
                          sid    = "SampleID"                  ,  # sample IDs (ctrl/stim.1234)
                          drop   = TRUE                        )  # drop all other colData columns

# store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk  <- length(kids <- levels(hahn.cp.sce.cl$cluster_id))
ns  <- length(sids <- levels(hahn.cp.sce.cl$sample_id))
names(kids) <- kids; names(sids) <- sids

# Aggregation of single-cell to pseudobulk data
pb <- aggregateData(hahn.cp.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# one list item per cell type
assayNames(pb)
# [1] "Astrocyte"                     "Ependymocyte"                  "Microglia_Macrophage"          "Neuron"                       
# [5] "Oligodendrocyte"               "OPC"                           "Pericyte_EndothelialCell"      "Tanycyte"                     
# [9] "VascularandLeptomeningealCell"

# Number of cells in each sample and cell type
cell.per.samp.tab <- t(table(hahn.cp.sce.cl$cluster_id, hahn.cp.sce.cl$sample_id))
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
# [1] "Astrocyte"       "Neuron"          "Oligodendrocyte"

# extract pseudobulk information for samples that pass the cell number cutoff
counts.pb.cp <- counts.pb.tmp[celltype.qc]



#### save counts
save(counts.pb.hippo, counts.pb.cp, file = paste0(Sys.Date(),"_Hahn_GSE212576_muscat_PB_objects_QC_Clean.RData"))
#############################################################################################

##############################################################################################
# 2. Use SVA to clean up batch effects on expression and DEseq2 for DE analysis

# clean up memory and reload only muscat PBs
load('2024-03-28_Hahn_GSE212576_muscat_PB_objects_QC_Clean.RData')

# will run SVA to clean up noise

######################################################
#######   DEG analysis   ++++   Hahn hippo     #######
######################################################

# import metadata and order it

colnames(counts.pb.hippo[[1]])
#  "HIP_O_1" "HIP_O_2" "HIP_Y_1" "HIP_Y_2"

my.meta           <- data.frame(matrix(0,4,4))
colnames(my.meta) <- c("SampleID","Age_Group", "Age_month", "Sex")
my.meta$Age_Group <- factor(c(rep("Y",2),rep("O",2)) , levels = c("Y","O"))
my.meta$Age_month <- c(rep(3,2),rep(21,2))
my.meta$Sex       <- rep("Mix",4)
my.meta$SampleID  <- c("HIP_Y_1", "HIP_Y_2", "HIP_O_1", "HIP_O_2")
rownames(my.meta) <- my.meta$SampleID
my.meta

# reorder count tables in sensical order
for  (i in 1:length(counts.pb.hippo)) {
  counts.pb.hippo[[i]] <- counts.pb.hippo[[i]][,my.meta$SampleID]
}


counts.pb.hippo <- counts.pb.hippo[sort(names(counts.pb.hippo))]

# Create list object to receive SVA normalized counts
sva.cts        <- vector(mode = "list", length = length(counts.pb.hippo))
names(sva.cts) <- names(counts.pb.hippo)

# Create list object to receive VST normalized counts
vst.cts        <- vector(mode = "list", length = length(counts.pb.hippo))
names(vst.cts) <- names(counts.pb.hippo)

# Create list object to receive DESeq2 results
deseq.res.list        <- vector(mode = "list", length = length(counts.pb.hippo))
names(deseq.res.list) <- names(counts.pb.hippo)

# loop over pseudobulk data
for  (i in 1:length(counts.pb.hippo)) {
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_Pseudobulk_Hahn_GSE212576_Hippocampus_",names(counts.pb.hippo)[[i]])
  
  ###################################
  #######       Run SVA      #######
  
  # build design matrix
  sva.dataDesign = data.frame( row.names = my.meta$SampleID , 
                               age       = my.meta$Age_month)
  
  # Set null and alternative models (ignore batch)
  mod1    = model.matrix(~ age , data = sva.dataDesign)
  n.sv.be = num.sv(counts.pb.hippo[[i]], mod1, method="be") 
  
  if (n.sv.be >0) {
    # apply SVAseq algortihm
    my.svseq = svaseq(as.matrix(counts.pb.hippo[[i]]), mod1, n.sv=n.sv.be, constant = 0.1)
    
    # remove RIN and SV, preserve age and sex
    my.clean <- removeBatchEffect(log2(counts.pb.hippo[[i]] + 0.1), 
                                  covariates = cbind(my.svseq$sv),
                                  design     = mod1)
    
    # delog and round data for DEseq2 processing
    my.filtered.sva <- round(2^my.clean-0.1)
  } else {
    my.filtered.sva <- counts.pb.hippo[[i]]
    
  }

  
  # keep only robustly expressed genes
  sva.cts[[i]] <- my.filtered.sva
  
  # legend
  my.cols  <- rep("",nrow(my.meta))
  my.cols[my.meta$Age_Group %in% "Y"] <- "purple"
  my.cols[my.meta$Age_Group %in% "O"] <- "purple4"
  
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
       main= paste0(names(counts.pb.hippo)[[i]]," (MDS)"),
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
    my.heatmap.title <- paste0(names(counts.pb.hippo)[[i]], " aging significant (FDR<5%), ", my.num.age, " genes")
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
my.rdata.age <- paste0(Sys.Date(),"_PB_Hahn_GSE212576_Hippocampus_DEseq2_objects.RData")
save(deseq.res.list, deseq.res.list, file = my.rdata.age)

my.vst.age <- paste0(Sys.Date(),"_PB_Hahn_GSE212576_Hippocampus_VST_data_objects.RData")
save(vst.cts, file = my.vst.age)


######################################################
#######   DEG analysis   ++++   Hahn Caudate     #######
######################################################

# import metadata and order it

colnames(counts.pb.cp[[1]])
#  "CP_O_1" "CP_O_2" "CP_Y_1" "CP_Y_2"

my.meta           <- data.frame(matrix(0,4,4))
colnames(my.meta) <- c("SampleID","Age_Group", "Age_month", "Sex")
my.meta$Age_Group <- factor(c(rep("Y",2),rep("O",2)) , levels = c("Y","O"))
my.meta$Age_month <- c(rep(3,2),rep(21,2))
my.meta$Sex       <- rep("Mix",4)
my.meta$SampleID  <- c("CP_Y_1", "CP_Y_2", "CP_O_1", "CP_O_2")
rownames(my.meta) <- my.meta$SampleID
my.meta

# reorder count tables in sensical order
for  (i in 1:length(counts.pb.cp)) {
  counts.pb.cp[[i]] <- counts.pb.cp[[i]][,my.meta$SampleID]
}


counts.pb.cp <- counts.pb.cp[sort(names(counts.pb.cp))]

# Create list object to receive SVA normalized counts
sva.cts        <- vector(mode = "list", length = length(counts.pb.cp))
names(sva.cts) <- names(counts.pb.cp)

# Create list object to receive VST normalized counts
vst.cts        <- vector(mode = "list", length = length(counts.pb.cp))
names(vst.cts) <- names(counts.pb.cp)

# Create list object to receive DESeq2 results
deseq.res.list        <- vector(mode = "list", length = length(counts.pb.cp))
names(deseq.res.list) <- names(counts.pb.cp)

# loop over pseudobulk data
for  (i in 1:length(counts.pb.cp)) {
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_Pseudobulk_Hahn_GSE212576_Caudate_Putamen_",names(counts.pb.cp)[[i]])
  
  ###################################
  #######       Run SVA      #######
  
  # build design matrix
  sva.dataDesign = data.frame( row.names = my.meta$SampleID , 
                               age       = my.meta$Age_month)
  
  # Set null and alternative models (ignore batch)
  mod1    = model.matrix(~ age , data = sva.dataDesign)
  n.sv.be = num.sv(counts.pb.cp[[i]], mod1, method="be") 
  
  if (n.sv.be >0) {
    # apply SVAseq algortihm
    my.svseq = svaseq(as.matrix(counts.pb.cp[[i]]), mod1, n.sv=n.sv.be, constant = 0.1)
    
    # remove RIN and SV, preserve age and sex
    my.clean <- removeBatchEffect(log2(counts.pb.cp[[i]] + 0.1), 
                                  covariates = cbind(my.svseq$sv),
                                  design     = mod1)
    
    # delog and round data for DEseq2 processing
    my.filtered.sva <- round(2^my.clean-0.1)
  } else {
    my.filtered.sva <- counts.pb.cp[[i]]
    
  }
  
  
  # keep only robustly expressed genes
  sva.cts[[i]] <- my.filtered.sva
  
  # legend
  my.cols  <- rep("",nrow(my.meta))
  my.cols[my.meta$Age_Group %in% "Y"] <- "purple"
  my.cols[my.meta$Age_Group %in% "O"] <- "purple4"
  
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
       main= paste0(names(counts.pb.cp)[[i]]," (MDS)"),
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
    my.heatmap.title <- paste0(names(counts.pb.cp)[[i]], " aging significant (FDR<5%), ", my.num.age, " genes")
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
my.rdata.age <- paste0(Sys.Date(),"_PB_Hahn_GSE212576_Caudate_Putamen_DEseq2_objects.RData")
save(deseq.res.list, deseq.res.list, file = my.rdata.age)

my.vst.age <- paste0(Sys.Date(),"_PB_Hahn_GSE212576_Caudate_Putamen_VST_data_objects.RData")
save(vst.cts, file = my.vst.age)
##############################################################################################

#######################
sink(file = paste0(Sys.Date(),"_Hahn_GSE212576_processing_session_Info.txt"))
sessionInfo()
sink()
