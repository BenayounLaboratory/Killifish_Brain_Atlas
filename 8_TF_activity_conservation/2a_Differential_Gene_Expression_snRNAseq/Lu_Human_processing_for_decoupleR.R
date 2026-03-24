setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Species_Comparison/Human_Datasets_for_comparison/snRNA/Muscat_DESeq2')
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

# 2025-07-04
# Parse/annotate Lu Brain data for DecoupleR

################################################################################################################################################################
#### 0. read data and metadata, create new Seurat and SCE objects
load('../GSE212606_RAW_Human/2025-07-04_Lu_Brain_Aging_Seurat_object.RData')
Lu.data.flt
# An object of class Seurat 
# 60706 features across 30801 samples within 1 assay 
# Active assay: RNA (60706 features, 0 variable features)
# 2 layers present: counts, data

table(Lu.data.flt$Region) # Hippocampus

table(Lu.data.flt@meta.data$Cell_type,Lu.data.flt@meta.data$Individual_ID)
#                                  1247 1304 1306 1311 5356 5459
# Astrocytes                        473  333  719  435  409  376
# Cortical projection neurons 1    1826  497 1003 1201 1018 1198
# Cortical projection neurons 2     110  103    3   87   72  100
# Cortical projection neurons 3      27    0    0   17   12   17
# Endothelial cells                 283  107  422  272  155  194
# Ependymal cells                     0    0  203    0    0    0
# Interneurons 1                    244   87   81  271  255  205
# Interneurons 2                    208  153  122  275  239  200
# Microglia                         306  668  807  485  385  473
# Oligodendrocyte progenitor cells  181  318  209  254  223  227
# Oligodendrocytes                 2621 1860 2938 1753 1851 1116
# Vascular leptomeningeal cells       0    0    0    0  113    1

# bring RNA as main assay fro processing
DefaultAssay(Lu.data.flt) <- "RNA"

# convert to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
Lu.data.flt.sce <- as.SingleCellExperiment(Lu.data.flt)
save(Lu.data.flt.sce, file = paste(Sys.Date(),"Lu_GSM6657986_SingleCellExperimnent_object.RData",sep = "_"))
###############################################################################################


###############################################################################################
# 1. Run muscat for pseudobulking and extraction of samples

######################################################
####### Data preparation   ++++   Lu   #######
######################################################

Lu.data.flt.sce.cl <- prepSCE(Lu.data.flt.sce, 
                            kid    = "Cell_type"     ,  # population assignments
                            gid    = "Condition"     ,  # group IDs (ctrl/stim)
                            sid    = "Individual_ID" ,  # sample IDs (ctrl/stim.1234)
                            drop   = TRUE           )  # drop all other colData columns

# store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk  <- length(kids <- levels(Lu.data.flt.sce.cl$cluster_id))
ns  <- length(sids <- levels(Lu.data.flt.sce.cl$sample_id))
names(kids) <- kids; names(sids) <- sids

# Aggregation of single-cell to pseudobulk data
pb <- aggregateData(Lu.data.flt.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# one list item per cell type
assayNames(pb)
# [1] "Astrocytes"                       "Cortical projection neurons 1"    "Cortical projection neurons 2"    "Cortical projection neurons 3"   
# [5] "Endothelial cells"                "Ependymal cells"                  "Interneurons 1"                   "Interneurons 2"                  
# [9] "Microglia"                        "Oligodendrocyte progenitor cells" "Oligodendrocytes"                 "Vascular leptomeningeal cells"   

# Number of cells in each sample and cell type
cell.per.samp.tab <- t(table(Lu.data.flt.sce.cl$cluster_id, Lu.data.flt.sce.cl$sample_id))
# Astrocytes Cortical projection neurons 1 Cortical projection neurons 2 Cortical projection neurons 3 Endothelial cells Ependymal cells Interneurons 1
# 1247        473                          1826                           110                            27               283               0            244
# 1304        333                           497                           103                             0               107               0             87
# 1306        719                          1003                             3                             0               422             203             81
# 1311        435                          1201                            87                            17               272               0            271
# 5356        409                          1018                            72                            12               155               0            255
# 5459        376                          1198                           100                            17               194               0            205
# 
# Interneurons 2 Microglia Oligodendrocyte progenitor cells Oligodendrocytes Vascular leptomeningeal cells
# 1247            208       306                              181             2621                             0
# 1304            153       668                              318             1860                             0
# 1306            122       807                              209             2938                             0
# 1311            275       485                              254             1753                             0
# 5356            239       385                              223             1851                           113
# 5459            200       473                              227             1116                             1

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
# [1] "Astrocytes"                       "Cortical projection neurons 1"    "Endothelial cells"                "Interneurons 1"                          
# [5] "Interneurons 2"                   "Microglia"                        "Oligodendrocyte progenitor cells" "Oligodendrocytes"         

# extract pseudobulk information for samples that pass the cell number cutoff
counts.pb <- counts.pb.tmp[celltype.qc]

#### save counts
save(counts.pb, file = paste0(Sys.Date(),"_Lu_GSM6657986_muscat_PB_objects_QC_Clean.RData"))
#############################################################################################

##############################################################################################
# 2. Use SVA to clean up batch effects on expression and DEseq2 for DE analysis

# clean up memory and reload only muscat PBs

# will run SVA to clean up noise

######################################################
#######   DEG analysis   ++++   Lu     #######
######################################################

# import metadata and order it

colnames(counts.pb[[1]])
# "1247" "1304" "1306" "1311" "5356" "5459"

my.meta           <- unique(Lu.data.flt@meta.data[,c("Individual_ID","Age","Sex" )])
rownames(my.meta) <- my.meta$Individual_ID
my.meta           <- my.meta[order(my.meta$Age),]
my.meta

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
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_Pseudobulk_Lu_GSM6657986_",names(counts.pb)[[i]])
  
  ###################################
  #######       Run SVA      #######
  
  # build design matrix
  sva.dataDesign = data.frame( row.names = my.meta$Individual_ID , 
                               age       = my.meta$Age,
                               sex       = as.factor(my.meta$Sex))
  
  # Set null and alternative models (ignore batch)
  mod1    = model.matrix(~ age + sex, data = sva.dataDesign)
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
  } else {
    my.filtered.sva <- counts.pb[[i]]
    
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
       main= paste0(names(counts.pb)[[i]]," (MDS)"),
       cex=3, col = "grey", pch = 16,
       cex.lab = 1.25,
       cex.axis = 1.25, las = 1)
  text(x,y, my.meta$Age)
  dev.off()
  
  # extract gene significance by DEseq2
  res.age <- results(dds.deseq, name = "Age") # FC per month
  
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
my.rdata.age <- paste0(Sys.Date(),"_PB_Lu_CellType_GSM6657986_DEseq2_objects.RData")
save(deseq.res.list, deseq.res.list, file = my.rdata.age)

my.vst.age <- paste0(Sys.Date(),"_PB_Lu_CellType_GSM6657986_VST_data_objects.RData")
save(vst.cts, file = my.vst.age)
##############################################################################################

#######################
sink(file = paste0(Sys.Date(),"_Lu_GSM6657986_processing_session_Info.txt"))
sessionInfo()
sink()
