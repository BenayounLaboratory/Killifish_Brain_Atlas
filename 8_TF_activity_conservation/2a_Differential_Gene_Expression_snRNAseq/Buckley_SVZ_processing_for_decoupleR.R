setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Mouse_Datasets_for_Comparison/Muscat_DESeq2')
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

# 2024-03-22
# Parse/annotate Buckley data for DecoupleR

################################################################################################################################################################
#### 0. read data and metadata, create new Seurat and SCE objects
svz.data <- readRDS('./Processed_Objects/Buckley_multi_intergrated_seurat_Dec2020.rds')
svz.data
# An object of class Seurat 
# 49678 features across 21458 samples within 3 assays 
# Active assay: SCT (18597 features, 3000 variable features)
# 2 other assays present: RNA, LMO
# 4 dimensional reductions calculated: pca, umap, harmony, umap_har

table(svz.data@meta.data$Celltype.LowRes,svz.data@meta.data$Age)
#                3.3 3.33 3.6 4.3 4.7 5.4 6.7 8.4 9.47 10.43 12.47 14.5 14.77 16.53 16.83 18.58 18.87 20.6 20.8 21.57 22.57 22.6 23.9 24.9 25.93  29
# Microglia      263   43  88 328  47 186 159 173  217   268   308  133    87    21   139    34   153    7  226   150    15  125  106  161   123 284
# Oligodendro    405   82 122 313 295 348 468 285  143   521   361  255    81    31   130    97   163   24  520   193    52  217   74  166   145 452
# Neuroblast     626   88 178 540 112 177 107 325   90   436   515  213    91    20   102    21   114    9   44   176    16  100   42   91    83  25
# Astrocyte_qNSC  91  115  39 156 260 215 349  95   57   143   172  136    63    29    67    76    72   49  145    50    69   65   34   78    50  81
# aNSC_NPC       350   66 155 315  91 100  91 165   60   183   205   63    67     5   100    12    79    3   36    74     6   40   22   25    53  30
# Endothelial     44   44  13  31  17 115  36  28  153    48    57   63     9    11    29    22    28    7   44    25     6   25    9   19    10  46
# Mural           40    9   6  12   8  36  18  22   67    26    39   28     4     8     6    10    12    6   28     3     5   17    0    6     3  15
# OPC             39   20  12  20  41  17  23   7    5     5     7    1     3     0     2     8     7    5   19     4     1    4    3    6     6   9
# Macrophage       9   13  14  11  19  34  22  10   14    15    16    5     5     5     3     4     6    4   22     6     1    7    1    7     9   6
# Neuron          11    5   2   8   2   6   2  12   22    17    13    6     4     1     3     1     7    1    6     2     1   12    3    5     2   2
# Ependymal        5   10   0   3  16  12  10   1    6     5     1    7     3     0     2     0     3    1    8     1     5    1    1    1     1   4
# Doublet          0    0   0   0   0   0   0   0    0     0     0    0     0     0     0     0     0    0    0     0     0    0    0    0     0   0

# bring RNA as main assay fro processing
DefaultAssay(svz.data) <- "RNA"

table(svz.data@meta.data$Age)
#  3.3  3.33   3.6   4.3   4.7   5.4   6.7   8.4  9.47 10.43 12.47  14.5 14.77 16.53 16.83 18.58 18.87  20.6  20.8 21.57 22.57  22.6  23.9  24.9 25.93 
# 1883   495   629  1737   908  1246  1285  1123   834  1667  1694   910   417   131   583   285   644   116  1098   684   177   613   295   565   485 
# 29 
# 954 

table(svz.data@meta.data$Age) > 200
#  3.3  3.33   3.6   4.3   4.7   5.4   6.7   8.4  9.47 10.43 12.47  14.5 14.77 16.53 16.83 18.58 18.87  20.6  20.8 21.57 22.57  22.6  23.9  24.9 25.93    29 
# TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE 

# Retain samples with at least 200 cells for downstream analysis
svz.data.cl <- subset(svz.data, subset = Age %in% c(16.53, 20.6, 22.57), invert = TRUE) # 21,034 cells
svz.data.cl

# convert to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
svz.data.sce <- as.SingleCellExperiment(svz.data.cl)
save(svz.data.sce, file = paste(Sys.Date(),"Buckley_SingleCellExperimnent_object_Clean.RData",sep = "_"))
###############################################################################################


###############################################################################################
# 1. Run muscat for pseudobulking and extraction of samples

######################################################
####### Data preparation   ++++   Buckley   #######
######################################################

svz.data.sce.cl <- prepSCE(svz.data.sce, 
                            kid    = "Celltype.LowRes",  # population assignments
                            gid    = "Phase"          ,  # group IDs (ctrl/stim)
                            sid    = "Age"            ,  # sample IDs (ctrl/stim.1234)
                            drop   = TRUE             )  # drop all other colData columns

# store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk  <- length(kids <- levels(svz.data.sce.cl$cluster_id))
ns  <- length(sids <- levels(svz.data.sce.cl$sample_id))
names(kids) <- kids; names(sids) <- sids

# Aggregation of single-cell to pseudobulk data
pb <- aggregateData(svz.data.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# one list item per cell type
assayNames(pb)
# [1] "Microglia"      "Oligodendro"    "Neuroblast"     "Astrocyte_qNSC" "aNSC_NPC"       "Endothelial"    "Mural"          "OPC"           
# [9] "Macrophage"     "Neuron"         "Ependymal"     

# Number of cells in each sample and cell type
cell.per.samp.tab <- t(table(svz.data.sce.cl$cluster_id, svz.data.sce.cl$sample_id))
#       Microglia Oligodendro Neuroblast Astrocyte_qNSC aNSC_NPC Endothelial Mural OPC Macrophage Neuron Ependymal Doublet
# 3.3         263         405        626             91      350          44    40  39          9     11         5       0
# 3.33         43          82         88            115       66          44     9  20         13      5        10       0
# 3.6          88         122        178             39      155          13     6  12         14      2         0       0
# 4.3         328         313        540            156      315          31    12  20         11      8         3       0
# 4.7          47         295        112            260       91          17     8  41         19      2        16       0
# 5.4         186         348        177            215      100         115    36  17         34      6        12       0
# 6.7         159         468        107            349       91          36    18  23         22      2        10       0
# 8.4         173         285        325             95      165          28    22   7         10     12         1       0
# 9.47        217         143         90             57       60         153    67   5         14     22         6       0
# 10.43       268         521        436            143      183          48    26   5         15     17         5       0
# 12.47       308         361        515            172      205          57    39   7         16     13         1       0
# 14.5        133         255        213            136       63          63    28   1          5      6         7       0
# 14.77        87          81         91             63       67           9     4   3          5      4         3       0
# 16.83       139         130        102             67      100          29     6   2          3      3         2       0
# 18.58        34          97         21             76       12          22    10   8          4      1         0       0
# 18.87       153         163        114             72       79          28    12   7          6      7         3       0
# 20.8        226         520         44            145       36          44    28  19         22      6         8       0
# 21.57       150         193        176             50       74          25     3   4          6      2         1       0
# 22.6        125         217        100             65       40          25    17   4          7     12         1       0
# 23.9        106          74         42             34       22           9     0   3          1      3         1       0
# 24.9        161         166         91             78       25          19     6   6          7      5         1       0
# 25.93       123         145         83             50       53          10     3   6          9      2         1       0
# 29          284         452         25             81       30          46    15   9          6      2         4       0


# extract pseudobulk information
counts.pb.tmp <- pb@assays@data

# get the genes with no reads in at least half the samples out, they mess up the algorithm
for (i in 1:length(counts.pb.tmp)) {
  my.good <- which(apply(counts.pb.tmp[[i]]>0, 1, sum) >= nrow(cell.per.samp.tab)/2) # see deseq2 vignette, need to remove too low genes
  counts.pb.tmp[[i]] <- counts.pb.tmp[[i]][my.good,]
}

# cell types with at least 10 cells in all samples (20 samples for GRZ)
celltype.qc <- colnames(cell.per.samp.tab)[colSums(cell.per.samp.tab  >= 10) == nrow(cell.per.samp.tab)]
celltype.qc
#  "Microglia"      "Oligodendro"    "Neuroblast"     "Astrocyte_qNSC" "aNSC_NPC"   

# extract pseudobulk information for samples that pass the cell number cutoff
counts.pb <- counts.pb.tmp[celltype.qc]

#### save counts
save(counts.pb, file = paste0(Sys.Date(),"_Buckley_muscat_PB_objects_QC_Clean.RData"))
#############################################################################################

##############################################################################################
# 2. Use SVA to clean up batch effects on expression and DEseq2 for DE analysis

# clean up memory and reload only muscat PBs
load('2024-03-25_Buckley_muscat_PB_objects_QC_Clean.RData')

# will run SVA to clean up noise

######################################################
#######   DEG analysis   ++++   Buckley     #######
######################################################

# import metadata and order it

colnames(counts.pb[[1]])
# [1] "3.3"   "3.33"  "3.6"   "4.3"   "4.7"   "5.4"   "6.7"   "8.4"   "9.47"  "10.43" "12.47" "14.5"  "14.77" "16.83" "18.58" "18.87" "20.8"  "21.57" "22.6" 
# [20] "23.9"  "24.9"  "25.93" "29"   

my.meta           <- data.frame(matrix(0, ncol(counts.pb[[1]]), 3))
colnames(my.meta) <- c("SampleID","Age_month", "Sex")
my.meta$Age_month <- as.numeric(colnames(counts.pb[[1]]))
my.meta$Sex       <- rep("M",ncol(counts.pb[[1]]))
my.meta$SampleID  <- paste0("Mouse_",colnames(counts.pb[[1]]))
rownames(my.meta) <- my.meta$SampleID
my.meta

# rename samples to not be numbers
for  (i in 1:length(counts.pb)) {
  colnames(counts.pb[[i]]) <- paste0("Mouse_",colnames(counts.pb[[i]]))
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
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_Pseudobulk_Buckley_",names(counts.pb)[[i]])
  
  ###################################
  #######       Run SVA      #######
  
  # build design matrix
  sva.dataDesign = data.frame( row.names = my.meta$SampleID , 
                               age       = my.meta$Age_month)
  
  # Set null and alternative models (ignore batch)
  mod1    = model.matrix(~ age , data = sva.dataDesign)
  n.sv.be = num.sv(counts.pb[[i]], mod1, method="be") 
  
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
  
  # legend
  my.cols  <-  colorRampPalette(c("deepskyblue","deepskyblue4"))(nrow(my.meta))
  
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
  
  # # output result tables of combined analysis to text files
  # my.out.ct.mat <- paste0(my.outprefix,"_AGING_VST_log2_counts_matrix.txt")
  # write.table(vst.cts[[i]], file = my.out.ct.mat , sep = "\t" , row.names = T, quote = F)
  # 
  # my.out.stats.age <- paste0(my.outprefix,"_AGING_all_genes_statistics.txt")
  # write.table(deseq.res.list[[i]], file = my.out.stats.age , sep = "\t" , row.names = T, quote = F)
  
}

# save R object with all DEseq2 results
my.rdata.age <- paste0(Sys.Date(),"_PB_Buckley_CellType_DEseq2_objects.RData")
save(deseq.res.list, deseq.res.list, file = my.rdata.age)

my.vst.age <- paste0(Sys.Date(),"_PB_Buckley_CellType_VST_data_objects.RData")
save(vst.cts, file = my.vst.age)
##############################################################################################

#######################
sink(file = paste0(Sys.Date(),"_Buckley_processing_session_Info.txt"))
sessionInfo()
sink()
