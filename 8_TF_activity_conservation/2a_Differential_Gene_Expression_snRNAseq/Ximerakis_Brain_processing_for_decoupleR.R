setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Mouse_Datasets_for_Comparison/Muscat_DESeq2/Ximerakis_Brain')
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

# 2024-04-11
# Parse/annotate Ximerakis Brain data for DecoupleR

################################################################################################################################################################
#### 0. read data and metadata, create new Seurat and SCE objects
load('../../Processed_Objects/Ximerakis/2024-04-09_Ximerakis_Brain_Aging_object.RData')
ximerakis.seurat
# An object of class Seurat 
# 29397 features across 37069 samples within 2 assays 
# Active assay: SCT (14698 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

table(ximerakis.seurat@meta.data$cell_type_full)
# Arachnoid_barrier_cells      Astrocyte_restricted_precursors                           Astrocytes      Choroid_plexus_epithelial_cells 
# 307                                  184                                 6747                                   89 
# Dendritic_cells                    Endothelial_cells                        Ependymocytes Hemoglobin_expressing_vascular_cells 
# 55                                 2413                                  274                                   81 
# Hypendymal_cells                     immature_Neurons                          Macrophages                       mature_Neurons 
# 12                                  162                                  377                                 5135 
# Microglia                            Monocytes                    Neural_stem_cells                 Neuroendocrine_cells 
# 3910                                   77                                  166                                  394 
# Neuronal_restricted_precursors                          Neutrophils           Olfactory_ensheathing_glia      Oligodendrocyte_precursor_cells 
# 82                                   29                                  892                                 2187 
# Oligodendrocytes                            Pericytes                            Tanycytes    Vascular_and_leptomeningeal_cells 
# 12384                                  735                                   29                                  105 
# Vascular_smooth_muscle_cells 
# 243 

# Simplify
ximerakis.seurat@meta.data$cell_type_simple <- ximerakis.seurat@meta.data$cell_type_full
ximerakis.seurat@meta.data$cell_type_simple[ximerakis.seurat@meta.data$cell_type_full %in% "Macrophages"] <- "Microglia"
ximerakis.seurat@meta.data$cell_type_simple[ximerakis.seurat@meta.data$cell_type_full %in% c("Hemoglobin_expressing_vascular_cells","Vascular_and_leptomeningeal_cells")] <- "Vascular_cells"
ximerakis.seurat@meta.data$cell_type_simple[ximerakis.seurat@meta.data$cell_type_full %in% "Oligodendrocyte_precursor_cells"] <- "OPCs"
ximerakis.seurat@meta.data$cell_type_simple[ximerakis.seurat@meta.data$cell_type_full %in% c("Neuronal_restricted_precursors","Neural_stem_cells")] <- "NSPCs"
ximerakis.seurat@meta.data$cell_type_simple[ximerakis.seurat@meta.data$cell_type_full %in% c("immature_Neurons","mature_Neurons","Neuroendocrine_cells")] <- "Neurons"

table(ximerakis.seurat@meta.data$cell_type_simple)
# Arachnoid_barrier_cells Astrocyte_restricted_precursors                      Astrocytes Choroid_plexus_epithelial_cells                 Dendritic_cells 
# 307                             184                            6747                              89                              55 
# Endothelial_cells                   Ependymocytes                Hypendymal_cells                       Microglia                       Monocytes 
# 2413                             274                              12                            4287                              77 
# Neurons                     Neutrophils                           NSPCs      Olfactory_ensheathing_glia                Oligodendrocytes 
# 5691                              29                             248                             892                           12384 
# OPCs                       Pericytes                       Tanycytes                  Vascular_cells    Vascular_smooth_muscle_cells 
# 2187                             735                              29                             186                             243 

# make an anima/age varaible
ximerakis.seurat@meta.data$Bio_Group <- paste0(ximerakis.seurat@meta.data$Age_Group,"_",ximerakis.seurat@meta.data$animal)

table(ximerakis.seurat@meta.data$Bio_Group,ximerakis.seurat@meta.data$cell_type_simple)
#            Arachnoid_barrier_cells Astrocyte_restricted_precursors Astrocytes Choroid_plexus_epithelial_cells Dendritic_cells Endothelial_cells Ependymocytes
# O_Mouse_33                      30                              32        761                              10               3               223            24
# O_Mouse_34                      36                              31        716                              13               5               222            31
# O_Mouse_37                      18                               6        348                               4               3               189            15
# O_Mouse_38                      29                               7        235                               2               4               126            15
# O_Mouse_39                      18                               7        300                               4               2               147            23
# O_Mouse_40                       5                               8        324                               5               1               124            12
# O_Mouse_43                      16                               8        333                               7               3               185            22
# O_Mouse_44                      33                              14        403                               6               4               186            28
# Y_Mouse_19                      16                               9        445                               9               7               118            16
# Y_Mouse_20                      12                               4        447                               4               3               139            19
# Y_Mouse_21                      26                              12        350                               2               2               115            16
# Y_Mouse_22                      10                              13        347                               2               2               103            10
# Y_Mouse_27                      12                              11        576                               8               5               173            10
# Y_Mouse_28                      29                              12        675                               5               3               196            12
# Y_Mouse_6                        5                               2        125                               5               4                79             4
# Y_Mouse_7                       12                               8        362                               3               4                88            17
# 
#            Hypendymal_cells Microglia Monocytes Neurons Neutrophils NSPCs Olfactory_ensheathing_glia Oligodendrocytes OPCs Pericytes Tanycytes Vascular_cells
# O_Mouse_33                2       466         8     465           0    20                        120             1482  209        71         4             22
# O_Mouse_34                0       416         4     457           2     8                         98             1453  211        54         6             12
# O_Mouse_37                0       178        11     303           0     8                         32              548   72        52         1             11
# O_Mouse_38                0       152         3     353           1     5                         54              402   27        49         0              9
# O_Mouse_39                0       259         6     376           3    17                         67              802   95        30         0              8
# O_Mouse_40                2       261         6     548           1    10                         51              706  115        42         2             10
# O_Mouse_43                1       317         1     394           1     5                         55              981  153        71         2              9
# O_Mouse_44                0       385         4     473           1    15                         36             1123  150        49         2              2
# Y_Mouse_19                2       200         3     264           0    16                         49              505  125        44         3              8
# Y_Mouse_20                1       253         5     328           6    24                         45              748  159        45         3             24
# Y_Mouse_21                0       212         2     277           6    31                         63              541  147        33         2             14
# Y_Mouse_22                0       172         1     380           0    22                         26              561  142        30         1             10
# Y_Mouse_27                1       323         3     317           2    21                         59              714  209        54         0             13
# Y_Mouse_28                0       333         4     253           2    25                        105              711  187        64         0             16
# Y_Mouse_6                 0       154         3     330           3     6                         14              528   74        21         2              3
# Y_Mouse_7                 3       206        13     173           1    15                         18              579  112        26         1             15
# 
#            Vascular_smooth_muscle_cells
# O_Mouse_33                           32
# O_Mouse_34                           28
# O_Mouse_37                           24
# O_Mouse_38                           16
# O_Mouse_39                           11
# O_Mouse_40                            9
# O_Mouse_43                           33
# O_Mouse_44                           14
# Y_Mouse_19                           12
# Y_Mouse_20                           10
# Y_Mouse_21                            5
# Y_Mouse_22                            7
# Y_Mouse_27                            8
# Y_Mouse_28                           16
# Y_Mouse_6                             8
# Y_Mouse_7                            10

# bring RNA as main assay fro processing
DefaultAssay(ximerakis.seurat) <- "RNA"

# convert to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
ximerakis.seurat.sce <- as.SingleCellExperiment(ximerakis.seurat)
save(ximerakis.seurat.sce, file = paste(Sys.Date(),"SCP263_Ximerakis_SingleCellExperimnent_object.RData",sep = "_"))
###############################################################################################


###############################################################################################
# 1. Run muscat for pseudobulking and extraction of samples

######################################################
####### Data preparation   ++++   Hajdarovic   #######
######################################################

ximerakis.seurat.sce.cl <- prepSCE(ximerakis.seurat.sce, 
                            kid    = "cell_type_simple" ,  # population assignments
                            gid    = "Age_Group"        ,  # group IDs (ctrl/stim)
                            sid    = "Bio_Group"        ,  # sample IDs (ctrl/stim.1234)
                            drop   = TRUE           )  # drop all other colData columns

# store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk  <- length(kids <- levels(ximerakis.seurat.sce.cl$cluster_id))
ns  <- length(sids <- levels(ximerakis.seurat.sce.cl$sample_id))
names(kids) <- kids; names(sids) <- sids

# Aggregation of single-cell to pseudobulk data
pb <- aggregateData(ximerakis.seurat.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# one list item per cell type
assayNames(pb)
# [1] "Arachnoid_barrier_cells"         "Astrocyte_restricted_precursors" "Astrocytes"                      "Choroid_plexus_epithelial_cells"
# [5] "Dendritic_cells"                 "Endothelial_cells"               "Ependymocytes"                   "Hypendymal_cells"               
# [9] "Microglia"                       "Monocytes"                       "Neurons"                         "Neutrophils"                    
# [13] "NSPCs"                           "Olfactory_ensheathing_glia"      "Oligodendrocytes"                "OPCs"                           
# [17] "Pericytes"                       "Tanycytes"                       "Vascular_cells"                  "Vascular_smooth_muscle_cells" 
# Number of cells in each sample and cell type

cell.per.samp.tab <- t(table(ximerakis.seurat.sce.cl$cluster_id, ximerakis.seurat.sce.cl$sample_id))

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
# [1] "Astrocytes"                 "Endothelial_cells"          "Microglia"                  "Neurons"                    "Olfactory_ensheathing_glia"
# [6] "Oligodendrocytes"           "OPCs"                       "Pericytes"      

# extract pseudobulk information for samples that pass the cell number cutoff
counts.pb <- counts.pb.tmp[celltype.qc]

#### save counts
save(counts.pb, file = paste0(Sys.Date(),"_SCP263_Ximerakis_muscat_PB_objects_QC_Clean.RData"))
#############################################################################################

##############################################################################################
# 2. Use SVA to clean up batch effects on expression and DEseq2 for DE analysis

# clean up memory and reload only muscat PBs
load('2024-04-11_SCP263_Ximerakis_muscat_PB_objects_QC_Clean.RData')

# will run SVA to clean up noise

######################################################
#######   DEG analysis   ++++   Hajdarovic     #######
######################################################

# import metadata and order it

colnames(counts.pb[[1]])
#  [1] "O_Mouse_33" "O_Mouse_34" "O_Mouse_37" "O_Mouse_38" "O_Mouse_39" "O_Mouse_40" "O_Mouse_43" "O_Mouse_44" "Y_Mouse_19" "Y_Mouse_20" "Y_Mouse_21" "Y_Mouse_22"
# [13] "Y_Mouse_27" "Y_Mouse_28" "Y_Mouse_6"  "Y_Mouse_7" 

my.meta           <- data.frame(matrix(0,16,4))
colnames(my.meta) <- c("SampleID","Age_Group", "Age_month", "Sex")
my.meta$Age_Group <- factor(c(rep("Y",8),rep("O",8)) , levels = c("Y","O"))
my.meta$Age_month <- c(rep(2.5,8),rep(21.5,8))
my.meta$Sex       <- rep("M",16)
my.meta$SampleID  <- rev(sort(colnames(counts.pb[[1]])))
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
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_Pseudobulk_SCP263_Ximerakis_",names(counts.pb)[[i]])
  
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
  my.cols  <- rep("",nrow(my.meta))
  my.cols[my.meta$Age_Group %in% "Y"] <- "deepskyblue"
  my.cols[my.meta$Age_Group %in% "O"] <- "deepskyblue4"

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
my.rdata.age <- paste0(Sys.Date(),"_PB_SCP263_Ximerakis_DEseq2_objects.RData")
save(deseq.res.list, deseq.res.list, file = my.rdata.age)

my.vst.age <- paste0(Sys.Date(),"_PB_SCP263_Ximerakis_VST_data_objects.RData")
save(vst.cts, file = my.vst.age)
##############################################################################################

#######################
sink(file = paste0(Sys.Date(),"_SCP263_Ximerakis_processing_session_Info.txt"))
sessionInfo()
sink()
