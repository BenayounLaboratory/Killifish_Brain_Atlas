setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Mouse_Datasets_for_Comparison/Processed_Objects/Hahn_GSE212576')
options(stringsAsFactors = F)

library('Seurat')
library(singleCellNet)

# 2024-03-26
# annotate hahn data using Hajdarovic object and SingleCellNet

#############################################################################################################
#### 1. Load hahn and Hajdarovic objects

# filtered singlets from hahn
load('2024-03-26_Hahn_Brain_Seurat_object_SINGLETS_ONLY.RData')
hahn.singlets
# An object of class Seurat 
# 58498 features across 75146 samples within 2 assays 
# Active assay: SCT (26213 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

# Hajdarovic object
hypo.seurat <- readRDS('../Hajdarovic_GSE188646_hypo.integrated.final.20210719.RDS')
hypo.seurat
# An object of class Seurat 
# 27135 features across 40064 samples within 2 assays 
# Active assay: integrated (2000 features, 2000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

table(hypo.seurat@meta.data$group)
# Oligodendrocyte                           Neuron                        Astrocyte                              OPC 
# 5027                            26940                             4494                             1093 
# Microglia/Macrophage                         Tanycyte                     Ependymocyte        Pericyte/Endothelial Cell 
# 1218                              688                              202                              244 
# Vascular and Leptomeningeal Cell 
# 158 

# clean cell type names
hypo.seurat@meta.data$CellType <- gsub(" ", "", gsub("/","_", hypo.seurat@meta.data$group))

DefaultAssay(hypo.seurat) <- "RNA"
DefaultAssay(hahn.singlets) <- "RNA"
#############################################################################################################

#############################################################################################################
#### 2. Train SCN model

########################## TRAINING ##############################
# Load query data

genes.hahn <- rownames(hahn.singlets)

# extract info from Seurat object for query data
# exp_type options can be: counts, normcounts, and logcounts, if they are available in your sce object
hypo.seuratfile    <- extractSeurat(hypo.seurat, exp_slot_name = "counts")
hypo.sampTab       <- hypo.seuratfile$sampTab
hypo.expDat        <- hypo.seuratfile$expDat
hypo.sampTab       <- droplevels(hypo.sampTab)
hypo.sampTab$cell  <- rownames(hypo.sampTab)

# Find genes in common to the data sets and limit analysis to these
commonGenes <- intersect(rownames(hypo.expDat), genes.hahn)
length(commonGenes)
# [1] 24573

hypo.expDat     <- hypo.expDat[commonGenes,]

# Split for training and assessment, and transform training data
set.seed(123456789)
stList   <- splitCommon(sampTab=hypo.sampTab, ncells=100, dLevel="CellType")
stTrain  <- stList[[1]]
expTrain <- hypo.expDat[,rownames(stTrain)]

# Train the classifier using the hypo data
class_info <- scn_train(stTrain = stTrain, expTrain = expTrain, 
                        nTopGenes = 100, nTopGenePairs = 50, nRand = 50, nTrees = 1000, 
                        dLevel = "CellType", colName_samp = "cell")
# There are 450 top gene pairs

########################## TESTING ##############################
# Assessing the classifier with heldout data Apply to held out data
stTestList <- splitCommon(sampTab=stList[[2]], ncells=100, dLevel="CellType") # normalize validation data so that the assessment is as fair as possible
stTest     <- stTestList[[1]]
expTest    <- hypo.expDat[commonGenes,rownames(stTest)]

# predict on held out data
classRes_val_all <- scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 100)

# Assess classifier performance on held out data
tm_heldoutassessment <- assess_comm(ct_scores  = classRes_val_all, 
                                    stTrain    = stTrain, 
                                    stQuery    = stTest, 
                                    dLevelSID  = "cell", 
                                    classTrain = "CellType", 
                                    classQuery = "CellType", 
                                    nRand = 100)

pdf(paste0(Sys.Date(),"_PR_curves_SingleCellNet_hypo_model_performance_hahn.pdf") )
plot_PRs(tm_heldoutassessment)
dev.off()

pdf(paste0(Sys.Date(),"_AUPRC_curves_SingleCellNet_hypo_model_performance_hahn.pdf") )
plot_metrics(tm_heldoutassessment)
dev.off()

# Classification result heatmap
# Create a name vector label used later in classification heatmap where the values are cell types/ clusters and names are the sample names
nrand          = 100
sla            = as.vector(stTest$CellType)
names(sla)     = as.vector(stTest$cell)
slaRand        = rep("rand", nrand)
names(slaRand) = paste("rand_", 1:nrand, sep='')
sla            = append(sla, slaRand) # include in the random cells profile created

# Attribution plot
pdf(paste0(Sys.Date(),"_Classification_result_barplot_SingleCellNet_model_performance_hahn.pdf"), height = 5, width = 10 )
plot_attr(classRes=classRes_val_all, sampTab=stTest, nrand=nrand, dLevel="CellType", sid="cell")
dev.off()
# pretty good

########################## QUERY ##############################
# extract info from Seurat object for query data
# exp_type options can be: counts, normcounts, and logcounts, if they are available in your sce object
hahn.seuratfile <- extractSeurat(hahn.singlets, exp_slot_name = "counts")
hahn.sampTab    <- hahn.seuratfile$sampTab
hahn.expDat     <- hahn.seuratfile$expDat

classRes_hahn <- scn_predict(class_info[['cnProc']], hahn.expDat, nrand=50)

# Classification annotation assignment
# This classifies a cell with  the catgory with the highest classification score or higher than a classification score threshold of your choosing.
# The annotation result can be found in a column named category in the query sample table.
sthahn      <- get_cate(classRes = classRes_hahn, sampTab = hahn.sampTab, dLevel = "Age_Group", sid = "SampleID", nrand = 50)

table(sthahn$category,sthahn$SampleID)
#                                CP_O_1 CP_O_2 CP_Y_1 CP_Y_2 HIP_O_1 HIP_O_2 HIP_Y_1 HIP_Y_2
# Astrocyte                        223     64    343     13    1921     680    1532     368
# Ependymocyte                       0      0      3      0       0       5       0       0
# Microglia_Macrophage              99      7    164    155      53     206     234     348
# Neuron                          6736   8784   6034   7824    3910    3304    5305    9623
# OPC                                1      5      1      0       1       1       0       3
# Oligodendrocyte                 3397   3628   2245   2365    1264    1127    1097    1424
# Pericyte_EndothelialCell          43      0     46     15      10      35     117     113
# Tanycyte                           1      0      4      0      22       9      20      17
# VascularandLeptomeningealCell      0      0      3      0       8      13      85      88


# put back into seurat object data
sum(rownames(hahn.singlets@meta.data) == rownames(sthahn)) # 75146

dim(hahn.singlets@meta.data) #  75146    20

hahn.singlets@meta.data$SingleCellNet_Hajdarovic          <- sthahn$category

# remove cells predicted as rand (1 cell)
hahn.singlets       <- subset(hahn.singlets,
                                subset = SingleCellNet_Hajdarovic %in% c("rand"),
                                invert = T)   # no rand cells, 75146 cells


pdf(paste0(Sys.Date(),"_hahn_SingleCellNet_hypo_UMAP_annotated.pdf"), height = 5, width = 6 )
DimPlot(hahn.singlets, reduction = "umap", group.by = "SingleCellNet_Hajdarovic")
dev.off()

save(hahn.singlets, file = paste(Sys.Date(),"hahn_Hip_PC_Brain_SingleCellNet_Annotated_Seurat_object.RData",sep = "_"))



#######################
sink(file = paste(Sys.Date(),"_hahn_Brain_Data_SCN_Annotation_session_Info.txt", sep =""))
sessionInfo()
sink()

