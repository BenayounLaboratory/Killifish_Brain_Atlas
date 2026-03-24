setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Mouse_Datasets_for_Comparison/Processed_Objects/Allen_GSE207848')
options(stringsAsFactors = F)

library('Seurat')
library(singleCellNet)

# 2024-03-26
# annotate Allen data using Hajdarovic object and SingleCellNet

#############################################################################################################
#### 1. Load Allen and Hajdarovic objects

# filtered singlets from Allen
load('2024-03-26_Allen_Seurat_object_SINGLETS_ONLY.RData')
allen.singlets
# An object of class Seurat 
# 58305 features across 63559 samples within 2 assays 
# Active assay: SCT (26020 features, 3000 variable features)
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
DefaultAssay(allen.singlets) <- "RNA"
#############################################################################################################

#############################################################################################################
#### 2. Train SCN model

########################## TRAINING ##############################
# Load query data

genes.allen <- rownames(allen.singlets)

# extract info from Seurat object for query data
# exp_type options can be: counts, normcounts, and logcounts, if they are available in your sce object
hypo.seuratfile    <- extractSeurat(hypo.seurat, exp_slot_name = "counts")
hypo.sampTab       <- hypo.seuratfile$sampTab
hypo.expDat        <- hypo.seuratfile$expDat
hypo.sampTab       <- droplevels(hypo.sampTab)
hypo.sampTab$cell  <- rownames(hypo.sampTab)

# Find genes in common to the data sets and limit analysis to these
commonGenes <- intersect(rownames(hypo.expDat), genes.allen)
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

pdf(paste0(Sys.Date(),"_PR_curves_SingleCellNet_hypo_model_performance_allen.pdf") )
plot_PRs(tm_heldoutassessment)
dev.off()

pdf(paste0(Sys.Date(),"_AUPRC_curves_SingleCellNet_hypo_model_performance_allen.pdf") )
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
pdf(paste0(Sys.Date(),"_Classification_result_barplot_SingleCellNet_model_performance_Allen.pdf"), height = 5, width = 10 )
plot_attr(classRes=classRes_val_all, sampTab=stTest, nrand=nrand, dLevel="CellType", sid="cell")
dev.off()
# pretty good

########################## QUERY ##############################
# extract info from Seurat object for query data
# exp_type options can be: counts, normcounts, and logcounts, if they are available in your sce object
allen.seuratfile <- extractSeurat(allen.singlets, exp_slot_name = "counts")
allen.sampTab    <- allen.seuratfile$sampTab
allen.expDat     <- allen.seuratfile$expDat

classRes_allen <- scn_predict(class_info[['cnProc']], allen.expDat, nrand=50)

# Classification annotation assignment
# This classifies a cell with  the catgory with the highest classification score or higher than a classification score threshold of your choosing.
# The annotation result can be found in a column named category in the query sample table.
stallen      <- get_cate(classRes = classRes_allen, sampTab = allen.sampTab, dLevel = "Age_Group", sid = "SampleID", nrand = 50)

table(stallen$category,stallen$SampleID)
#                                   O1   O2   O3   O4   Y1   Y2   Y3   Y4
#   Astrocyte                      597  583  827  928 1048  983  998 1034
#   Ependymocyte                     0    0    0    0    0    1    1    0
#   Microglia_Macrophage           692  676  496  620  538  557  363  428
#   Neuron                        2362 1930 1935 2518 5710 4723 3637 3949
#   OPC                            324   96  359   62   37   52  395  432
#   Oligodendrocyte               3626 3433 2631 2839 1384 1285 1034 1222
#   Pericyte_EndothelialCell        53   47   79   77    1   10   99  101
#   Tanycyte                        11    6    5   14  113  165   83   79
#   VascularandLeptomeningealCell   25   24   45   58    1    0   20   13

# put back into seurat object data
sum(rownames(allen.singlets@meta.data) == rownames(stallen)) # 58474

dim(allen.singlets@meta.data) #  58474    20

allen.singlets@meta.data$SingleCellNet_Hajdarovic          <- stallen$category

# remove cells predicted as rand (1 cell)
allen.singlets       <- subset(allen.singlets,
                                subset = SingleCellNet_Hajdarovic %in% c("rand"),
                                invert = T)   # no rand cells, 63559 cells


pdf(paste0(Sys.Date(),"_Allen_SingleCellNet_hypo_UMAP_annotated.pdf"), height = 5, width = 6 )
DimPlot(allen.singlets, reduction = "umap", group.by = "SingleCellNet_Hajdarovic")
dev.off()

save(allen.singlets, file = paste(Sys.Date(),"Allen_PFC_Brain_SingleCellNet_Annotated_Seurat_object.RData",sep = "_"))



#######################
sink(file = paste(Sys.Date(),"_Allen_Brain_Data_SCN_Annotation_session_Info.txt", sep =""))
sessionInfo()
sink()

