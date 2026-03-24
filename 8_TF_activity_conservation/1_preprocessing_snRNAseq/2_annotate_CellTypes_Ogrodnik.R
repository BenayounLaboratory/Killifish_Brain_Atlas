setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Mouse_Datasets_for_Comparison/Processed_Objects/Ogrodnik_GSE161340/')
options(stringsAsFactors = F)

library('Seurat')
library(singleCellNet)

# 2024-04-12
# annotate data from Ogrognik et al, 2021

#############################################################################################################
#### 1. Load Ogrodnik and Hajdarovic objects

# filtered singlets from Ogrodnik
load('2024-04-12_Ogrodnik_Brain_Seurat_object_SINGLETS_ONLY.RData')
Ogrodnik.singlets
# An object of class Seurat 
# 52162 features across 16060 samples within 2 assays 
# Active assay: SCT (21109 features, 3000 variable features)
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
DefaultAssay(Ogrodnik.singlets) <- "RNA"
#############################################################################################################

#############################################################################################################
#### 2. Train SCN model

########################## TRAINING ##############################
# Load query data

genes.Ogrodnik <- rownames(Ogrodnik.singlets)

# extract info from Seurat object for query data
# exp_type options can be: counts, normcounts, and logcounts, if they are available in your sce object
hypo.seuratfile    <- extractSeurat(hypo.seurat, exp_slot_name = "counts")
hypo.sampTab       <- hypo.seuratfile$sampTab
hypo.expDat        <- hypo.seuratfile$expDat
hypo.sampTab       <- droplevels(hypo.sampTab)
hypo.sampTab$cell  <- rownames(hypo.sampTab)

# Find genes in common to the data sets and limit analysis to these
commonGenes <- intersect(rownames(hypo.expDat), genes.Ogrodnik)
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

pdf(paste0(Sys.Date(),"_PR_curves_SingleCellNet_hypo_model_performance_Ogrodnik.pdf") )
plot_PRs(tm_heldoutassessment)
dev.off()

pdf(paste0(Sys.Date(),"_AUPRC_curves_SingleCellNet_hypo_model_performance_Ogrodnik.pdf") )
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
pdf(paste0(Sys.Date(),"_Classification_result_barplot_SingleCellNet_model_performance_Ogrodnik.pdf"), height = 5, width = 10 )
plot_attr(classRes=classRes_val_all, sampTab=stTest, nrand=nrand, dLevel="CellType", sid="cell")
dev.off()
# pretty good

########################## QUERY ##############################
# extract info from Seurat object for query data
# exp_type options can be: counts, normcounts, and logcounts, if they are available in your sce object
Ogrodnik.seuratfile <- extractSeurat(Ogrodnik.singlets, exp_slot_name = "counts")
Ogrodnik.sampTab    <- Ogrodnik.seuratfile$sampTab
Ogrodnik.expDat     <- Ogrodnik.seuratfile$expDat

classRes_Ogrodnik <- scn_predict(class_info[['cnProc']], Ogrodnik.expDat, nrand=50)

# Classification annotation assignment
# This classifies a cell with  the catgory with the highest classification score or higher than a classification score threshold of your choosing.
# The annotation result can be found in a column named category in the query sample table.
stOgrodnik      <- get_cate(classRes = classRes_Ogrodnik, sampTab = Ogrodnik.sampTab, dLevel = "Age_Group", sid = "SampleID", nrand = 50)

table(stOgrodnik$category,stOgrodnik$SampleID)
#                               HIP_O_1 HIP_O_2 HIP_Y_1 HIP_Y_2
# Astrocyte                         287     551     322     273
# Ependymocyte                        5      27       7       4
# Microglia_Macrophage              388     433     190     178
# Neuron                           2298    3137    1881    1841
# Oligodendrocyte                  1163    1425     677     521
# OPC                                48      53      23      29
# Pericyte_EndothelialCell           15      18       1      14
# rand                                0       1       0       0
# Tanycyte                           36     102      22      21
# VascularandLeptomeningealCell      21      20      12      16


# put back into seurat object data
sum(rownames(Ogrodnik.singlets@meta.data) == rownames(stOgrodnik)) # 16060

dim(Ogrodnik.singlets@meta.data) # 16060    21

Ogrodnik.singlets@meta.data$SingleCellNet_Hajdarovic          <- stOgrodnik$category

# remove cells predicted as rand (1 cell)
Ogrodnik.singlets       <- subset(Ogrodnik.singlets,
                                subset = SingleCellNet_Hajdarovic %in% c("rand"),
                                invert = T)   # no rand cells, 16059 cells


pdf(paste0(Sys.Date(),"_Ogrodnik_SingleCellNet_hypo_UMAP_annotated.pdf"), height = 5, width = 6 )
DimPlot(Ogrodnik.singlets, reduction = "umap", group.by = "SingleCellNet_Hajdarovic")
dev.off()

save(Ogrodnik.singlets, file = paste(Sys.Date(),"Ogrodnik_Hippocampus_SingleCellNet_Annotated_Seurat_object.RData",sep = "_"))



#######################
sink(file = paste(Sys.Date(),"_Ogrodnik_Brain_Data_SCN_Annotation_session_Info.txt", sep =""))
sessionInfo()
sink()

