setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/UCell_Geneset')
options(stringsAsFactors = F)

# Packages
library(Seurat)            #
library(ggplot2)           #
library(harmony)           #
library(Hmisc)             #

library(ComplexHeatmap)    #
library(circlize)          #
library(viridis)           #
library(RColorBrewer)      #
library(gridExtra)

library(readxl)
library(bitops)

library(UCell)

# 2024-03-08
# Run Ucell on desired relevant gene lists for important brain cell types

# 2024-03-20
# Updated gene sets

# 2024-04-12
# Clean/simplify output

####################################################################################################################
# 0. Load annotated Seurat Objects 

# # Import final annotation
# load('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Atlas_Analysis/2023-08-23_Seurat_object_with_Manual_annotation_FINAL.RData')
# killi.brain.clean
# # An object of class Seurat 
# # 21160 features across 209939 samples within 1 assay 
# # Active assay: RNA (21160 features, 5000 variable features)
# # 3 dimensional reductions calculated: pca, umap, harmony
# 
# # Add in a sample ID
# killi.brain.clean@meta.data$SampleID <- paste0(killi.brain.clean@meta.data$Group,"_",killi.brain.clean@meta.data$Batch)
# 
# # subset per strains
# killi.brain.grz <- subset(killi.brain.clean, subset = Strain %in% "GRZ")    # 104,665 cells
# killi.brain.zmz <- subset(killi.brain.clean, subset = Strain %in% "ZMZ")    # 105,274 cells
# save(killi.brain.grz, killi.brain.zmz, 
#      file = paste(Sys.Date(),"Seurat_objects_SPLIT_PER_STRAIN_withSampleID_with_Manual_annotation_FINAL.RData",sep = "_"))
#
# load("../Differential_Expression/2024-02-16_Seurat_objects_SPLIT_PER_STRAIN_withSampleID_with_Manual_annotation_FINAL.RData")
# 
# # subset QC cell types (from DEG analysis)
# both.celltype.qc <- c("Astrocytes_Radial_Glia","GABAergic_neurons","Granule_Excitatory_Neurons", "Microglia","NSPCs","Oligodendrocytes","OPCs","PV_interneurons"  )
# 
# # subset per strains
# killi.brain.grz <- subset(killi.brain.grz, subset = Cell_Identity %in% both.celltype.qc)    # 46,222 cells
# killi.brain.zmz <- subset(killi.brain.zmz, subset = Cell_Identity %in% both.celltype.qc)    # 48,123 cells
# save(killi.brain.grz, killi.brain.zmz, 
#      file = paste(Sys.Date(),"Seurat_objects_SPLIT_PER_STRAIN_Only_QC_CellTypes.RData",sep = "_"))
####################################################################################################################

####################################################################################################################
# 1. Run UCell on gene lists of interest

# load single cell
load('2024-03-08_Seurat_objects_SPLIT_PER_STRAIN_Only_QC_CellTypes.RData')

# load gene sets
load('../Pathways_GeneSets/2024-04-12_Aging_Gene_lists_killified.RData')

##################################################
####################   GRZ   #####################
##################################################

killi.brain.grz <- AddModuleScore_UCell(killi.brain.grz, features=aging.glists.killi, name=NULL)
killi.brain.grz$Age_Group <- factor(killi.brain.grz$Age_Group, levels = c("Y", "M", "O"))

# split seurat object by cell type
grz.list <- SplitObject(killi.brain.grz, split.by = "Cell_Identity")

my.sort.names <- c("SenMayo", "Hahn_CAG", "GTEx_UP","GTEx_DWN" )

# create pval matrix
grz.stats           <- data.frame(matrix(NA,8,5))
colnames(grz.stats) <- c("Cell_Type","Wilcoxon_P_SenMayo","Wilcoxon_P_Hahn_CAG", "Wilcoxon_P_GTEx_UP", "Wilcoxon_P_GTEx_DWN")
grz.stats$Cell_Type <- names(grz.list)

my.plot.list        <- vector(mode = "list", length = length(grz.list))
names(my.plot.list) <- names(grz.list)

for (i in 1:length(grz.list)) {
  
  grz.dot <- DotPlot(grz.list[[i]], features = my.sort.names, group.by = "Age_Group", col.min = -0.5,col.max = 0.5) 
  grz.dot <- grz.dot + scale_colour_gradient2(low = "#333399", mid = "lightgrey", high = "#CC3333")
  grz.dot <- grz.dot + coord_flip() + scale_size(range = c(3, 10)) + scale_size_area(limits=c(30,100))
  grz.dot <- grz.dot + ggtitle(names(grz.list)[i]) 
  # grz.dot
  
  my.plot.list[[i]] <- grz.dot
  
  pdf(paste0(Sys.Date(),"_DotPLot_Aging_GeneLists_UCell_Scores_GRZ_",names(grz.list)[i],".pdf"), width = 5, height = 4)
  print(grz.dot)
  dev.off()
  
  hahn.test     <- wilcox.test(Hahn_CAG   ~  Age_Group, data = grz.list[[i]]@meta.data[grz.list[[i]]$Age_Group %in% c("Y", "O"),]) 
  SenMayo.test  <- wilcox.test(SenMayo    ~  Age_Group, data = grz.list[[i]]@meta.data[grz.list[[i]]$Age_Group %in% c("Y", "O"),]) 
  gtex.up.test  <- wilcox.test(GTEx_UP    ~  Age_Group, data = grz.list[[i]]@meta.data[grz.list[[i]]$Age_Group %in% c("Y", "O"),]) 
  gtex.dwn.test <- wilcox.test(GTEx_DWN   ~  Age_Group, data = grz.list[[i]]@meta.data[grz.list[[i]]$Age_Group %in% c("Y", "O"),]) 

  grz.stats[i,]$Wilcoxon_P_Hahn_CAG     <- hahn.test$p.value    
  grz.stats[i,]$Wilcoxon_P_SenMayo     <- SenMayo.test$p.value  
  grz.stats[i,]$Wilcoxon_P_GTEx_UP     <- gtex.up.test$p.value  
  grz.stats[i,]$Wilcoxon_P_GTEx_DWN    <- gtex.dwn.test$p.value  

}

write.table(t(grz.stats), file = paste0(Sys.Date(),"_Enrichment_Wilcoxon_Aging_GeneLists_UCell_Scores_GRZ.txt"), sep = "\t", col.names = F, quote = F)

pdf(paste0(Sys.Date(),"_DotPLot_Aging_GeneLists_UCell_Scores_GRZ_ALL_files.pdf"), width = 35, height = 4 )
gridExtra::grid.arrange(grobs = my.plot.list[sort(names(grz.list))], nrow = 1)
dev.off()


##################################################
####################   ZMZ   #####################
##################################################

killi.brain.zmz <- AddModuleScore_UCell(killi.brain.zmz, features=aging.glists.killi, name=NULL)
killi.brain.zmz$Age_Group <- factor(killi.brain.zmz$Age_Group, levels = c("Y", "M", "O", "G"))

# split seurat object by cell type
zmz.list <- SplitObject(killi.brain.zmz, split.by = "Cell_Identity")

# create pval matrix
zmz.stats           <- data.frame(matrix(0,8,5))
colnames(zmz.stats) <- c("Cell_Type","Wilcoxon_P_SenMayo","Wilcoxon_P_Hahn_CAG", "Wilcoxon_P_GTEx_UP", "Wilcoxon_P_GTEx_DWN")
zmz.stats$Cell_Type <- names(zmz.list)

my.plot.list        <- vector(mode = "list", length = length(grz.list))
names(my.plot.list) <- names(zmz.list)

for (i in 1:length(zmz.list)) {
  
  zmz.dot <- DotPlot(zmz.list[[i]], features = my.sort.names, group.by = "Age_Group", col.min = -0.5,col.max = 0.5) 
  zmz.dot <- zmz.dot + scale_colour_gradient2(low = "#333399", mid = "lightgrey", high = "#CC3333")
  zmz.dot <- zmz.dot + coord_flip() + scale_size(range = c(3, 10)) + scale_size_area(limits=c(30,100))
  zmz.dot <- zmz.dot + ggtitle(names(zmz.list)[i]) 
  # zmz.dot
  
  my.plot.list[[i]] <- zmz.dot
  
  pdf(paste0(Sys.Date(),"_DotPLot_Aging_GeneLists_UCell_Scores_ZMZ_",names(zmz.list)[i],".pdf"), width = 5, height = 4)
  print(zmz.dot)
  dev.off()
  
  hahn.test     <- wilcox.test(Hahn_CAG   ~  Age_Group, data = zmz.list[[i]]@meta.data[zmz.list[[i]]$Age_Group %in% c("Y", "G"),]) 
  SenMayo.test  <- wilcox.test(SenMayo    ~  Age_Group, data = zmz.list[[i]]@meta.data[zmz.list[[i]]$Age_Group %in% c("Y", "G"),]) 
  gtex.up.test  <- wilcox.test(GTEx_UP    ~  Age_Group, data = zmz.list[[i]]@meta.data[zmz.list[[i]]$Age_Group %in% c("Y", "G"),]) 
  gtex.dwn.test <- wilcox.test(GTEx_DWN   ~  Age_Group, data = zmz.list[[i]]@meta.data[zmz.list[[i]]$Age_Group %in% c("Y", "G"),]) 

  zmz.stats[i,]$Wilcoxon_P_Hahn_CAG     <- hahn.test$p.value    
  zmz.stats[i,]$Wilcoxon_P_SenMayo     <- SenMayo.test$p.value  
  zmz.stats[i,]$Wilcoxon_P_GTEx_UP     <- gtex.up.test$p.value  
  zmz.stats[i,]$Wilcoxon_P_GTEx_DWN    <- gtex.dwn.test$p.value  

}


write.table(t(zmz.stats), file = paste0(Sys.Date(),"_Enrichment_Wilcoxon_Aging_GeneLists_UCell_Scores_ZMZ.txt"), sep = "\t", col.names = F, quote = F)

pdf(paste0(Sys.Date(),"_DotPLot_Aging_GeneLists_UCell_Scores_ZMZ_ALL_files.pdf"), width = 35, height = 4 )
gridExtra::grid.arrange(grobs = my.plot.list[sort(names(zmz.list))], nrow = 1)
dev.off()

#############################################################


#############################################################
sink(file = paste(Sys.Date(),"_Ucell_AgingLists_Scoring_session_Info.txt", sep =""))
sessionInfo()
sink()