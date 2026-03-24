setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/AUGUR/')
options(stringsAsFactors=FALSE)

# scRNA tools
library(Seurat)
library(Augur)
library(viridis)

# 2024-04-19
# run AUGUR
# run vanilla R due to memory constraints


## ########################################################################################################################
## ## Import/split annotated Seurat object
## 
## # Import final annotation
## load('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Atlas_Analysis/2023-08-23_Seurat_object_with_Manual_annotation_FINAL.RData')
## killi.brain.clean
## # An object of class Seurat 
## # 21160 features across 209939 samples within 1 assay 
## # Active assay: RNA (21160 features, 5000 variable features)
## # 3 dimensional reductions calculated: pca, umap, harmony
## 
## # subset per strains
## killi.brain.grz <- subset(killi.brain.clean, subset = Strain %in% "GRZ")    # 104,665 cells
## killi.brain.zmz <- subset(killi.brain.clean, subset = Strain %in% "ZMZ")    # 105,274 cells
## save(killi.brain.grz, killi.brain.zmz, 
##      file = paste(Sys.Date(),"Seurat_objects_SPLIT_PER_STRAIN_with_Manual_annotation_FINAL.RData",sep = "_"))
## 
## # and by sex
## killi.brain.grz.f <- subset(killi.brain.grz, subset = Sex %in% "F")    # 52,003 cells
## killi.brain.grz.m <- subset(killi.brain.grz, subset = Sex %in% "M")    # 52,662 cells
## killi.brain.zmz.f <- subset(killi.brain.zmz, subset = Sex %in% "F")    # 51,399 cells
## killi.brain.zmz.m <- subset(killi.brain.zmz, subset = Sex %in% "M")    # 53,875 cells
## save(killi.brain.grz.f, killi.brain.grz.m,
##      killi.brain.zmz.f, killi.brain.zmz.m, 
##      file = paste(Sys.Date(),"Seurat_objects_SPLIT_PER_STRAIN_and_SEX_with_Manual_annotation_FINAL.RData",sep = "_"))
## 
## # rm(killi.brain.clean, killi.brain.grz, killi.brain.zmz)
## ########################################################################################################################

########################################################################################################################
# Run AUGUR based on transferred cell type labels

load('2024-02-16_Seurat_objects_SPLIT_PER_STRAIN_and_SEX_with_Manual_annotation_FINAL.RData')

####################################
############     GRZ    ############
####################################

#augur.brain.grz.f <-  calculate_auc(killi.brain.grz.f,
#                                  cell_type_col = "Cell_Identity", 
#                                  label_col = "Group",
#                                  n_threads = 1,
#                                  min_cells = 100)
#
#augur.brain.grz.m <-  calculate_auc(killi.brain.grz.m,
#                                  cell_type_col = "Cell_Identity", 
#                                  label_col = "Group",
#                                  n_threads = 1,
#                                  min_cells = 100)
#
#
#save(augur.brain.grz.f, augur.brain.grz.m, file = paste0(Sys.Date(),"_Augur_Killi_Brain_scRNAseq_GRZ_by_Sex.RData"))

load('2024-04-19_Augur_Killi_Brain_scRNAseq_GRZ_by_Sex.RData')

pdf(paste0(Sys.Date(),"_Augur_Killi_Brain_Aging_scRNAseq_GRZ_Female.pdf"), width = 3, height = 3)
plot_umap(augur.brain.grz.f, killi.brain.grz.f, cell_type_col = "Cell_Identity", palette = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","lightgrey","#CCCCFF","#9999FF","#333399")))(50))
dev.off()

pdf(paste0(Sys.Date(),"_Augur_Killi_Brain_Aging_scRNAseq_GRZ_Male.pdf"), width = 3, height = 3)
plot_umap(augur.brain.grz.m, killi.brain.grz.m, cell_type_col = "Cell_Identity", palette = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","lightgrey","#CCCCFF","#9999FF","#333399")))(50))
dev.off()

pdf(paste0(Sys.Date(),"_Augur_Killi_Brain_Aging_scRNAseq_GRZ_Female_Lollipop.pdf"), width = 3, height = 3)
plot_lollipop(augur.brain.grz.f )
dev.off()

pdf(paste0(Sys.Date(),"_Augur_Killi_Brain_Aging_scRNAseq_GRZ_Male_Lollipop.pdf"), width = 3, height = 3)
plot_lollipop(augur.brain.grz.m)
dev.off()

####################################
############     ZMZ    ############
####################################

augur.brain.zmz.f <-  calculate_auc(killi.brain.zmz.f,
                                  cell_type_col = "Cell_Identity", 
                                  label_col = "Group",
                                  n_threads = 1,
                                  min_cells = 100)

augur.brain.zmz.m <-  calculate_auc(killi.brain.zmz.m,
                                  cell_type_col = "Cell_Identity", 
                                  label_col = "Group",
                                  n_threads = 1,
                                  min_cells = 100)


save(augur.brain.zmz.f, augur.brain.zmz.m, file = paste0(Sys.Date(),"_Augur_Killi_Brain_scRNAseq_ZMZ_by_Sex.RData"))

pdf(paste0(Sys.Date(),"_Augur_Killi_Brain_Aging_scRNAseq_ZMZ_Female.pdf"), width = 3, height = 3)
plot_umap(augur.brain.zmz.f, killi.brain.zmz.f, cell_type_col = "Cell_Identity", palette = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","lightgrey","#CCCCFF","#9999FF","#333399")))(50))
dev.off()

pdf(paste0(Sys.Date(),"_Augur_Killi_Brain_Aging_scRNAseq_ZMZ_Male.pdf"), width = 3, height = 3)
plot_umap(augur.brain.zmz.m, killi.brain.zmz.m, cell_type_col = "Cell_Identity", palette = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","lightgrey","#CCCCFF","#9999FF","#333399")))(50))
dev.off()

pdf(paste0(Sys.Date(),"_Augur_Killi_Brain_Aging_scRNAseq_ZMZ_Female_Lollipop.pdf"), width = 3, height = 3)
plot_lollipop(augur.brain.zmz.f)
dev.off()

pdf(paste0(Sys.Date(),"_Augur_Killi_Brain_Aging_scRNAseq_ZMZ_Male_Lollipop.pdf"), width = 3, height = 3)
plot_lollipop(augur.brain.zmz.m)
dev.off()

###################################################################################################################


#######################
sink(file = paste(Sys.Date(),"_Killi_brain_Aging_scRNAseq_AUGUR_session_Info.txt", sep =""))
sessionInfo()
sink()
