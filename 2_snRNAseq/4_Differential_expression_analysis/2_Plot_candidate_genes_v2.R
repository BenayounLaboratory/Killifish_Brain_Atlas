setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Differential_Expression_Age/DGE_Analysis')
options(stringsAsFactors = F)

#### Packages
library('Seurat')         # 
library(sctransform)      # 
library("singleCellTK")   # 

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

# 2024-04-25
# Plot genes in Nr3c1/LOC107375066 pathway
# https://www.ncbi.nlm.nih.gov/gene/?term=LOC107375066

# 2025-03-18
#  plot top genes selected in RNAscope
# include PB heatmap
##### - LOC107374900	optineurin
##### - helz2	helicase with zinc finger domain 2
##### - stat1	signal transducer and activator of transcription 1
##### - znfx1	zinc finger, NFX1-type containing 1
##### - LOC107387991	irf3- interferon regulatory factor 3
##### - ifih1	interferon induced with helicase C domain 1
##########################################################################################################################################

##########################################################################################################################################
### 1. prepare data
load('2024-02-16_Seurat_objects_SPLIT_PER_STRAIN_withSampleID_with_Manual_annotation_FINAL.RData')
load('2024-02-28_pseudobulk_killi_cell_types_AGING_GRZ_DEseq2_objects.RData')
load('2024-02-28_pseudobulk_killi_cell_types_AGING_ZMZ_DEseq2_objects.RData')
load('2024-02-28_pseudobulk_killi_cell_types_AGING_GRZ_VST_data_objects.RData')
load('2024-02-28_pseudobulk_killi_cell_types_AGING_ZMZ_VST_data_objects.RData')
qc.celltypes <- names(deseq.res.list.genes.grz)


# subset to QC cell types
killi.brain.grz <- subset(killi.brain.grz, subset = Cell_Identity %in% qc.celltypes)
killi.brain.zmz <- subset(killi.brain.zmz, subset = Cell_Identity %in% qc.celltypes)

killi.brain.grz$Group <- factor(killi.brain.grz$Group, levels = c("GRZ_Y_F","GRZ_Y_M","GRZ_M_F","GRZ_M_M","GRZ_O_F","GRZ_O_M"))
killi.brain.zmz$Group <- factor(killi.brain.zmz$Group, levels = c("ZMZ_Y_F","ZMZ_Y_M",
                                                                          "ZMZ_M_F","ZMZ_M_M",
                                                                          "ZMZ_O_F","ZMZ_O_M",
                                                                          "ZMZ_G_F","ZMZ_G_M"))

killi.brain.grz$Age_Group <- factor(killi.brain.grz$Age_Group, levels = c("Y","M","O"))
killi.brain.zmz$Age_Group <- factor(killi.brain.zmz$Age_Group, levels = c("Y","M","O", "G"))

killi.brain.grz$Age_weeks <- as.numeric(killi.brain.grz$Age_weeks)
killi.brain.zmz$Age_weeks <- as.numeric(killi.brain.zmz$Age_weeks)


##########################################################################################################################################
### 2. plot expression


pdf(paste0(Sys.Date(),"_Seurat_violinPlot_expression_Nr3c1_perBioGroup_GRZ.pdf"), height = 4, width = 8)
VlnPlot(killi.brain.grz, 
        features = "LOC107375066", 
        split.by = "Group", 
        cols = c("deeppink" ,"deepskyblue",
                 "deeppink3","deepskyblue3",
                 "deeppink4","deepskyblue4"),
        pt.size = 0.5, assay = 'RNA', flip = T)
dev.off()  

pdf(paste0(Sys.Date(),"_Seurat_violinPlot_expression_Nr3c1_perBioGroup_ZMZ.pdf"), height = 4, width = 9)
VlnPlot(killi.brain.zmz, 
        features = "LOC107375066", 
        split.by = "Group", 
        cols = c("deeppink" ,"deepskyblue",
                 "deeppink3","deepskyblue3",
                 "deeppink4","deepskyblue4",
                 "magenta4","royalblue4"),
        pt.size = 0.5, assay = 'RNA', flip = T)
dev.off()  

##########################################################################################################################################

#######################
sink(file = paste(Sys.Date(),"_MuscatDEseq2_PB_DESeq2_GSEA_scRNAseq_BrainAtlas_session_Info.txt", sep =""))
sessionInfo()
sink()