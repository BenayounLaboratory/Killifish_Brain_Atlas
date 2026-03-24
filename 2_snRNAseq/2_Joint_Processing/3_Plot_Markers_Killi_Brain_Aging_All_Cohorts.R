setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Atlas_Analysis')
options(stringsAsFactors = F)
options (future.globals.maxSize = 32000 * 1024^2)

# Load packages
library('Seurat')    # 
library(sctransform) # 
library(clustree)    # 
library(scales)      # 
library(dplyr)      # 
library(readxl)
library(Polychrome)


#####################################################################################################################
# 2024-03-11
# Plot for manuscript QC
#
#####################################################################################################################

#####################################################################################################################
#### 1. Load Cleaned up annotated Seurat Object

# Import final annotation
load('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Atlas_Analysis/2023-08-23_Seurat_object_with_Manual_annotation_FINAL.RData')
killi.brain.clean
# An object of class Seurat
# 21160 features across 209939 samples within 1 assay
# Active assay: RNA (21160 features, 5000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony


# astrocyte/radial glia
# ependymal cells               ### LOC107386535 (ependymin-2-like), LOC107383970	(serotransferrin-like)
# Erythrocytes (hemoglobin)
# mph/microglia                 ### LOC107387973 (itgam/cd11b), LOC107379395	(apoeb)
# NSPCs                         ### gli3, efna2, LOC107385497 (msi1)
# Oligodendrocyte               ### LOC107394899 (mbpa), LOC107386530 (plp)
# OPCs
# Vascular smooth muscle        ### LOC107375895 (kcnq5)

# GABAergic neurons             ### LOC107384443 (gad1l), LOC107391088 (scl6a1; GABA transporter)
# Granule excitatory neuron     ### LOC107392205 (qka), LOC107376653 (Slc17a7/VGLUT1)
# Neurons_misc_1                ### LOC107392873 (fcho1),LOC107397057 (gjc1),LOC107387921 (ngfr)
# Neurons_misc_2                ### LOC107385736 (grm2), LOC107385201 (camkv)
# Neurons_misc_3                ### LOC107380370 (elnl)
# Neurons_misc_4                ### LOC107375999 (ebf3)
# Purkinje cells                ### LOC107379976 (aldoca)
# PV interneurons
# mature neuron markers

# LOC107378372 hemoglobin subunit beta-A hbba


pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_Marker_UMAP_postHarmony_FINAL_s100b_Astrocytes.pdf", sep = "_"), height = 4, width = 5)
FeaturePlot(killi.brain.clean, features = c("s100b")         , cols = c("lightgrey", "darkred"), pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_Marker_UMAP_postHarmony_FINAL_LOC107386535_ependymin-2-like_ependymal_cells.pdf", sep = "_"), height = 4, width = 5)
FeaturePlot(killi.brain.clean, features = c("LOC107386535")  , cols = c("lightgrey", "darkred"), pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_Marker_UMAP_postHarmony_FINAL_LOC107378372_hbba_Erythrocytes.pdf", sep = "_"), height = 4, width = 5)
FeaturePlot(killi.brain.clean, features = c("LOC107378372")  , cols = c("lightgrey", "darkred"), pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_Marker_UMAP_postHarmony_FINAL_LOC107387973_itgam_microglia.pdf", sep = "_"), height = 4, width = 5)
FeaturePlot(killi.brain.clean, features = c("LOC107387973")  , cols = c("lightgrey", "darkred"), pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_Marker_UMAP_postHarmony_FINAL_LOC107379395_apoeb_microglia.pdf", sep = "_"), height = 4, width = 5)
FeaturePlot(killi.brain.clean, features = c("LOC107379395")  , cols = c("lightgrey", "darkred"), pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_Marker_UMAP_postHarmony_FINAL_gli3_NSPCs.pdf", sep = "_"), height = 4, width = 5)
FeaturePlot(killi.brain.clean, features = c("gli3")         , cols = c("lightgrey", "darkred"), pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_Marker_UMAP_postHarmony_FINAL_mpz_Oligodendrocyte.pdf", sep = "_"), height = 4, width = 5)
FeaturePlot(killi.brain.clean, features = c("mpz")          , cols = c("lightgrey", "darkred"), pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_Marker_UMAP_postHarmony_FINAL_sema5a_OPCs.pdf", sep = "_"), height = 4, width = 5)
FeaturePlot(killi.brain.clean, features = c("sema5a")       , cols = c("lightgrey", "darkred"), pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_Marker_UMAP_postHarmony_FINAL_tpm1_VascularSmoothMuscle.pdf", sep = "_"), height = 4, width = 5)
FeaturePlot(killi.brain.clean, features = c("tpm1")        , cols = c("lightgrey", "darkred"), pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_Marker_UMAP_postHarmony_FINAL_LOC107384443_gad1l_GABAergic.pdf", sep = "_"), height = 4, width = 5)
FeaturePlot(killi.brain.clean, features = c("LOC107384443")  , cols = c("lightgrey", "darkred"), pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_Marker_UMAP_postHarmony_FINAL_fat2_GranuleExcitatoryNeuron.pdf", sep = "_"), height = 4, width = 5)
FeaturePlot(killi.brain.clean, features = c("fat2")          , cols = c("lightgrey", "darkred"), pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_Marker_UMAP_postHarmony_FINAL_itpr1_PurkinjeCells.pdf", sep = "_"), height = 4, width = 5)
FeaturePlot(killi.brain.clean, features = c("itpr1")          , cols = c("lightgrey", "darkred"), pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_Marker_UMAP_postHarmony_FINAL_pvalb_PVInterneuron.pdf", sep = "_"), height = 4, width = 5)
FeaturePlot(killi.brain.clean, features = c("pvalb")          , cols = c("lightgrey", "darkred"), pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
dev.off()


pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_Marker_UMAP_postHarmony_FINAL_map2_Neuron.pdf", sep = "_"), height = 4, width = 5)
FeaturePlot(killi.brain.clean, features = c("map2")          , cols = c("lightgrey", "darkred"), pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
dev.off()


###############################################################################################################################################################
sink(file = paste(Sys.Date(),"_Killi_Brain_Atlas_MarkerPlotting_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()

