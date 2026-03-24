setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Atlas_Analysis')
options(stringsAsFactors = F)
options (future.globals.maxSize = 32000 * 1024^2)

# Load packages
library('Seurat')         # 
library(sctransform)      # 
library("singleCellTK")   # 
library("scater")         # 


#####################################################################################################################
# 2024-03-11
# Plot variance explained
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

killi.brain.clean.sce <- as.SingleCellExperiment(killi.brain.clean)

# Computing variance explained on the log-counts
# https://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html
vars <- getVarianceExplained(killi.brain.clean.sce,
                             variables=c("Sex", "Age_weeks", "Strain", "Phase", "Cell_Identity"))
head(vars)
#                       Sex  Age_weeks       Strain       Phase Cell_Identity
# LOC107382895 0.0006632060 0.03796509 0.0018246477 0.004513751    0.05509551
# LOC107382813 0.0020425507 0.06207491 0.0101245555 0.005050294    0.41687373
# LOC107382745 0.0007047965 0.01142275 0.0005713469 0.003909957    0.02389697
# fmr1         0.0007845892 0.08983037 0.0009244118 0.008738652    0.06527242
# aff2         0.0076241724 0.37342631 0.1140107019 0.085174348    5.99048551
# ids          0.0001684893 0.05665763 0.0122698040 0.003662576    0.61604589


pdf(paste0(Sys.Date(),"_Variance_Explained_plot_Killi_Brain_Aging_Dataset.pdf"), width = 4, height = 4)
plotExplanatoryVariables(vars[,c("Sex", "Age_weeks", "Strain", "Cell_Identity")])
dev.off()
###############################################################################################################################################################
sink(file = paste(Sys.Date(),"_Killi_Brain_Atlas_MarkerPlotting_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()
