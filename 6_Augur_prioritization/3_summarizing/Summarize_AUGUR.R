setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/Integration_snRNA_snATAC/AUGUR')
options(stringsAsFactors = F)
options (future.globals.maxSize = 32000 * 1024^2)

#### Packages
library(Augur)
library(viridis)

library(ggplot2)          # 
library(scales)           # 
library(ComplexHeatmap)   #
library(circlize)         #

library(datawizard)       # rank transformation

theme_set(theme_bw()) 

# 2024-04-22
# run AUGUR
# run vanilla R due to memory constraints

# 2025-06-25
# Summarize AUGUR results across modalities

###########################################################################################
# 1. Load AUGUR results

load('../../scATACseq/snATAC_Brain_Aging_Meta/Downstream_Analyses/AUGUR/2025-06-25_Augur_Killi_brain_ATAC_2Cohorts_GeneActivityScore_AGING_ONLY.RData')
augur.brain.atac.f
augur.brain.atac.m

load('../../scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/AUGUR/2024-04-19_Augur_Killi_Brain_scRNAseq_GRZ_by_Sex.RData')
augur.brain.grz.f
augur.brain.grz.m

load('../../scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/AUGUR/2024-04-19_Augur_Killi_Brain_scRNAseq_ZMZ_by_Sex.RData')
augur.brain.zmz.f
augur.brain.zmz.m
###################################################################################################################

###########################################################################################
# 2. Merge AUGUR results across modalities

## stepwise merge
my.auc.tmp.1 <- merge(augur.brain.atac.f$AUC, augur.brain.atac.m$AUC, by = "cell_type", suffixes = c(".atac.f",".atac.m"))
my.auc.tmp.2 <- merge(augur.brain.grz.f$AUC , augur.brain.grz.m$AUC , by = "cell_type", suffixes = c(".rna.grz.f",".rna.grz.m"))
my.auc.tmp.3 <- merge(augur.brain.zmz.f$AUC , augur.brain.zmz.m$AUC , by = "cell_type", suffixes = c(".rna.zmz.f",".rna.zmz.m"))
my.auc.tmp.4 <- merge(my.auc.tmp.2, my.auc.tmp.3, by = "cell_type")
my.auc.merge <- merge(my.auc.tmp.4, my.auc.tmp.1, by = "cell_type")
rownames(my.auc.merge) <- my.auc.merge$cell_type


auc.ranks <- apply(my.auc.merge[,-1], 2, ranktransform)
auc.ranks.srt <- auc.ranks[sort(apply(auc.ranks,1,prod), decreasing = T, index.return = T)$ix,]

pdf(paste0(Sys.Date(),"_AUGUR_AUC_Ranks_GRZ_ZMZ_snRNA_GRZ_snATAC.pdf"), width = 10, height = 5)
Heatmap(as.matrix(auc.ranks.srt), 
        col = colorRamp2(c(0, 6,  12), c("darkgreen", "white", "darkorchid"), transparency = 0, space = "LAB"),
        heatmap_legend_param = list(title = "Rank", at = c(0, 6, 12), labels = c("Low", "middle", "High")),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
        cluster_rows = F,
        cluster_columns = F,
        column_title = "AUGUR AUCs ranks",
        width  = 6*unit(7, "mm"),
        height = nrow(my.auc.merge)*unit(7, "mm") )
dev.off()

###########################################################################################

#######################
sink(file = paste(Sys.Date(),"_Killi_brain_ATAC_by_sex_AUGUR_session_Info.txt", sep =""))
sessionInfo()
sink()
