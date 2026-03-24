setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Microglia/Machine_Learning')
options(stringsAsFactors = F)

library(Seurat)          #

library(caret)           #
library(randomForest)    # random forest
library(glmnet)          # Elasticnet
library(MLeval)          #

library(ggplot2)         #
library(scales)          #
library(vioplot)         #
library(matrixStats)     # for column multiplication

library(ComplexHeatmap)  
library(circlize)         #

# 2025-07-01
# ML to predict biological age from microglia across strains
# Microglia is the cell type with clearest changes across 'omic' layers
#
# Regression is too difficult of a training problem with *only* 3 time points
# Try a classification young/old and use probability of old
# similar to the CellBiAge paper from Webb lab
# see Yu et al, Cell Reports, 2023 [CellBiAge]

####################################################################################################################
### 1. feature importance
#### https://topepo.github.io/caret/variable-importance.html

load('2025-07-04_ElasticNet_model_Aging_Microglia.RData')
load('2025-07-04_RF_models_Aging_Microglia.RData')

############ Caculate/extract feature importance values
# Create df to receive var imps

# feature importance on coefficients
my.coeffs.enet           <- as.data.frame(coef(my.enet.fit$finalModel, s=my.enet.fit$finalModel$lambdaOpt))
my.coeffs.enet$abs_coeff <- abs(my.coeffs.enet$s1)

# remove intercept
my.coeffs.enet <- my.coeffs.enet[rownames(my.coeffs.enet) != "(Intercept)",]


my.var.Imps.RF          <- as.data.frame(importance(my.rf.fit$finalModel))[,-c(1:2)]

# combine
my.var.imp.merged <- merge(my.coeffs.enet, my.var.Imps.RF, by = "row.names")

# Save feature importance to files
write.table(my.coeffs.enet   , file = paste0(Sys.Date(),"_Coefficient_Parsing_ENET_perCellType_CARET.txt"), quote = F, row.names = T, sep = "\t")
write.table(my.var.Imps.RF   , file = paste0(Sys.Date(),"_Variable_Importance_Parsing_RF_perCellType_CARET.txt") , quote = F, row.names = T, sep = "\t")
write.table(my.var.imp.merged,file = paste0(Sys.Date(),"_Variable_Importance_Coefficient_Parsing_ENET_RF.txt"), quote = F, row.names = F, sep = "\t")



############ Caculate ranks of VarImp
# Since high importance is better, use rank(-x) to sort from largest to smallest
my.var.Imps.RF$Rank_Gini_RF   <- rank(-my.var.Imps.RF$MeanDecreaseGini)
my.coeffs.enet$Rank_absC_ENET <- rank(-as.numeric(my.coeffs.enet$abs_coeff))

Rank.varImps.merged <- merge(my.var.Imps.RF, my.coeffs.enet, by = "row.names")

Rank.varImps.merged$RankProd      <- rowProds(cbind("enet" = Rank.varImps.merged$Rank_absC_ENET, "rf" =  Rank.varImps.merged$Rank_Gini_RF))
Rank.varImps.merged$RankProd_Rank <- rank(Rank.varImps.merged$RankProd )

# Save combined feature importance and rank to files
write.table(Rank.varImps.merged, file = paste0(Sys.Date(),"_Variable_Importance_Parsing_perCellType_CARET_ENET_RF_summary.txt") , quote = F, row.names = F, sep = "\t")

save(Rank.varImps.merged, file = paste0(Sys.Date(),"_CombinedVarImps.RData"))

# grab top 15 based on rank product in that cell type
my.top.vars <- Rank.varImps.merged[Rank.varImps.merged$RankProd_Rank %in% 1:15,]
my.top.vars <- my.top.vars[order(my.top.vars$RankProd_Rank),]
rownames(my.top.vars) <- my.top.vars$Row.names


### heatmap of ranks
pdf(paste0(Sys.Date(),"_Microglai_ML_top15_features_Importance_heatmap.pdf"), width = 12, height = 5)
Heatmap(as.matrix(my.top.vars[,c("Rank_absC_ENET", "Rank_Gini_RF")]), 
                    col = colorRamp2(c(0, 25, 50, 100), c( "darkslategrey", "darkslategray3","lightcyan2","white" ), transparency = 0, space = "LAB"),
                    border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
                    cluster_rows = F,
                    cluster_columns = F,
                    column_title = paste0("Microglia ", "(GRZ)"),
                    width  = 2*unit(6, "mm"),
                    height = 15*unit(6, "mm"))
dev.off()

pdf(paste0(Sys.Date(),"_Microglai_ML_top15_features_ENET_coefficients.pdf"), width = 4, height = 5)
par(oma = c(0.1, 4, 0.1, 0.1))
barplot(my.top.vars$s1, horiz = T, names.arg = my.top.vars$Row.names, las = 1, xlim = c(-0.75, 0.75), col = "black", main = "Top 15  - enet Coefficients")
box()
dev.off()
####################################################################################################################

####################################################################################################################
### 2. expression of top features with aging

load('../../Differential_Expression_Age/DGE_Analysis/2024-02-28_pseudobulk_killi_cell_types_AGING_GRZ_VST_data_objects.RData')
load('../../Differential_Expression_Age/DGE_Analysis/2024-02-28_pseudobulk_killi_cell_types_AGING_ZMZ_VST_data_objects.RData')

grz.mglia <- vst.cts.grz[["Microglia"]]
zmz.mglia <- vst.cts.zmz[["Microglia"]]

top.features <- my.top.vars$Row.names

grz.av     <- data.frame(row.names = top.features)
grz.av$YF  <- apply(grz.mglia[top.features,c(1:4)],1,median)
grz.av$MF  <- apply(grz.mglia[top.features,c(9:11)],1,median)
grz.av$OF  <- apply(grz.mglia[top.features,c(15:17)],1,median)
grz.av$YM  <- apply(grz.mglia[top.features,c(5:8)],1,median)
grz.av$MM  <- apply(grz.mglia[top.features,c(12:14)],1,median)
grz.av$OM  <- apply(grz.mglia[top.features,c(18:20)],1,median)

zmz.av     <- data.frame(row.names = top.features)
zmz.av$YF  <- apply(zmz.mglia[top.features,c(1:3)],1,median)
zmz.av$MF  <- apply(zmz.mglia[top.features,c(7:9)],1,median)
zmz.av$OF  <- apply(zmz.mglia[top.features,c(13:15)],1,median)
zmz.av$GF  <- apply(zmz.mglia[top.features,c(19:21)],1,median)
zmz.av$YM  <- apply(zmz.mglia[top.features,c(4:6)],1,median)
zmz.av$MM  <- apply(zmz.mglia[top.features,c(10:12)],1,median)
zmz.av$OM  <- apply(zmz.mglia[top.features,c(16:18)],1,median)
zmz.av$GM  <- apply(zmz.mglia[top.features,c(22:24)],1,median)

grz.heat <- Heatmap(t(scale(t(grz.av))), 
                    col = colorRamp2(c(-2,0, 2), c("darkblue", "white", "#CC3333"), transparency = 0, space = "LAB"),
                    border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
                    cluster_rows = F,
                    cluster_columns = F,
                    column_title = paste0("Microglia ", "(GRZ)"),
                    width  = ncol(grz.av)*unit(7, "mm"),
                    height = nrow(grz.av)*unit(7, "mm"))

zmz.heat <- Heatmap(t(scale(t(zmz.av))), 
                    col = colorRamp2(c(-2,0, 2), c("darkblue", "white", "#CC3333"), transparency = 0, space = "LAB"),
                    border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
                    cluster_rows = F,
                    cluster_columns = F,
                    column_title = paste0("Microglia ", "(ZMZ1001)"),
                    width  = ncol(zmz.av)*unit(7, "mm"),
                    height = nrow(zmz.av)*unit(7, "mm"))


pdf(paste0(Sys.Date(),"_MICROGLIAset_aging_GRZ_ZMZ_top15features_median_heatmap.pdf"), width = 12, height = 5)
print(grz.heat + zmz.heat)
dev.off()


#######################
sink(file = paste(Sys.Date(),"Machine_Learning_R_session_Info_KilliBrain_Aging_Per_Cell_Type.txt", sep =""))
sessionInfo()
sink()


