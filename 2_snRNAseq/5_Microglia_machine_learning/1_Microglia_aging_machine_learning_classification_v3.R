setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Microglia/Machine_Learning')
options(stringsAsFactors = F)

library(Seurat)          #

library(caret)           #
library(randomForest)    # random forest
library(glmnet)          # Elasticnet
library(gbm)             # GBM
library(MLeval)          #

library(ggplot2)         #
library(scales)          #
library(vioplot)         #
library(matrixStats)     # for column multiplication

# 2025-07-01
# ML to predict biological age from microglia across strains
# Microglia is the cell type with clearest changes across 'omic' layers
#
# Regression is too difficult of a training problem with *only* 3 time points
# Try a classification young/old and use probability of old
# similar to the CellBiAge paper from Webb lab
# see Yu et al, Cell Reports, 2023 [CellBiAge]

# 2025-07-04
# Training/Testing/Validation split works well, but seems confusing to non ML peeps
# for the purpose of the paper, we will rerun doing only Training/Testing 
# This will also improve power

#####################################################################################################################
#### 1. Load annotated microglia Seurat object
# this is the scaled/normalized datasets on microglia only
load('../UCell/2024-02-26_Seurat_objects_MICROGLIA_SPLIT_PER_STRAIN.RData')
killi.brain.microglia
# An object of class Seurat 
# 21160 features across 3466 samples within 1 assay 
# Active assay: RNA (21160 features, 5000 variable features)
# 3 layers present: counts, data, scale.data
# 3 dimensional reductions calculated: pca, umap, harmony

# scale microglia data gene expression (see Yu et al, Cell Reports, 2023 [CellBiAge])
killi.brain.microglia <- ScaleData(object = killi.brain.microglia, vars.to.regress = c("nCount_RNA","nFeature_RNA", "percent.mito", "Phase", "Batch"), features = rownames(killi.brain.microglia))
killi.brain.microglia
# An object of class Seurat 
# 21160 features across 3466 samples within 1 assay 
# Active assay: RNA (21160 features, 5000 variable features)
# 3 layers present: counts, data, scale.data
# 3 dimensional reductions calculated: pca, umap, harmony

save(killi.brain.microglia, file = paste0(Sys.Date(),"_Microglia_SCALED_Seurat_objects_Killi_Brain_for_ML.RData"))

#### For training/testing partition
## We will use GRZ pilot, set 2 and 3 cohorts for training/testing vs. set 4 validation set (independent animals)
## This allows for all groups to be represented in each partition derived from independent animals/cohorts
## ZMZ will only be used as experimental check (for potentially slower aging)

load('2025-07-02_Microglia_SCALED_Seurat_objects_Killi_Brain_for_ML.RData')

# subset per strains
killi.mglia.grz           <- subset(killi.brain.microglia, subset = Strain %in% "GRZ") # 1661 cells
killi.mglia.zmz           <- subset(killi.brain.microglia, subset = Strain %in% "ZMZ") # 1805 cells

# check cohorts
table(killi.mglia.grz@meta.data$Batch)
# Pilot Set_2 Set_3 Set_4 
#    95   505   735   326 

save(killi.mglia.zmz, killi.mglia.grz, file = paste0(Sys.Date(),"_Microglia_Seurat_objects_GRZ_ZMZ__Killi_Brain_for_ML.RData"))
#####################################################################################################################


#####################################################################################################################
#### 2. Prepare data for ML

# for ML training: will train a RF model for each cell type, need to sample cells in each age-sex group/cell type
# get 100 cells per group
my.BioGroupForSampling                            <- unique(killi.mglia.grz@meta.data$Group)
length(my.BioGroupForSampling) # [1] 6

table(killi.mglia.grz@meta.data$Group)
# GRZ_M_F GRZ_M_M GRZ_O_F GRZ_O_M GRZ_Y_F GRZ_Y_M 
#     200     192     284     691     153     141 

# exclude middle time point
my.middle <- grep("_M_",my.BioGroupForSampling)
my.BioGroupForSampling <- my.BioGroupForSampling[-my.middle]
length(my.BioGroupForSampling) # [1] 4

# create object to store info of randomly sampled cells from each group
train.cell.list        <- vector(mode = "list", length = length(my.BioGroupForSampling))
names(train.cell.list) <- my.BioGroupForSampling

# set seed for reproducibility
set.seed(2017)

# get training cells
for (i in 1:length(my.BioGroupForSampling)){
  
  # identify cells from the group
  # only get cells from young and old (not middle) group
  cell.names <- rownames(killi.mglia.grz@meta.data)[killi.mglia.grz@meta.data$Group %in% my.BioGroupForSampling[i]]
  
  # sample desired cell number from each group (125)
  samp.cells <- sample(cell.names, size = 125)
  
  # extract training and testing
  train.cell.list[[i]] <- samp.cells
}

length(unlist(train.cell.list)) # 500

# get training and testing:
train.cells    <- unlist(train.cell.list)                               #  500
test.cells.all <- setdiff(colnames(killi.mglia.grz),train.cells)        # 1161                
length(grep("GRZ_Y_",test.cells.all)) # 44
length(grep("GRZ_M_",test.cells.all)) # 392
length(grep("GRZ_O_",test.cells.all)) # 725

# create subsetted Seurat objects by extracting cells
killi.mglia.grz.train    <- subset(killi.mglia.grz , cells = train.cells    ) # 21160 features across  500 samples within 1 assay 
killi.mglia.grz.test.all <- subset(killi.mglia.grz , cells = test.cells.all ) # 21160 features across 1161 samples within 1 assay 

# select most variable BUT robustly expressed genes in **training set** (so feature selection remains agnostic to testing/validation data)
killi.mglia.grz.train    <- FindVariableFeatures(killi.mglia.grz.train, nfeatures = 3000)
killi.mglia.variable     <- VariableFeatures(killi.mglia.grz.train)

# get average expression over groups
Idents(killi.mglia.grz.train) <- "Group"
killi.mglia.grz.train.av      <- AverageExpression(killi.mglia.grz.train, assays = "RNA")

# select only genes robustly detected in all cell types
my.robust <-  apply(killi.mglia.grz.train.av$RNA>0.25,1,sum) == 4 # detected in groups
sum(my.robust) # 5125

# genes used as features are both variable AND robustly expressed
select.genes <- intersect(killi.mglia.variable,names(my.robust)[my.robust])
select.genes <- select.genes[!grepl("NotFur", select.genes)]                 # only genes, not TEs
length(select.genes) # 826

# binarize gene expression (to make more robust to batch effects)
# do it on scaled data
#    see Yu et al, Cell Reports, 2023 [CellBiAge]
tmp.gene.grz.train    <- apply(as.matrix(killi.mglia.grz.train@assays$RNA@scale.data   [select.genes,] ) > 0, 1, as.numeric)
tmp.gene.grz.test.all <- apply(as.matrix(killi.mglia.grz.test.all@assays$RNA@scale.data[select.genes,] ) > 0, 1, as.numeric)
tmp.gene.grz.test.yo  <- apply(as.matrix(killi.mglia.grz.test.yo@assays$RNA@scale.data [select.genes,] ) > 0, 1, as.numeric)
tmp.gene.zmz          <- apply(as.matrix(killi.mglia.zmz@assays$RNA@scale.data         [select.genes,] ) > 0, 1, as.numeric)

# attach age/sex covariate
featmat.train    <- cbind(killi.mglia.grz.train@meta.data   [,c("Age_Group" ,"Sex")], tmp.gene.grz.train    )
featmat.test.all <- cbind(killi.mglia.grz.test.all@meta.data[,c("Age_Group" ,"Sex")], tmp.gene.grz.test.all )
featmat.test.yo  <- cbind(killi.mglia.grz.test.yo@meta.data [,c("Age_Group" ,"Sex")], tmp.gene.grz.test.yo  )
featmat.zmz      <- cbind(killi.mglia.zmz@meta.data         [,c("Age_Group" ,"Sex")], tmp.gene.zmz          )

featmat.train   $Age_Group <- as.factor(featmat.train   $Age_Group)
featmat.test.all$Age_Group <- as.factor(featmat.test.all$Age_Group)
featmat.test.yo $Age_Group <- as.factor(featmat.test.yo $Age_Group)
featmat.zmz     $Age_Group <- as.factor(featmat.zmz     $Age_Group)

featmat.train   $Sex <- as.factor(featmat.train   $Sex )
featmat.test.all$Sex <- as.factor(featmat.test.all$Sex )
featmat.test.yo $Sex <- as.factor(featmat.test.yo $Sex )
featmat.zmz     $Sex <- as.factor(featmat.zmz     $Sex )

dim(featmat.train)
# 500 828
dim(featmat.test.all)
# 1161  828
dim(featmat.test.yo)
# 974 828
dim(featmat.zmz)
# 1805  828


# save feature matrices
save(featmat.train   ,
     featmat.test.all,
     featmat.zmz     ,
     file = paste0(Sys.Date(),"_Training_Testing_ZMZ_Data_Killi_Microglia_Aging_for_ML.RData"))
#####################################################################################################################


#####################################################################################################################
#### 3. Run  ML

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%% A. run elastic net model

# set seed for reproducibility
set.seed(2017)

# use 10-fold cross-validation to build the model
my.ctrl.opt.enet <- trainControl(method          = "cv",
                                 number          = 10,
                                 allowParallel   = TRUE,
                                 verbose         = F,
                                 classProbs      = TRUE,
                                 savePredictions = TRUE,
                                 summaryFunction = twoClassSummary)

# Create parameter search grid
fineGrid.enet <- expand.grid(alpha    = c(0.01,0.05,seq(0.1,1,0.3)),
                             lambda   = c(0.0001,0.001,0.01, 0.1))

# train model with caret train function
my.enet.fit       <- train( Age_Group ~ .,
                            data       =  featmat.train,
                            method     = "glmnet",
                            importance = TRUE,
                            trControl  = my.ctrl.opt.enet,
                            tuneGrid   = fineGrid.enet,
                            metric     = "ROC")
my.enet.fit
# ROC was used to select the optimal model using the largest value.
# The final values used for the model were alpha = 0.1 and lambda = 0.1.
#   alpha  lambda  ROC       Sens   Spec 
#     0.10   1e-01   0.85120  0.744  0.784

save(my.enet.fit, file = paste0(Sys.Date(),"_ElasticNet_model_Aging_Microglia.RData"))


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%% B. run RF model

# set seed for reproducibility
set.seed(2017)

# use 10-fold cross-validation to build the model
my.ctrl.opt.rf <- trainControl(method          = "cv",
                               number          = 10,
                               allowParallel   = TRUE,
                               verbose         = F,
                               classProbs      = TRUE,
                               savePredictions = TRUE,
                               summaryFunction = twoClassSummary)

# Create parameter search grid
fineGrid.rf <- expand.grid(mtry = c(25, 50, 100, 200, 300))

# train model with caret train function
my.rf.fit       <- train( Age_Group ~ .,
                          data       = featmat.train,
                          method     = "rf",
                          importance = TRUE,
                          trControl  =  my.ctrl.opt.rf,
                          tuneGrid   =  fineGrid.rf,
                          metric     = "ROC")
my.rf.fit
# ROC was used to select the optimal model using the largest value.
# The final value used for the model was mtry = 50
# mtry  ROC      Sens   Spec 
#      50   0.85032  0.732  0.788

save(my.rf.fit, file = paste0(Sys.Date(),"_RF_models_Aging_Microglia.RData"))
#####################################################################################################################


#####################################################################################################################
#### 4. OOB performance analysis

# Grab 10-fold CV OOB training accuracies

results.ml  <- resamples(list("ENET" = my.enet.fit,
                              "RF"   = my.rf.fit))

# summary of model differences
my.model.summaries <- summary(results.ml)

my.10cv.data     <- my.model.summaries$values
my.10cv.data.roc <- my.10cv.data[,grep("ROC",colnames(my.10cv.data))]

my.median.roc     <- apply(my.10cv.data.roc,2,median)
my.median.roc
# ENET~ROC   RF~ROC 
#   0.8472   0.8664 

pdf(paste0(Sys.Date(),"_Machine_Learning_OOB_Accuracy_boxplots_10CV_beeswarm_ENET_RF.pdf"), height = 3.5, width = 5)
par(oma=c(0.1,4,0.1,0.1))
boxplot(rev(my.10cv.data.roc), las = 1, ylim = c(0.4,1), horizontal = T, xlab = "10 Fold-CV AUROC", outline = F, col = c("firebrick","coral"), main = "OOB AUROC")
beeswarm::beeswarm(rev(my.10cv.data.roc), add = T, horizontal = T, pch = 16, cex = 0.7)
abline(v = 0.5, col = "red", lty = "dashed")
dev.off()

# get ROC curves
# run roc calculation on final model
ROC.analysis <- evalm(list(my.enet.fit,my.rf.fit),
                      gnames=c('eNET','RF'),
                      rlinethick=0.8,fsize=10,
                      plots='r',
                      title = names(my.rf.fit)[i], 
                      cols = c("firebrick","coral"))

pdf(paste0(Sys.Date(),"_ROC_Microglia_ENET_RF_training.pdf") )
print(ROC.analysis$roc)
dev.off()

#####################################################################################################################


#####################################################################################################################
#### 5. testing and validation performance analysis

######################################################################
# create output object (for balanced accuracy)
perf.mat           <- matrix(NA,2,2)
colnames(perf.mat) <- c("ENET", "RF")
rownames(perf.mat) <- c("Training_OOB","Testing")

# grab only young/old from testing data
tmp.test           <- featmat.test.all[featmat.test.all$Age_Group != "M",]
tmp.test$Age_Group <- factor(tmp.test$Age_Group)

################ ENET ################
# get training (OOB over folds)
enet.confus.train                <- confusionMatrix(my.enet.fit$pred$pred, my.enet.fit$pred$obs)
perf.mat["Training_OOB", "ENET"] <- enet.confus.train$byClass[11]
# Balanced Accuracy 
# 0.7470833 

# get testing
enet.predict.age.test            <- predict(my.enet.fit, tmp.test)
enet.confus.test                 <- confusionMatrix(enet.predict.age.test, as.factor(tmp.test$Age_Group))
perf.mat["Testing", "ENET"]      <- enet.confus.test$byClass[11]
# Balanced Accuracy 
# 0.8180251 


################ RF ################
# get training (OOB over folds)
rf.confus.train                <- confusionMatrix(my.rf.fit$pred$pred, my.rf.fit$pred$obs)
perf.mat["Training_OOB", "RF"] <- rf.confus.train$byClass[11]
# Balanced Accuracy 
# 0.7636 

# get testing
rf.predict.age.test            <- predict(my.rf.fit, tmp.test)
rf.confus.test                 <- confusionMatrix(rf.predict.age.test, as.factor(tmp.test$Age_Group))
perf.mat["Testing", "RF"]      <- rf.confus.test$byClass[11]
# Balanced Accuracy 
# 0.7666928 

pdf(paste(Sys.Date(),"Machine_Learning_Training_Testing_Accuracy_dotchart_ENET_RF.pdf", sep =""), height = 4, width = 4)
dotchart(t(perf.mat), xlim = c(0.5,1), pch = 16, pt.cex = 2, main ="ML Sex / Balanced Accuracy", col = c("firebrick","coral"))
abline(v = 0.5, col = "red", lty = "dashed")
dev.off()
#####################################################################################################################


#####################################################################################################################
#### 6. get class probabilities

my.grz.cols <- c("deeppink"    ,
                 "deeppink3"   ,
                 "deeppink4"   ,
                 "deepskyblue" ,
                 "deepskyblue3",
                 "deepskyblue4")

my.zmz.cols <- c("deeppink"    ,
                 "deeppink3"   ,
                 "deeppink4"   ,
                 "magenta4",
                 "deepskyblue" ,
                 "deepskyblue3",
                 "deepskyblue4",
                 "royalblue4")

# grab testing and validation data + zmz data
data.test           <- featmat.test.all
data.zmz            <- featmat.zmz


# use models to get class prob
enet.predict.age.test           <- predict(my.enet.fit, data.test, type = "prob")
rf.predict.age.test             <- predict(my.rf.fit, data.test, type = "prob")
rownames(enet.predict.age.test) <- rownames(data.test)
data.test$prob_Old_enet         <- enet.predict.age.test$O
data.test$prob_Old_rf           <- rf.predict.age.test$O
data.test$Group                 <- factor(paste0(data.test$Age_Group,data.test$Sex), levels = c("YF", "MF", "OF", "YM", "MM", "OM"))

# get ZMZ data
enet.predict.age.zmz            <- predict(my.enet.fit, data.zmz, type = "prob")
rf.predict.age.zmz              <- predict(my.rf.fit, data.zmz, type = "prob")
rownames(enet.predict.age.zmz)  <- rownames(data.zmz)
data.zmz$prob_Old_enet          <- enet.predict.age.zmz$O
data.zmz$prob_Old_rf            <- rf.predict.age.zmz$O
data.zmz$Group                  <- factor(paste0(data.zmz$Age_Group,data.zmz$Sex), levels = c("YF", "MF", "OF","GF", "YM", "MM", "OM","GM"))


data.test$Age_Group  <- factor(data.test$Age_Group , levels = c("Y","M","O"))
data.zmz$Age_Group   <- factor(data.zmz$Age_Group  , levels = c("Y","M","O", "G"))


#### Plot violins of class probs
pdf(paste0(Sys.Date(),"_Microglia_OClass_probability_violin_ENET.pdf"), height = 4, width = 7)
par(mfrow = c(1,2))
vioplot(prob_Old_enet  ~ Group, 
        data = data.test, 
        col  = my.grz.cols,
        main = "Testing (ENET)",
        ylim = c(0,1),
        ylab = "Predicted probability of 'old' age",las = 1)
abline(h = 0.5, col = "red", lty = "dashed")

vioplot(prob_Old_enet  ~ Group, 
        data = data.zmz, 
        col  = my.zmz.cols,
        main = "ZMZ (ENET)",
        ylim = c(0,1),
        ylab = "Predicted probability of 'old' age",las = 1)
abline(h = 0.5, col = "red", lty = "dashed")
dev.off()


pdf(paste0(Sys.Date(),"_Microglia_OClass_probability_violin_RF.pdf"), height = 4, width = 7)
par(mfrow = c(1,2))
vioplot(prob_Old_rf  ~ Group, 
        data = data.test, 
        col  = my.grz.cols,
        main = "Testing (RF)",
        ylim = c(0,1),
        ylab = "Predicted probability of 'old' age",las = 1)
abline(h = 0.5, col = "red", lty = "dashed")

vioplot(prob_Old_rf  ~ Group, 
        data = data.zmz, 
        col  = my.zmz.cols,
        main =  "ZMZ (RF)",
        ylim = c(0,1),
        ylab = "Predicted probability of 'old' age",las = 1)
abline(h = 0.5, col = "red", lty = "dashed")
dev.off()

#####################################################################################################################

#####################################################################################################################
#### 7. test differences in class probabilities

### ENET
krusk.test.enet <- kruskal.test(prob_Old_enet  ~ Group , data = data.test)
dunn.test.enet  <- FSA::dunnTest(prob_Old_enet  ~ Group, data = data.test, method = "holm")

krusk.zmz.enet <- kruskal.test(prob_Old_enet  ~ Group , data = data.zmz)
dunn.zmz.enet  <- FSA::dunnTest(prob_Old_enet  ~ Group, data = data.zmz, method = "holm")

# RF
krusk.test.rf <- kruskal.test(prob_Old_rf  ~ Group , data = data.test)
dunn.test.rf  <- FSA::dunnTest(prob_Old_rf  ~ Group, data = data.test, method = "holm")

krusk.zmz.rf <- kruskal.test(prob_Old_rf  ~ Group , data = data.zmz)
dunn.zmz.rf  <- FSA::dunnTest(prob_Old_rf  ~ Group, data = data.zmz, method = "holm")

# create output object 
stats.mat           <- data.frame(matrix(NA,2,2))
colnames(stats.mat) <- c("ENET", "RF")
rownames(stats.mat) <- c("Testing", "ZMZ")

stats.mat["Testing"   ,"ENET"]  <- krusk.test.enet$p.value
stats.mat["ZMZ"       ,"ENET"]  <- krusk.zmz.enet$p.value 

stats.mat["Testing"   ,"RF"]    <- krusk.test.rf$p.value
stats.mat["ZMZ"       ,"RF"]    <- krusk.zmz.rf$p.value 

stats.mat
#                 ENET           RF
# Testing 3.879894e-51 4.591291e-44
# ZMZ     9.236497e-89 3.138487e-82

write.table( stats.mat ,
             file = paste0(Sys.Date(),"_Kruskal_Wallis_Test_results_class_probs_by_model_and_set.txt"),
             sep = "\t", quote = F)

### format Dunn's
dunn.test.enet $res$Model <- "ENET"
dunn.zmz.enet  $res$Model <- "ENET"
dunn.test.rf   $res$Model <- "RF"
dunn.zmz.rf    $res$Model <- "RF"

dunn.test.enet $res$Set <- "Testing"
dunn.zmz.enet  $res$Set <- "ZMZ-1001"
dunn.test.rf   $res$Set <- "Testing"
dunn.zmz.rf    $res$Set <- "ZMZ-1001"


write.table( rbind(dunn.test.enet $res,
                   dunn.zmz.enet  $res,
                   dunn.test.rf   $res,
                   dunn.zmz.rf    $res)  ,
             file = paste0(Sys.Date(),"_Dunn_Test_results_Holm_class_probs_by_model_and_set.txt"),
             sep = "\t", quote = F, row.names = F)


################# get group medians
# create output object 
medians.mat.grz           <- data.frame(matrix(NA,6,2))
colnames(medians.mat.grz) <- c("ENET", "RF")
rownames(medians.mat.grz) <- levels(data.test$Group)

for (g in rownames(medians.mat.grz)) {
  
  medians.mat.grz[g, "ENET"] <- median(data.test$prob_Old_enet[data.test$Group %in% g])
  medians.mat.grz[g, "RF"  ] <- median(data.test$prob_Old_rf[data.test$Group %in% g])
  
}

medians.mat.grz
#         ENET    RF
# YF 0.2892790 0.386
# MF 0.5989945 0.535
# OF 0.6585276 0.578
# YM 0.3131016 0.421
# MM 0.3769478 0.442
# OM 0.7634600 0.620

###
medians.mat.zmz           <- data.frame(matrix(NA,8,2))
colnames(medians.mat.zmz) <- c("ENET", "RF")
rownames(medians.mat.zmz) <- levels(data.zmz$Group)


for (g in rownames(medians.mat.zmz)) {
  
  medians.mat.zmz[g, "ENET"] <- median(data.zmz$prob_Old_enet[data.zmz$Group %in% g])
  medians.mat.zmz[g, "RF"  ] <- median(data.zmz$prob_Old_rf[data.zmz$Group %in% g])
  
}

medians.mat.zmz
#         ENET    RF
# YF 0.2922394 0.398
# MF 0.4252286 0.427
# OF 0.5299102 0.498
# GF 0.7018567 0.602
# YM 0.2768653 0.406
# MM 0.5395775 0.502
# OM 0.7312132 0.604
# GM 0.6123025 0.538

#####################################################################################################################

#######################
sink(file = paste(Sys.Date(),"Machine_Learning_R_session_Info_KilliBrain_Aging_Microglia.txt", sep =""))
sessionInfo()
sink()


