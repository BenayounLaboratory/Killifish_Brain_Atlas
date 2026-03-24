setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Atlas_Analysis')
options(stringsAsFactors = F)
options (future.globals.maxSize = 32000 * 1024^2)

# Load packages
library('Seurat')    # 
library(sctransform) # 
library(clustree)    # 
library(scales)      # 
library(harmony)     #


##################################################################################
##### 2023-07-04
# Perform joint integration/analysis across cohorts
# use decontX cleaned up data per cohort (with 0.25 score filter)
# - Pilot cohort (GRZ 6w F & M) (had the manual annotation)
# - Brain aging cohort # 1 (ZMZ only)
# - Brain aging cohort # 2 (GRZ, ZMZ)
# - Brain aging cohort # 3 (GRZ, ZMZ)
# - Brain aging cohort # 4 (GRZ only)
#
# Try Harmony for cross dataset integration
# https://hbctraining.github.io/scRNA-seq_online/lessons/06a_integration_harmony.html
# https://portals.broadinstitute.org/harmony/SeuratV3.html
# 
# SCT keeps crashing even on Vanilla R (runs out of memory)
# Will have to use more "basic" NormalizeData version to handle dataset
# see: https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/64
##################################################################################



#####################################################################################################################
#### 1. Load Cleaned up Seurat Objects and merge data

# load up cleaned up singlets
load("../Preprocessing/DecontX/Pilot_Cohort/2023-07-03_Killifish_Brain_Pilot_Cohort_Seurat_object_SINGLETS_ONLY.RData")
load("../Preprocessing/DecontX/Cohort_1/2023-07-03_Killifish_Brain_Cohort1_Seurat_object_SINGLETS_ONLY.RData")
load("../Preprocessing/DecontX/Cohort_2/2023-07-03_Killifish_Brain_Cohort2_Seurat_object_SINGLETS_ONLY.RData")
load("../Preprocessing/DecontX/Cohort_3/2023-07-03_Killifish_Brain_Cohort3_Seurat_object_SINGLETS_ONLY.RData")
load("../Preprocessing/DecontX/Cohort_4/2023-07-04_Killifish_Brain_Cohort4_Seurat_object_SINGLETS_ONLY.RData")

# Merge Seurat objects
killi.brain.sg.all <- merge(killi.singlets.cp,
                            y =  c(killi.singlets.c1,
                                   killi.singlets.c2,
                                   killi.singlets.c3,
                                   killi.singlets.c4),
                            project = "10x_killi_brain")
killi.brain.sg.all
# An object of class Seurat 
# 53290 features across 209939 samples within 2 assays 
# Active assay: SCT (25456 features, 0 variable features)
#  1 other assay present: RNA

# bring RNA as main assay again
DefaultAssay(killi.brain.sg.all) <- "RNA"

table(killi.brain.sg.all@meta.data$Batch)
# Pilot Set_1 Set_2 Set_3 Set_4 
# 11864 37963 70599 66352 23161 

# will need to clean up and rerun SCT on merged data
# https://github.com/satijalab/seurat/issues/2662
killi.brain.sg.all[['SCT']] <- NULL
killi.brain.sg.all
# An object of class Seurat 
# 27834 features across 209939 samples within 1 assay 
# Active assay: RNA (27834 features, 0 variable features)

rm(killi.singlets.cp, killi.singlets.c1, killi.singlets.c2, killi.singlets.c3, killi.singlets.c4)

#### Filter genes with expression that is too sparse across all batches
# https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day2/scRNA_Workshop-PART2.html
min.value = 0
min.cells = 250
genes.use <- rownames(killi.brain.sg.all@assays$RNA)
num.cells <- Matrix::rowSums(killi.brain.sg.all@assays$RNA@counts > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])

# remove low/null genes
killi.brain.sg.all <- subset(killi.brain.sg.all, features = genes.use)
killi.brain.sg.all
# An object of class Seurat 
# 21160 features across 209939 samples within 1 assay 
# Active assay: RNA (21160 features, 0 variable features)

# remove irrelevant columns (previous SCT on unmerged objects, DoubletFinder column and old clustering info)
sct.cols <- grep("SCT",colnames(killi.brain.sg.all@meta.data))
rm.cols  <- c(grep("seurat_clusters",colnames(killi.brain.sg.all@meta.data))  ,
              grep("DoubletFinder"  ,colnames(killi.brain.sg.all@meta.data))  ,
              grep("scds_hybrid"    ,colnames(killi.brain.sg.all@meta.data))  )

killi.brain.sg.all@meta.data <- killi.brain.sg.all@meta.data[,-c(sct.cols,rm.cols)]

save(killi.brain.sg.all, file = paste0(Sys.Date(),"_Killi_Fish_AgingBrain_All_Cohorts_Seurat_object_Merge_preNorm.RData"))
################################################################################################################################################################


################################################################################################################################################################
#### 2. Scaling the data and removing unwanted sources of variation

# Since SCT runs out of memory, will use the basic alternative pipeline
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
killi.brain.sg.all <- NormalizeData(killi.brain.sg.all, normalization.method = "LogNormalize", scale.factor = 10000)
killi.brain.sg.all <- FindVariableFeatures(killi.brain.sg.all, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
killi.brain.var.top20 <- head(VariableFeatures(killi.brain.sg.all), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(killi.brain.sg.all)
plot2 <- LabelPoints(plot = plot1, points = killi.brain.var.top20, repel = TRUE)

pdf(paste0(Sys.Date(),"_Killi_Fish_AgingBrain_All_Cohorts_top20_variable_genes_QC_plot.pdf"), width = 12, height = 5)
plot1 + plot2
dev.off()

# Scaling the data (for all genes => runs out of memory, let's only do the variable genes)
killi.brain.clean <- ScaleData(object = killi.brain.sg.all, vars.to.regress = c("nCount_RNA","nFeature_RNA", "percent.mito", "Phase", "Batch"))
killi.brain.clean
# An object of class Seurat 
# 21160 features across 209939 samples within 1 assay 
# Active assay: RNA (21160 features, 5000 variable features)

table(killi.brain.clean@meta.data$Batch, killi.brain.clean@meta.data$Group)
#         GRZ_M_F GRZ_M_M GRZ_O_F GRZ_O_M GRZ_Y_F GRZ_Y_M ZMZ_G_F ZMZ_G_M ZMZ_M_F ZMZ_M_M ZMZ_O_F ZMZ_O_M ZMZ_Y_F ZMZ_Y_M
#   Pilot       0       0       0       0    5905    5959       0       0       0       0       0       0       0       0
#   Set_1       0       0       0       0       0       0    3947    4472    5108    4567    4358    5036    4736    5739
#   Set_2    6574    4894    4701    6390    5880    6400    4760    4253    4360    3958    4508    4364    4366    5191
#   Set_3    5990    5451    5225    5551    5988    6596    1498    1734    5143    7205    4468    2121    4147    5235
#   Set_4    4289    3950    3665    4142    3786    3329       0       0       0       0       0       0       0       0

table(killi.brain.clean@meta.data$Group)
# GRZ_M_F GRZ_M_M GRZ_O_F GRZ_O_M GRZ_Y_F GRZ_Y_M ZMZ_G_F ZMZ_G_M ZMZ_M_F ZMZ_M_M ZMZ_O_F ZMZ_O_M ZMZ_Y_F ZMZ_Y_M 
#     16853   14295   13591   16083   21559   22284   10205   10459   14611   15730   13334   11521   13249   16165 

save(killi.brain.clean, file = paste0(Sys.Date(),"_Killi_Fish_AgingBrain_AllCohorts_Seurat_object_logNorm.RData"))

# clean up unnecessary objects from memory
rm(killi.brain.sg.all)
###############################################################################################################################################################


###############################################################################################################################################################
##### 3. Run dimensionality reduction

# Run dimensionality reduction PCA
killi.brain.clean <- RunPCA(killi.brain.clean, npcs = 40)

# Determine the ‘dimensionality’ of the dataset
# Approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
pdf(paste0(Sys.Date(), "_Killi_Fish_AgingBrain_AllCohorts_elbowplot.pdf"), height = 5, width= 6)
ElbowPlot(killi.brain.clean, ndims = 40)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- killi.brain.clean[["pca"]]@stdev / sum(killi.brain.clean[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 34

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 19

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 19

# Based on these metrics, first 14 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_Killi_Fish_AgingBrain_AllCohorts_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()
###############################################################################


# Calculate UMAP
killi.brain.clean <- RunUMAP(killi.brain.clean, dims = 1:pcs)
killi.brain.clean
# An object of class Seurat 
# 21160 features across 209939 samples within 1 assay 
# Active assay: RNA (21160 features, 5000 variable features)
#  2 dimensional reductions calculated: pca, umap

save(killi.brain.clean, file = paste0(Sys.Date(),"_Killi_Fish_AgingBrain_AllCohorts_Seurat_object_logNorm_with_UMAP.RData"))

################  Plot summary UMAPs  ################

# QC UMAPs
pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_color_by_Sex.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "Sex", pt.size	= 2, shuffle = T,
        cols = c(alpha("deeppink", alpha = 0.5 ), alpha("deepskyblue", alpha = 0.5 ) ) , raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_color_by_Strain.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "Strain", pt.size	= 2, shuffle = T,
        cols = c(alpha("goldenrod1", alpha = 0.5 ), alpha("forestgreen", alpha = 0.5 ) ) , raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_color_by_Age.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "Age_Group", pt.size	= 2, shuffle = T,
        cols = c(alpha("firebrick3", alpha = 0.5 ), alpha("darkturquoise", alpha = 0.5 ), alpha("darkorange", alpha = 0.5 ), alpha("darkolivegreen2", alpha = 0.5 ) ),
        raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_color_by_Group.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "Group", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024) )
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_color_by_Batch.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "Batch", pt.size	= 2, shuffle = T, raster = T , raster.dpi = c(1024, 1024) )
dev.off()

#### Plot markers
pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_Markers_UMAP.pdf", sep = "_"), height = 5, width = 6.5)
FeaturePlot(killi.brain.clean, features = c("olig1","olig2","mpz")                               , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
FeaturePlot(killi.brain.clean, features = c("marco","csf1r", "ptprc")                            , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # mph/microglia
FeaturePlot(killi.brain.clean, features = c("s100b", "slc1a2")                                   , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # astrocyte/radial glia
FeaturePlot(killi.brain.clean, features = c("rbfox3", "map2", "eno2", "ncam1")                   , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # mature neuron
FeaturePlot(killi.brain.clean, features = c("dcx")                                               , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # immature neuron
FeaturePlot(killi.brain.clean, features = c( "fat2", "neurod1", "eomes", "pax6")                 , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Granule excitatory neuron
FeaturePlot(killi.brain.clean, features = c("pvalb")                                             , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # pvalb interneurons
FeaturePlot(killi.brain.clean, features = c("gad1", "gad2",  "LOC107384443", "LOC107391088")     , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # GABAergic neurons
FeaturePlot(killi.brain.clean, features = c("LOC107386767")                                      , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # dopaminergic neurons
FeaturePlot(killi.brain.clean, features = c("slc17a6","LOC107381463")                            , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # glutamatergic neurons
FeaturePlot(killi.brain.clean, features = c("LOC107386535",  "clu")                              , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # ependymal cells
FeaturePlot(killi.brain.clean, features = c("sox5",  "sox2")                                     , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # NSPCs
dev.off()
###############################################################################################################################################################


###############################################################################################################################################################
##### 4. Harmony integration to lessen impact of batch effects
# https://portals.broadinstitute.org/harmony/SeuratV3.html
# https://hbctraining.github.io/scRNA-seq_online/lessons/06a_integration_harmony.html

p1 <- DimPlot(object = killi.brain.clean, reduction = "pca", pt.size = 2, group.by = "Batch", raster = T , raster.dpi = c(1024, 1024))
p2 <- VlnPlot(object = killi.brain.clean, features = "PC_1", group.by = "Batch",  pt.size = 0, raster = F)

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Pre_HarmonyIntegration_QC_Plots.pdf", sep = "_"), height = 5, width = 8)
p1+p2
dev.off()

# Run Harmony
# Use RNA not SCT, since SCT cannot run on an object of this size
# https://github.com/immunogenomics/harmony/issues/24
killi.brain.clean <- RunHarmony(killi.brain.clean,  
                                group.by.vars    = c("Batch","Strain"),
                                reduction        = "pca", 
                                assay.use        = "RNA", 
                                reduction.save   = "harmony",
                                theta            = c(3,3),
                                max.iter.harmony = 20) 
# Harmony converged after 8 iterations


###################################
# Check batches are integrated in the first 2 dimensions after Harmony.
p1 <- DimPlot(object = killi.brain.clean, reduction = "harmony", pt.size = 2, group.by = "Batch", raster = T , raster.dpi = c(1024, 1024))
p2 <- VlnPlot(object = killi.brain.clean, features = "harmony_1", group.by = "Batch",  pt.size = 0, raster = F)

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_POST_HarmonyIntegration_QC_Plots.pdf", sep = "_"), height = 5, width = 8)
p1+p2
dev.off()


####################################################################################
# Generate a UMAP/clustering derived from harmony embeddings instead of PCs:
killi.brain.clean <- RunUMAP(killi.brain.clean, reduction = "harmony", assay = "RNA", dims = 1:19)
killi.brain.clean <- FindNeighbors(object = killi.brain.clean, reduction = "harmony")
killi.brain.clean <- FindClusters(killi.brain.clean, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0))
# Number of communities: 14
# Number of communities: 21
# Number of communities: 24
# Number of communities: 28
# Number of communities: 31

killi.brain.clean
# An object of class Seurat 
# 21160 features across 209939 samples within 1 assay 
# Active assay: RNA (21160 features, 5000 variable features)
#  3 dimensional reductions calculated: pca, umap, harmony


################  Plot summary UMAPs post harmony ################

# QC UMAPs
pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_color_by_Sex_postHarmony.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "Sex", pt.size	= 2, shuffle = T,
        cols = c(alpha("deeppink", alpha = 0.5 ), alpha("deepskyblue", alpha = 0.5 ) ), raster = T, raster.dpi = c(1024, 1024) )
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_color_by_Strain_postHarmony.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "Strain", pt.size	= 2, shuffle = T,
        cols = c(alpha("goldenrod1", alpha = 0.5 ), alpha("forestgreen", alpha = 0.5 ) ), raster = T, raster.dpi = c(1024, 1024) )
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_color_by_Age_postHarmony.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "Age_Group", pt.size	= 2, shuffle = T,
        cols = c(alpha("firebrick3", alpha = 0.5 ), alpha("darkturquoise", alpha = 0.5 ), alpha("darkorange", alpha = 0.5 ), alpha("darkolivegreen2", alpha = 0.5 ) ),
raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_color_by_Group_postHarmony.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "Group", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_color_by_Batch_postHarmony.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "Batch", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024))
dev.off()

#### Plot clusters
pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_postHarmony_snn_0.2.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "RNA_snn_res.0.2", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_postHarmony_snn_0.4.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "RNA_snn_res.0.4", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_postHarmony_snn_0.6.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "RNA_snn_res.0.6", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_postHarmony_snn_0.8.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "RNA_snn_res.0.8", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_postHarmony_snn_1.0.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "RNA_snn_res.1", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024))
dev.off()


#### Plot markers
pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_Markers_UMAP_postHarmony.pdf", sep = "_"), height = 5, width = 6.5)
FeaturePlot(killi.brain.clean, features = c("olig1","olig2","mpz")                               , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
FeaturePlot(killi.brain.clean, features = c("marco","csf1r", "ptprc")                            , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # mph/microglia
FeaturePlot(killi.brain.clean, features = c("s100b", "slc1a2")                                   , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # astrocyte/radial glia
FeaturePlot(killi.brain.clean, features = c("rbfox3", "map2", "eno2", "ncam1")                   , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # mature neuron
FeaturePlot(killi.brain.clean, features = c("dcx")                                               , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # immature neuron
FeaturePlot(killi.brain.clean, features = c("fat2", "neurod1", "eomes", "pax6")                  , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Granule excitatory neuron
FeaturePlot(killi.brain.clean, features = c("pvalb")                                             , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # pvalb interneurons
FeaturePlot(killi.brain.clean, features = c("gad1", "gad2",  "LOC107384443", "LOC107391088")     , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # GABAergic neurons
FeaturePlot(killi.brain.clean, features = c("LOC107386767")                                      , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # dopaminergic neurons
FeaturePlot(killi.brain.clean, features = c("slc17a6","LOC107381463")                            , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # glutamatergic neurons
FeaturePlot(killi.brain.clean, features = c("LOC107386535",  "clu")                              , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # ependymal cells
FeaturePlot(killi.brain.clean, features = c("sox5",  "sox2")                                     , pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # NSPCs
dev.off()

save(killi.brain.clean, file = paste0(Sys.Date(),"_Killi_Fish_AgingBrain_AllCohorts_Seurat_object_logNorm_with_UMAP_Post_Harmony.RData"))
###############################################################################################################################################################

###############################################################################################################################################################
sink(file = paste(Sys.Date(),"_Killi_Brain_Atlas_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()

