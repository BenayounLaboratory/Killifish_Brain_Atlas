setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Microglia/UCell')
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

library(readxl)
library(bitops)

library(UCell)

# 2024-02-26
# Run Ucell on desired gene lists for microglia

# 2024-03-11
# Add additional microglia lists
# remove more general aging/senescence lists (will be used for All QC cell types)

# 2024-07-01
# Add TIM (terminally inflammatory microglia) Millet et al 2024

# 2025-03-18
# Add mouse M1/M2 enriched
# Add Kang, J Neuro, 2024 mouse microglia aging lists

####################################################################################################################
# 0. Load annotated Seurat Objects and extract microglia

# Import final annotation
load('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Atlas_Analysis/2023-08-23_Seurat_object_with_Manual_annotation_FINAL.RData')
killi.brain.clean
# An object of class Seurat 
# 21160 features across 209939 samples within 1 assay 
# Active assay: RNA (21160 features, 5000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony

# subset microglia
killi.brain.microglia <- subset(killi.brain.clean, subset = Cell_Identity %in% "Microglia")    
killi.brain.microglia
# An object of class Seurat 
# 21160 features across 3466 samples within 1 assay 
# Active assay: RNA (21160 features, 5000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony

save(killi.brain.microglia, file = paste(Sys.Date(),"Seurat_object_MICROGLIA.RData",sep = "_"))
####################################################################################################################


####################################################################################################################
# 1.Process/recluster microglia cells

# load('2024-02-19_Seurat_object_MICROGLIA.RData')

# use the same pipeline as whole object (including harmony integration)
# for consistency but with just microglia as input

# Run dimensionality reduction PCA
killi.brain.microglia <- RunPCA(killi.brain.microglia, npcs = 40)

# Determine the ‘dimensionality’ of the dataset
# Approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
pdf(paste0(Sys.Date(), "_Killi_Fish_AgingBrain_MICROGLIA_AllCohorts_elbowplot.pdf"), height = 5, width= 6)
ElbowPlot(killi.brain.microglia, ndims = 40)
dev.off()

################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# To give us an idea of the number of PCs needed to be included:
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - The principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- killi.brain.microglia[["pca"]]@stdev / sum(killi.brain.microglia[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 35

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 11

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 11

# Based on these metrics, first 14 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_Killi_Fish_AgingBrain_MICROGLIA_AllCohorts_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

##############
# Run Harmony integration
# https://portals.broadinstitute.org/harmony/SeuratV3.html
# https://hbctraining.github.io/scRNA-seq_online/lessons/06a_integration_harmony.html

p1 <- DimPlot(object = killi.brain.microglia, reduction = "pca", pt.size = 2, group.by = "Batch", raster = T , raster.dpi = c(1024, 1024))
p2 <- VlnPlot(object = killi.brain.microglia, features = "PC_1", group.by = "Batch",  pt.size = 0, raster = F)

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_MICROGLIA_AllCohorts_Pre_HarmonyIntegration_QC_Plots.pdf", sep = "_"), height = 5, width = 8)
p1+p2
dev.off()

# Use RNA 
# https://github.com/immunogenomics/harmony/issues/24
killi.brain.microglia <- RunHarmony(killi.brain.microglia,  
                                    group.by.vars    = c("Batch","Strain", "Sex"),
                                    reduction        = "pca", 
                                    assay.use        = "RNA", 
                                    reduction.save   = "harmony",
                                    theta            = c(3,3,3),
                                    max.iter.harmony = 10) 
# Harmony converged after 3 iterations

# Check batches are integrated in the first 2 dimensions after Harmony.
p1 <- DimPlot(object = killi.brain.microglia, reduction = "harmony", pt.size = 2, group.by = "Batch", raster = T , raster.dpi = c(1024, 1024))
p2 <- VlnPlot(object = killi.brain.microglia, features = "harmony_1", group.by = "Batch",  pt.size = 0, raster = F)

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_MICROGLIA_AllCohorts_POST_HarmonyIntegration_QC_Plots.pdf", sep = "_"), height = 5, width = 8)
p1+p2
dev.off()

# run umap and clustering
set.seed(1234)
killi.brain.microglia <- RunUMAP(killi.brain.microglia, reduction = "harmony", assay = "RNA", dims = 1:11)
killi.brain.microglia <- FindNeighbors(object = killi.brain.microglia, reduction = "harmony")
killi.brain.microglia <- FindClusters(killi.brain.microglia, resolution = c(0.2, 0.4, 0.6))
# Number of communities: 4
# Number of communities: 7
# Number of communities: 9


################  Plot summary UMAPs post harmony ################

# QC UMAPs
pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_MICROGLIA_AllCohorts_Singlets_UMAP_color_by_Sex_postHarmony.pdf", sep = "_"), height = 5, width = 6)
DimPlot(killi.brain.microglia, reduction = "umap", group.by = "Sex", pt.size	= 3, shuffle = T,
        cols = c(alpha("deeppink", alpha = 0.5 ), alpha("deepskyblue", alpha = 0.5 ) ), raster = T, raster.dpi = c(1024, 1024) )
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_MICROGLIA_AllCohorts_Singlets_UMAP_color_by_Strain_postHarmony.pdf", sep = "_"), height = 5, width = 6)
DimPlot(killi.brain.microglia, reduction = "umap", group.by = "Strain", pt.size	= 3, shuffle = T,
        cols = c(alpha("goldenrod1", alpha = 0.5 ), alpha("forestgreen", alpha = 0.5 ) ), raster = T, raster.dpi = c(1024, 1024) )
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_MICROGLIA_AllCohorts_Singlets_UMAP_color_by_Age_postHarmony.pdf", sep = "_"), height = 5, width = 6)
DimPlot(killi.brain.microglia, reduction = "umap", group.by = "Age_Group", pt.size	= 3, shuffle = T,
        cols = c(alpha("firebrick3", alpha = 0.5 ), alpha("darkturquoise", alpha = 0.5 ), alpha("darkorange", alpha = 0.5 ), alpha("darkolivegreen2", alpha = 0.5 ) ),
        raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_MICROGLIA_AllCohorts_Singlets_UMAP_color_by_Batch_postHarmony.pdf", sep = "_"), height = 5, width = 6)
DimPlot(killi.brain.microglia, reduction = "umap", group.by = "Batch", pt.size	= 3, shuffle = T, raster = T, raster.dpi = c(1024, 1024))
dev.off()

#### Plot clusters
pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_MICROGLIA_AllCohorts_Singlets_UMAP_postHarmony_snn_0.2.pdf", sep = "_"), height = 5, width = 6)
DimPlot(killi.brain.microglia, reduction = "umap", group.by = "RNA_snn_res.0.2", pt.size	= 3, shuffle = T, raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_MICROGLIA_AllCohorts_Singlets_UMAP_postHarmony_snn_0.4.pdf", sep = "_"), height = 5, width = 6)
DimPlot(killi.brain.microglia, reduction = "umap", group.by = "RNA_snn_res.0.4", pt.size	= 3, shuffle = T, raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_MICROGLIA_AllCohorts_Singlets_UMAP_postHarmony_snn_0.6.pdf", sep = "_"), height = 5, width = 6)
DimPlot(killi.brain.microglia, reduction = "umap", group.by = "RNA_snn_res.0.6", pt.size	= 3, shuffle = T, raster = T, raster.dpi = c(1024, 1024))
dev.off()

# subset per strains
killi.microglia.grz <- subset(killi.brain.microglia, subset = Strain %in% "GRZ")    # 1,661 cells
killi.microglia.zmz <- subset(killi.brain.microglia, subset = Strain %in% "ZMZ")    # 1,805 cells
save(killi.microglia.grz, killi.microglia.zmz, killi.brain.microglia,
     file = paste(Sys.Date(),"Seurat_objects_MICROGLIA_SPLIT_PER_STRAIN.RData",sep = "_"))
####################################################################################################################


####################################################################################################################
# 2.Prepare gene lists of interest
# 
load("2024-02-19_Seurat_object_MICROGLIA.RData")
load("2024-02-26_Seurat_objects_MICROGLIA_SPLIT_PER_STRAIN.RData")

# Read in BLAST homology file for killifish/mouse (best mouse hit to killifish to get conversion)
mouse.homol <- read.csv("../../../Mouse_alignment/2022-10-11_Mouse_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-3_Minimal_HOMOLOGY_TABLE_REV.txt", sep = "\t", header = T)


### read gene lists
mglia.aging.Hajdarovic <- read.csv('../GeneLists/Hajdarovic-2022-TableS4.csv')
mglia.aging.Ocanas     <- read_xlsx('../GeneLists/Ocanas-2023-12974_2023_2870_MOESM8_ESM.xlsx')
mglia.lists.Ocanas     <- read_xlsx('../GeneLists/Ocanas-2023-12974_2023_2870_MOESM10_ESM.xlsx', skip = 1)
mglia.millet.TIM       <- read.table('../GeneLists/Millet_2024_TIM_markers.txt', sep = "\t", header = T)
mph.polar              <- qusage::read.gmt('/Volumes/BB_HQ_1/Immune_sex_dimorphism_Aging/Macrophages/Macrophage_datasets_for_Comparison/M1_M2_macrophages/GSE103958_Bulk/DEseq2/2021-11-18_GSE103958_BMDM_polarization_genes_lists.gmt')
mglia.lists.Kang       <- read_xlsx('../GeneLists/Kang2024_12974_2024_3130_MOESM3_ESM.xlsx', skip = 2, sheet = "Additional File 1b")

NAMES.mglia.glists <- c("mouse_Aging_UP"       ,
                        "mouse_Aging_DOWN"     ,
                        "mouse_M1_enriched"    ,
                        "mouse_M2_enriched"    ,
                        "Homeostatic"          ,
                        "DAM"                  ,
                        "IRM"                  ,
                        "LDAM"                 ,
                        "TIM"                  )

mglia.glists.mouse        <- vector(mode = "list", length = length(NAMES.mglia.glists))
names(mglia.glists.mouse) <- NAMES.mglia.glists

## sumarize microglia aging
Aging_up_Hajdarovic   <- mglia.aging.Hajdarovic[bitAnd(mglia.aging.Hajdarovic$cell_type %in% "Microglia", mglia.aging.Hajdarovic$updown %in% "Upregulated with age")>0,]$gene
Aging_down_Hajdarovic <- mglia.aging.Hajdarovic[bitAnd(mglia.aging.Hajdarovic$cell_type %in% "Microglia", mglia.aging.Hajdarovic$updown %in% "Downregulated with age")>0,]$gene
Aging_up_Ocanas       <- mglia.aging.Ocanas[mglia.aging.Ocanas$`up with age (1), down with age (-1), discordant (0)` == 1,]$SYMBOL
Aging_down_Ocanas     <- mglia.aging.Ocanas[mglia.aging.Ocanas$`up with age (1), down with age (-1), discordant (0)` == -1,]$SYMBOL
Aging_up_Kang         <- mglia.lists.Kang[mglia.lists.Kang$log2FoldChange > 0,]$Gene
Aging_down_Kang       <- mglia.lists.Kang[mglia.lists.Kang$log2FoldChange < 0,]$Gene

mglia.up.summary <- aggregate(c(Aging_up_Hajdarovic, Aging_up_Ocanas, Aging_up_Kang), 
                              by = list("up_genes" = c(Aging_up_Hajdarovic, Aging_up_Ocanas, Aging_up_Kang)),
                              FUN = length)
mglia.dn.summary <- aggregate(c(Aging_down_Hajdarovic, Aging_down_Ocanas, Aging_down_Kang), 
                              by = list("dn_genes" = c(Aging_down_Hajdarovic, Aging_down_Ocanas, Aging_down_Kang)),
                              FUN = length)

mglia.glists.mouse$mouse_Aging_UP       <- mglia.up.summary$up_genes[mglia.up.summary$x >= 2] # 246
mglia.glists.mouse$mouse_Aging_DOWN     <- mglia.dn.summary$dn_genes[mglia.dn.summary$x >= 2] # 35


## macrophage polarization
mglia.glists.mouse$mouse_M1_enriched     <- mph.polar$BMDM_M1_enriched
mglia.glists.mouse$mouse_M2_enriched     <- mph.polar$BMDM_M2_enriched


## microglia state signatures
mglia.glists.mouse$Homeostatic           <- mglia.lists.Ocanas$Homeostatic
mglia.glists.mouse$DAM                   <- mglia.lists.Ocanas$`DAM/MgND/ARM`
mglia.glists.mouse$IRM                   <- mglia.lists.Ocanas$IRM
mglia.glists.mouse$LDAM                  <- mglia.lists.Ocanas$LDAM
mglia.glists.mouse$TIM                   <- mglia.millet.TIM$gene



### function for conversion
killify_gsets <- function(my.gset.list, homol.table = mouse.homol) {
  
  # prepare list
  my.gset.list.killi <- vector(length = length(my.gset.list), mode = "list")
  names(my.gset.list.killi) <- names(my.gset.list) 
  
  # grab all killi homologs to mouse genes in gene set
  for (i in 1:length(my.gset.list)) {
    my.gset.list.killi[[i]] <- unique(homol.table$Nfur_Symbol[homol.table$Mmu_Symbol %in% my.gset.list[[i]]])
  }
  return(my.gset.list.killi)
}

mglia.glists.killi <- killify_gsets(mglia.glists.mouse)

save(mglia.glists.killi,
     file = paste(Sys.Date(),"Microglia_Gene_lists_killified.RData",sep = "_"))
####################################################################################################################


####################################################################################################################
# 3. Run UCell on gene lists of interest

# killi.microglia.grz
# killi.microglia.zmz

##################################################
####################   GRZ   #####################
##################################################

killi.microglia.grz <- AddModuleScore_UCell(killi.microglia.grz, features=mglia.glists.killi, name=NULL)
killi.microglia.grz$Age_Group <- factor(killi.microglia.grz$Age_Group, levels = c("Y", "M", "O"))

grz.dot <- DotPlot(killi.microglia.grz, features = rev(names(mglia.glists.killi)), group.by = "Age_Group", col.min = -0.5,col.max = 0.5) 
grz.dot <- grz.dot + scale_colour_gradient2(low = "#333399", mid = "lightgrey", high = "#CC3333")
grz.dot <- grz.dot + coord_flip() + scale_size(range = c(3, 10)) + scale_size_area(limits=c(25,100))
grz.dot

pdf(paste0(Sys.Date(),"_DotPLot_Microglia_GeneLists_UCell_Scores_GRZ_Microglia.pdf"), width = 5, height = 4)
grz.dot
dev.off()

grz.stats           <- data.frame(matrix(0,9,2))
colnames(grz.stats) <- c("Gene_List","Wilcoxon_P")
grz.stats$Gene_List <- names(mglia.glists.killi)

test.mouse_Aging_UP        <- wilcox.test(mouse_Aging_UP        ~  Age_Group, data = killi.microglia.grz@meta.data[killi.microglia.grz$Age_Group %in% c("Y", "O"),]) # p-value = 1.207e-11
test.mouse_Aging_DOWN      <- wilcox.test(mouse_Aging_DOWN      ~  Age_Group, data = killi.microglia.grz@meta.data[killi.microglia.grz$Age_Group %in% c("Y", "O"),]) # p-value = 1.882e-05
test.mouse_M1_enriched     <- wilcox.test(mouse_M1_enriched     ~  Age_Group, data = killi.microglia.grz@meta.data[killi.microglia.grz$Age_Group %in% c("Y", "O"),]) # p-value = 2.186e-05
test.mouse_M2_enriched     <- wilcox.test(mouse_M2_enriched     ~  Age_Group, data = killi.microglia.grz@meta.data[killi.microglia.grz$Age_Group %in% c("Y", "O"),]) # p-value = 5.987e-07
test.Homeostatic           <- wilcox.test(Homeostatic           ~  Age_Group, data = killi.microglia.grz@meta.data[killi.microglia.grz$Age_Group %in% c("Y", "O"),]) # p-value = 0.01849
test.DAM                   <- wilcox.test(DAM                   ~  Age_Group, data = killi.microglia.grz@meta.data[killi.microglia.grz$Age_Group %in% c("Y", "O"),]) # p-value = 7.819e-12
test.IRM                   <- wilcox.test(IRM                   ~  Age_Group, data = killi.microglia.grz@meta.data[killi.microglia.grz$Age_Group %in% c("Y", "O"),]) # p-value < 2.2e-16
test.LDAM                  <- wilcox.test(LDAM                  ~  Age_Group, data = killi.microglia.grz@meta.data[killi.microglia.grz$Age_Group %in% c("Y", "O"),]) # p-value = 0.000344
test.TIM                   <- wilcox.test(TIM                   ~  Age_Group, data = killi.microglia.grz@meta.data[killi.microglia.grz$Age_Group %in% c("Y", "O"),]) # p-value = 0.000344

grz.stats[1,]$Wilcoxon_P <- test.mouse_Aging_UP   $p.value
grz.stats[2,]$Wilcoxon_P <- test.mouse_Aging_DOWN $p.value
grz.stats[3,]$Wilcoxon_P <- test.mouse_M1_enriched$p.value
grz.stats[4,]$Wilcoxon_P <- test.mouse_M2_enriched$p.value
grz.stats[5,]$Wilcoxon_P <- test.Homeostatic      $p.value
grz.stats[6,]$Wilcoxon_P <- test.DAM              $p.value
grz.stats[7,]$Wilcoxon_P <- test.IRM              $p.value
grz.stats[8,]$Wilcoxon_P <- test.LDAM             $p.value
grz.stats[9,]$Wilcoxon_P <- test.TIM              $p.value

write.table(grz.stats, file = paste0(Sys.Date(),"_Enrichment_Wilcoxon_MICROGLIA_GeneLists_UCell_Scores_GRZ.txt"), sep = "\t", row.names = F, quote = F)


##################################################
####################   ZMZ   #####################
##################################################

killi.microglia.zmz <- AddModuleScore_UCell(killi.microglia.zmz, features=mglia.glists.killi, name=NULL)
killi.microglia.zmz$Age_Group <- factor(killi.microglia.zmz$Age_Group, levels = c("Y", "M", "O", "G"))

zmz.dot <- DotPlot(killi.microglia.zmz, features = rev(names(mglia.glists.killi)), group.by = "Age_Group", col.min = -0.5,col.max = 0.5) 
zmz.dot <- zmz.dot + scale_colour_gradient2(low = "#333399", mid = "lightgrey", high = "#CC3333")
zmz.dot <- zmz.dot + coord_flip() + scale_size(range = c(3, 10)) + scale_size_area(limits=c(25,100))
zmz.dot

pdf(paste0(Sys.Date(),"_DotPLot_Microglia_GeneLists_UCell_Scores_ZMZ_Microglia.pdf"), width = 5.5, height = 4)
zmz.dot
dev.off()

zmz.stats           <- data.frame(matrix(0,9,2))
colnames(zmz.stats) <- c("Gene_List","Wilcoxon_P")
zmz.stats$Gene_List <- names(mglia.glists.killi)

test.mouse_Aging_UP        <- wilcox.test(mouse_Aging_UP        ~  Age_Group, data = killi.microglia.zmz@meta.data[killi.microglia.zmz$Age_Group %in% c("Y", "G"),]) # p-value < 2.2e-16
test.mouse_Aging_DOWN      <- wilcox.test(mouse_Aging_DOWN      ~  Age_Group, data = killi.microglia.zmz@meta.data[killi.microglia.zmz$Age_Group %in% c("Y", "G"),]) # p-value = 0.2215
test.mouse_M1_enriched     <- wilcox.test(mouse_M1_enriched     ~  Age_Group, data = killi.microglia.zmz@meta.data[killi.microglia.zmz$Age_Group %in% c("Y", "G"),]) # p-value =  3.022e-13
test.mouse_M2_enriched     <- wilcox.test(mouse_M2_enriched     ~  Age_Group, data = killi.microglia.zmz@meta.data[killi.microglia.zmz$Age_Group %in% c("Y", "G"),]) # p-value = 9.443e-11
test.Homeostatic           <- wilcox.test(Homeostatic           ~  Age_Group, data = killi.microglia.zmz@meta.data[killi.microglia.zmz$Age_Group %in% c("Y", "G"),]) # p-value = 0.08374226
test.DAM                   <- wilcox.test(DAM                   ~  Age_Group, data = killi.microglia.zmz@meta.data[killi.microglia.zmz$Age_Group %in% c("Y", "G"),]) # p-value = 6.503306e-15
test.IRM                   <- wilcox.test(IRM                   ~  Age_Group, data = killi.microglia.zmz@meta.data[killi.microglia.zmz$Age_Group %in% c("Y", "G"),]) # p-value < 2.2e-16
test.LDAM                  <- wilcox.test(LDAM                  ~  Age_Group, data = killi.microglia.zmz@meta.data[killi.microglia.zmz$Age_Group %in% c("Y", "G"),]) # p-value = 2.671e-07
test.TIM                   <- wilcox.test(TIM                   ~  Age_Group, data = killi.microglia.zmz@meta.data[killi.microglia.zmz$Age_Group %in% c("Y", "G"),]) # p-value = 0.0009139

zmz.stats[1,]$Wilcoxon_P <- test.mouse_Aging_UP   $p.value
zmz.stats[2,]$Wilcoxon_P <- test.mouse_Aging_DOWN $p.value
zmz.stats[3,]$Wilcoxon_P <- test.mouse_M1_enriched$p.value
zmz.stats[4,]$Wilcoxon_P <- test.mouse_M2_enriched$p.value
zmz.stats[5,]$Wilcoxon_P <- test.Homeostatic      $p.value
zmz.stats[6,]$Wilcoxon_P <- test.DAM              $p.value
zmz.stats[7,]$Wilcoxon_P <- test.IRM              $p.value
zmz.stats[8,]$Wilcoxon_P <- test.LDAM             $p.value
zmz.stats[9,]$Wilcoxon_P <- test.TIM              $p.value

write.table(zmz.stats, file = paste0(Sys.Date(),"_Enrichment_Wilcoxon_MICROGLIA_GeneLists_UCell_Scores_ZMZ.txt"), sep = "\t", row.names = F, quote = F)
#############################################################

#############################################################
# 4. on microglia subclusters
##################################################
####################   GRZ   #####################
##################################################

killi.microglia.clust <- AddModuleScore_UCell(killi.brain.microglia, features=mglia.glists.killi, name=NULL)

clust.dot <- DotPlot(killi.microglia.clust, features = rev(names(mglia.glists.killi)), group.by = "RNA_snn_res.0.2", col.min = 0,col.max = 1, scale = F) 
clust.dot <- clust.dot + scale_colour_gradient2(low = "white", mid = "lightgrey", high = "#CC3333")
clust.dot <- clust.dot + coord_flip() + scale_size(range = c(3, 10)) + scale_size_area(limits=c(10,100))
clust.dot

pdf(paste0(Sys.Date(),"_DotPLot_Microglia_GeneLists_UCell_Scores_SubClusters_Microglia.pdf"), width = 5, height = 4)
clust.dot
dev.off()

#############################################################
sink(file = paste(Sys.Date(),"_Ucell_Microglia_Scoring_session_Info.txt", sep =""))
sessionInfo()
sink()