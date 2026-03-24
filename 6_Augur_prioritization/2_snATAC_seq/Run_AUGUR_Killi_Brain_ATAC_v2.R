setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/snATAC_Brain_Aging_Meta/Downstream_Analyses/AUGUR')
options(stringsAsFactors = F)
options (future.globals.maxSize = 32000 * 1024^2)


#### Packages
library(Seurat)        # single cell general package
library(Signac)        # scATAC processing

library(GenomeInfoDb)  # for genome info
library(GenomicRanges) # for genome info
library(seqinr)        # for genome info

library(Augur)
library(viridis)

# 2024-04-22
# run AUGUR
# run vanilla R due to memory constraints

# 2025-06-25
# run separately aging split by sex

######################################################################
# Load annotated object
load("../../Signac/2023-08-23_Killi_brain_ATAC_2Cohorts_FiltNorm_Singlets_ANNOTATEDfromRNA.RData")
brain.atac.singlets
# An object of class Seurat
# 148827 features across 27171 samples within 2 assays
# Active assay: RNA (22135 features, 0 variable features)
# 1 other assay present: ATAC
# 2 dimensional reductions calculated: lsi, umap

brain.atac.singlets.f <- subset(brain.atac.singlets, subset = sex %in% "F")
brain.atac.singlets.m <- subset(brain.atac.singlets, subset = sex %in% "M")

brain.atac.singlets.f
# An object of class Seurat 
# 148827 features across 14698 samples within 2 assays 
# Active assay: RNA (22135 features, 0 variable features)
# 2 layers present: counts, data
# 1 other assay present: ATAC
# 2 dimensional reductions calculated: lsi, umap

brain.atac.singlets.m
# An object of class Seurat 
# 148827 features across 12473 samples within 2 assays 
# Active assay: RNA (22135 features, 0 variable features)
# 2 layers present: counts, data
# 1 other assay present: ATAC
# 2 dimensional reductions calculated: lsi, umap

brain.atac.singlets.f$age_group <- ifelse(brain.atac.singlets.f$age_weeks > 10,"Old","Young")
brain.atac.singlets.m$age_group <- ifelse(brain.atac.singlets.m$age_weeks > 10,"Old","Young")

###########################################################################################


###########################################################################################
# Run AUGUR based on transferred cell type labels

###### A. based on gene activity scores [peak scores doesn't run]

augur.brain.atac.f <-  calculate_auc(brain.atac.singlets.f,
                                     cell_type_col = "predicted.id", 
                                     label_col = "age_group",
                                     n_threads = 1,
                                     min_cells = 25)

augur.brain.atac.f
## # A tibble: 14 × 2
##    cell_type                      auc
##    <chr>                        <dbl>
##  1 OPCs                         0.600
##  2 Vascular_smooth_muscle_cells 0.572
##  3 Neurons_misc_2               0.571
##  4 PV_interneurons              0.565
##  5 Neurons_misc_4               0.557
##  6 GABAergic_neurons            0.540
##  7 Oligodendrocytes             0.536
##  8 Granule_Excitatory_Neurons   0.525
##  9 Microglia                    0.524
## 10 Neurons_misc_3               0.523
## 11 NSPCs                        0.504
## 12 Neurons_misc_1               0.503
## 13 Astrocytes_Radial_Glia       0.492
## 14 Ependymal_cells              0.490


augur.brain.atac.m <-  calculate_auc(brain.atac.singlets.m,
                                     cell_type_col = "predicted.id", 
                                     label_col = "age_group",
                                     n_threads = 1,
                                     min_cells = 25)

augur.brain.atac.m
## $AUC
## # A tibble: 13 × 2
##    cell_type                      auc
##    <chr>                        <dbl>
##  1 OPCs                         0.624
##  2 PV_interneurons              0.582
##  3 Neurons_misc_4               0.580
##  4 GABAergic_neurons            0.562
##  5 Microglia                    0.540
##  6 Granule_Excitatory_Neurons   0.540
##  7 Vascular_smooth_muscle_cells 0.538
##  8 Oligodendrocytes             0.529
##  9 NSPCs                        0.518
## 10 Astrocytes_Radial_Glia       0.517
## 11 Ependymal_cells              0.501
## 12 Neurons_misc_2               0.489
## 13 Neurons_misc_3               0.477

save(augur.brain.atac.f,augur.brain.atac.m, file = paste0(Sys.Date(),"_Augur_Killi_brain_ATAC_2Cohorts_GeneActivityScore_AGING_ONLY.RData"))



load('2025-06-25_Augur_Killi_brain_ATAC_2Cohorts_GeneActivityScore_AGING_ONLY.RData')

# pdf(paste0(Sys.Date(),"_Augur_Killi_brain_ATAC_2Cohorts_FEMALE_UMAP_Red_Blue_GeneActivityScore.pdf"), width = 3, height = 3)
# plot_umap(augur.brain.atac.f, brain.atac.singlets.f, cell_type_col = "Cell_Identity", palette = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","lightgrey","#CCCCFF","#9999FF","#333399")))(50))
# dev.off()
# 
# pdf(paste0(Sys.Date(),"_Augur_Killi_brain_ATAC_2Cohorts_MALE_UMAP_Red_Blue_GeneActivityScore.pdf"), width = 3, height = 3)
# plot_umap(augur.brain.atac.m, brain.atac.singlets.m, cell_type_col = "Cell_Identity", palette = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","lightgrey","#CCCCFF","#9999FF","#333399")))(50))
# dev.off()

pdf(paste0(Sys.Date(),"_Augur_Killi_brain_ATAC_2Cohorts_Lollipop_GeneActivityScore_AGING_FEMALE.pdf"), width = 3, height = 3)
plot_lollipop(augur.brain.atac.f)
dev.off()

pdf(paste0(Sys.Date(),"_Augur_Killi_brain_ATAC_2Cohorts_Lollipop_GeneActivityScore_AGING_MALE.pdf"), width = 3, height = 3)
plot_lollipop(augur.brain.atac.m)
dev.off()
###################################################################################################################


#######################
sink(file = paste(Sys.Date(),"_Killi_brain_ATAC_by_sex_AUGUR_session_Info.txt", sep =""))
sessionInfo()
sink()
