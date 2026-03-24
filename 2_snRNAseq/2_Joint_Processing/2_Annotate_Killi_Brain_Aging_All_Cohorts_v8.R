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

##########  Cell identity annotation packages ##########  
# https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html
library(scSorter)

# # https://github.com/pcahan1/singleCellNet
# library(singleCellNet)
# SingleCellNet runs out of memory

#https://github.com/TianLab-Bioinfo/scMAGIC
library(scMAGIC)
library(sparseMatrixStats) ### used to mousify the gene expression table for annotation
########################################################

#####################################################################################################################
# 2023-07-28
# annotate cells with scSorter/PanglaoDB
#
# 2023-08-03
# Try splitting by cohort to run scMAGIC
#
# 2023-08-07
# reparse other annotations from output (do not rerun as the expression has not changed)
#
# 2023-08-18
# Run additional clustering resolutions and cluster tree to help annotation (helps gad2+ cells to stand apart)
# run annotation using Raj et al 2018
#
# 2023-08-21
# recompile all intermediates using parsed output
# run marker genes on new resolution
#####################################################################################################################


#####################################################################################################################
#### 1. Load Cleaned up Seurat Objects and merge data

############ load integrated object
load('2023-07-04_Killi_Fish_AgingBrain_AllCohorts_Seurat_object_logNorm_with_UMAP_Post_Harmony.RData')
killi.brain.clean
# An object of class Seurat 
# 21160 features across 209939 samples within 1 assay 
# Active assay: RNA (21160 features, 5000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony


############ Add clustering granularity to get all important cell types resolved
### Add 2 additional resolutions
killi.brain.clean <- FindClusters(killi.brain.clean, resolution = c(1.2,1.5))
# Number of communities: 36
# Number of communities: 42

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_postHarmony_snn_1.2.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "RNA_snn_res.1.2", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_postHarmony_snn_1.5.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "RNA_snn_res.1.5", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_postHarmony_snn_1.2_LABEL.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "RNA_snn_res.1.2", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024), label = T)
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_postHarmony_snn_1.5_LABEL.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "RNA_snn_res.1.5", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024), label = T)
dev.off()

pdf(paste(Sys.Date(),"Killi_Fish_AgingBrain_AllCohorts_Singlets_UMAP_postHarmony_snn_0.6_LABEL.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "RNA_snn_res.0.6", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024), label = T)
dev.off()

############ Dotplot of markers by cluster

# res 0.6
pdf(paste(Sys.Date(),"Killifish_Brain_Dotplot_Known_Cell_type_markers_by_RNA_snn_res.0.6_CLEAN.pdf", sep = "_"), height = 7, width = 14)
DotPlot(killi.brain.clean,
        features = c("olig1","olig2","mpz",                           # Oligodendrocyte/OPC
                     "marco","csf1r", "ptprc",                        # mph/microglia
                     "s100b", "slc1a2",                               # astrocyte/radial glia
                     "rbfox3", "map2", "eno2", "ncam1",               # mature neuron
                     "dcx",                                           # immature neuron
                     "fat2", "neurod1", "eomes", "pax6",              # Granule excitatory neuron
                     "pvalb",                                         # pvalb interneurons
                     "LOC107373010",                                  # SST interneurons ### LOC107373010 sst1.1 somatostatin 1, tandem duplicate 1 
                     "gad1", "gad2",  "LOC107384443", "LOC107391088", # GABAergic neurons
                     "LOC107386767",                                  # dopaminergic neurons
                     "slc17a6","LOC107381463",                        # glutamatergic neurons
                     "LOC107386535",  "clu",                          # ependymal cells
                     "sox5",  "sox2"),                                # NSPCs
        group.by = "RNA_snn_res.0.6") + RotatedAxis()
dev.off()

# res 1.0
pdf(paste(Sys.Date(),"Killifish_Brain_Dotplot_Known_Cell_type_markers_by_RNA_snn_res.1_CLEAN.pdf", sep = "_"), height = 7, width = 14)
DotPlot(killi.brain.clean,
        features = c("olig1","olig2","mpz",                           # Oligodendrocyte/OPC
                     "marco","csf1r", "ptprc",                        # mph/microglia
                     "s100b", "slc1a2",                               # astrocyte/radial glia
                     "rbfox3", "map2", "eno2", "ncam1",               # mature neuron
                     "dcx",                                           # immature neuron
                     "fat2", "neurod1", "eomes", "pax6",              # Granule excitatory neuron
                     "pvalb",                                         # pvalb interneurons
                     "LOC107373010",                                  # SST interneurons ### LOC107373010 sst1.1 somatostatin 1, tandem duplicate 1 
                     "gad1", "gad2",  "LOC107384443", "LOC107391088", # GABAergic neurons
                     "LOC107386767",                                  # dopaminergic neurons
                     "slc17a6","LOC107381463",                        # glutamatergic neurons
                     "LOC107386535",  "clu",                          # ependymal cells
                     "sox5",  "sox2"),                                # NSPCs
        group.by = "RNA_snn_res.1") + RotatedAxis()
dev.off()

# res 1.2
pdf(paste(Sys.Date(),"Killifish_Brain_Dotplot_Known_Cell_type_markers_by_RNA_snn_res.1.2_CLEAN.pdf", sep = "_"), height = 8, width = 14)
DotPlot(killi.brain.clean,
        features = c("olig1","olig2","mpz",                           # Oligodendrocyte/OPC
                     "marco","csf1r", "ptprc",                        # mph/microglia
                     "s100b", "slc1a2",                               # astrocyte/radial glia
                     "rbfox3", "map2", "eno2", "ncam1",               # mature neuron
                     "dcx",                                           # immature neuron
                     "fat2", "neurod1", "eomes", "pax6",              # Granule excitatory neuron
                     "pvalb",                                         # pvalb interneurons
                     "LOC107373010",                                  # SST interneurons ### LOC107373010 sst1.1 somatostatin 1, tandem duplicate 1 
                     "gad1", "gad2",  "LOC107384443", "LOC107391088", # GABAergic neurons
                     "LOC107386767",                                  # dopaminergic neurons
                     "slc17a6","LOC107381463",                        # glutamatergic neurons
                     "LOC107386535",  "clu",                          # ependymal cells
                     "sox5",  "sox2"),                                # NSPCs
        group.by = "RNA_snn_res.1.2") + RotatedAxis()
dev.off()

# res 1.5
pdf(paste(Sys.Date(),"Killifish_Brain_Dotplot_Known_Cell_type_markers_by_RNA_snn_res.1.5_CLEAN.pdf", sep = "_"), height = 8, width = 14)
DotPlot(killi.brain.clean,
        features = c("olig1","olig2","mpz",                           # Oligodendrocyte/OPC
                     "marco","csf1r", "ptprc",                        # mph/microglia
                     "s100b", "slc1a2",                               # astrocyte/radial glia
                     "rbfox3", "map2", "eno2", "ncam1",               # mature neuron
                     "dcx",                                           # immature neuron
                     "fat2", "neurod1", "eomes", "pax6",              # Granule excitatory neuron
                     "pvalb",                                         # pvalb interneurons
                     "LOC107373010",                                  # SST interneurons ### LOC107373010 sst1.1 somatostatin 1, tandem duplicate 1 
                     "gad1", "gad2",  "LOC107384443", "LOC107391088", # GABAergic neurons
                     "LOC107386767",                                  # dopaminergic neurons
                     "slc17a6","LOC107381463",                        # glutamatergic neurons
                     "LOC107386535",  "clu",                          # ependymal cells
                     "sox5",  "sox2"),                                # NSPCs
        group.by = "RNA_snn_res.1.5") + RotatedAxis()
dev.off()
################################################################################################################################################################



################################################################################################################################################################
##### 2. Try using scSorter "gates" to annotate cell types
# https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html

# Use mouse markers (PanglaoDB) and Zebrafish markers (Raj et al 2018)

####################################################################################################################################################
#####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       Mouse markers from PanglaoDB        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#####
####################################################################################################################################################

# Load Mouse/Killifish homology table
# (best mouse hit to killifish to get conversion)
mouse.homol <- read.csv("../Mouse_alignment/2022-10-11_Mouse_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-3_Minimal_HOMOLOGY_TABLE_REV.txt", sep = "\t", header = T)

# ###%%%%%%%%%%%%%%%%%%%### A. Filter and prep PanglaoDB marker data  ###%%%%%%%%%%%%%%%%%###

############ Section 1 - Generate annotation object for scSorter ############
# Prepare the annotation for scSorter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Read in marker file
panglao.markers <- read.csv("../Atlas_Info/Mouse_markers/PanglaoDB_markers_27_Mar_2020_BRAIN.txt", sep = "\t", header = T)
panglao.markers$official.gene.symbol <-  firstup(tolower( panglao.markers$official.gene.symbol )) 

panglao.cell.types <- unique(panglao.markers$cell.type)

# create marker list object to summarize all marker genes
my.marker.list         <- vector(mode = "list",length = length(panglao.cell.types))
names(my.marker.list)  <- panglao.cell.types

for (i in 1:length(panglao.cell.types)) {
  
  # get unique markers from panglaoDB
  my.markers <- unique(panglao.markers$official.gene.symbol[panglao.markers$cell.type %in% panglao.cell.types[i]])
  
  # get closest Nfur corresponding gene (unique as well)
  my.marker.list[[i]] <- unique(mouse.homol$Nfur_Symbol[mouse.homol$Mmu_Symbol %in% my.markers])
  
}

names(my.marker.list)
# [1] "Adrenergic neurons"               "Anterior pituitary gland cells"   "Astrocytes"                       "Bergmann glia"                   
# [5] "Cajal-Retzius cells"              "Cholinergic neurons"              "Choroid plexus cells"             "Dopaminergic neurons"            
# [9] "Ependymal cells"                  "GABAergic neurons"                "Glutaminergic neurons"            "Glycinergic neurons"             
# [13] "Immature neurons"                 "Interneurons"                     "Meningeal cells"                  "Microglia"                       
# [17] "Motor neurons"                    "Neural stem/precursor cells"      "Neuroblasts"                      "Neuroendocrine cells"            
# [21] "Neurons"                          "Noradrenergic neurons"            "Oligodendrocyte progenitor cells" "Oligodendrocytes"                
# [25] "Pinealocytes"                     "Purkinje neurons"                 "Pyramidal cells"                  "Radial glia cells"               
# [29] "Retinal ganglion cells"           "Satellite glial cells"            "Schwann cells"                    "Serotonergic neurons"            
# [33] "Tanycytes"                        "Trigeminal neurons"   


# overlap with genes detected in our data
my.detected <- rownames(killi.brain.clean)
my.marker.list.v2 <- vector(mode = "list",length = length(panglao.cell.types))

for (i in 1:length(my.marker.list)) {
  my.marker.list.v2[[i]] <- intersect(my.marker.list[[i]],my.detected)
}
names(my.marker.list.v2) <- names(my.marker.list)

# check marker lists
lapply(my.marker.list.v2, length)


### Generate anno table
# initialize
anno           <- data.frame(cbind(rep(names(my.marker.list.v2)[1],length(my.marker.list.v2[[1]])), my.marker.list.v2[[1]], rep(2,length(my.marker.list.v2[[1]]))))
colnames(anno) <- c("Type", "Marker", "Weight")

for (i in 2:length(my.marker.list.v2)) {
  anno <- rbind(anno,
                cbind("Type" = rep(names(my.marker.list.v2)[i],length(my.marker.list.v2[[i]])),
                      "Marker" = my.marker.list.v2[[i]],
                      "Weight" = rep(2,length(my.marker.list.v2[[i]]))))
  
}

############ Section 2 - Pre-processing the data for scSorter ############

# get highly variable genes and filter out genes with non-zero expression in less than 5% of total cells.
topgenes        <- head(VariableFeatures(killi.brain.clean), 2500)
expr            <- GetAssayData(killi.brain.clean)
topgene_filter  <- rowSums(expr[topgenes, ]!=0) > ncol(expr)*.05
topgenes        <- topgenes[topgene_filter]

# At last, we subset the preprocessed expression data.
# Keep all genes that can be markers, even if not for cell types in that tissue
picked_genes    <- unique(c(anno$Marker, topgenes))
print(length(picked_genes)) # anno  ### 1861
expr.sctype     <- expr[rownames(expr) %in% picked_genes, ]

# Now, we are ready to run scSorter.


############ Section 3 - Running scSorter ############
# run scSorter
rts.brain  <- scSorter(expr.sctype , anno)
save(rts.brain, file = paste0(Sys.Date(),"_scSorter_out_Brain_PanglaoDB.RData"))

############ Section 4 - gate ScSorter calls back to Seurat  ############
# get calls from scSorter calls
killi.brain.clean@meta.data$ScSorter_PanglaoDB <- "other" # initialize
killi.brain.clean@meta.data[colnames(expr.sctype ), ]$ScSorter_PanglaoDB <- rts.brain$Pred_Type


pdf(paste(Sys.Date(),"Killifish_Brain_UMAP_color_by_ScSorter_PanglaoDB_call.pdf", sep = "_"), height = 5, width = 10)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "ScSorter_PanglaoDB", raster = F)
dev.off()

pdf(paste(Sys.Date(),"Killifish_Brain_UMAP_color_by_ScSorter_PanglaoDB_call_RASTER.pdf", sep = "_"), height = 5, width = 10)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "ScSorter_PanglaoDB", raster = T, raster.dpi = c(600,600))
dev.off()

# save annotated object
save(killi.brain.clean, file = paste0(Sys.Date(),"_Seurat_object_with_ScSorter_PanglaoDB_Raw_output.RData"))
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###



##############################################################################################################################
#####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Zebrafish markers from Raj et al 2018    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%########
##############################################################################################################################

# Load Zebrafish/Killifish homology table
# (best Zebrafish hit to killifish to get conversion)
zebrafish.homol <- read.csv("../zebrafish_alignment/2022-03-15_Zebrafish_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-5_Minimal_HOMOLOGY_TABLE_REV.txt", sep = "\t", header = T)


# ###%%%%%%%%%%%%%%%%%%%%%### A. Filter and prep PanglaoDB marker data  ###%%%%%%%%%%%%%%%###

############ Section 1 - Generate annotation object for scSorter ############
# Prepare the annotation for scSorter
# Read in marker file
raj.markers <- read_xlsx("../Atlas_Info/Zebrafish/41587_2018_BFnbt4103_MOESM56_ESM.xlsx", n_max = 63)

# Clean up cell type names
raj.markers$cell.type <- raj.markers$Identity
raj.markers$cell.type <- gsub("+","" ,raj.markers$cell.type, fixed = TRUE) # remove + signs
raj.markers$cell.type <- gsub("[[:space:]]*$", "",raj.markers$cell.type)   # remove trailing spaces
raj.markers$cell.type <- gsub("(","" ,raj.markers$cell.type, fixed = TRUE)
raj.markers$cell.type <- gsub(")","" ,raj.markers$cell.type, fixed = TRUE)
raj.markers$cell.type <- gsub(" ","_",raj.markers$cell.type, fixed = TRUE)
raj.markers$cell.type <- gsub("/","_",raj.markers$cell.type, fixed = TRUE)
raj.markers$cell.type <- gsub(",","_",raj.markers$cell.type, fixed = TRUE)
raj.markers$cell.type <- gsub("_+","_",raj.markers$cell.type)
raj.markers$cell.type <- gsub("__","_",raj.markers$cell.type, fixed = TRUE)

# harminize for the cell types spelled 2 different ways
raj.markers$cell.type[raj.markers$cell.type %in% "Granule_cell" ] <- "Granule_cells" 

# get cell types
raj.cell.types <- sort(setdiff(unique(raj.markers$cell.type),"-"))

# create marker list object to summarize all marker genes
my.marker.list         <- vector(mode = "list",length = length(raj.cell.types))
names(my.marker.list)  <- raj.cell.types

for (i in 1:length(raj.cell.types)) {
  
  # get unique markers from Raj et al 2018
  my.markers <-  unique(unlist(strsplit(raj.markers$Markers[raj.markers$cell.type %in% raj.cell.types[i]],", ")))
  
  # get closest Nfur corresponding gene (unique as well)
  my.marker.list[[i]] <- unique(zebrafish.homol$Nfur_Symbol[zebrafish.homol$DanRer_Symbol %in% my.markers])
  
}

names(my.marker.list)
# [1] "Cholinergic_cells"               "Cranial_ganglion"                "Diencephalon"                    "Differentiating_Oligodendrocyte"
# [5] "Dorsal_Habenula"                 "Endothelial_cells"               "Endothelial_cells_proliferative" "Ependymal_cells"                
# [9] "Erythrocytes"                    "Gaba_Hindbrain"                  "Granule_cells"                   "Granule_cells_Glut"             
# [13] "Hyp"                             "Hyp_Gaba"                        "Hyp_POA"                         "Microglia"                      
# [17] "Midbrain"                        "Midbrain_Glut"                   "Midbrain_Thalamus"               "Nascent_neurons"                
# [21] "Oligodendrocyte"                 "OPC"                             "Optic_Tectum"                    "Optic_Tectum_Gaba"              
# [25] "Optic_Tectum_Radial_Glia"        "Pallium"                         "Pallium_Glut"                    "Perivascular-FGP_cells"         
# [29] "Phox2_Hindbrain_cells"           "Pineal_gland"                    "Progenitor"                      "Purkinje_cells"                 
# [33] "Radial_Glia"                     "Radial_Glia_Bergmann_Glia"       "Statoacoustic_ganglion"          "Subpallium_Gaba"                
# [37] "Tcell_immune_cell"               "Telencephalon_Glut"              "Torus_Longitudinalis"            "URL_Progenitor"                 
# [41] "Ventral_Forebrain_Gaba_Glut"     "Ventral_Habenula"         

# overlap with genes detected in our data
my.detected <- rownames(killi.brain.clean)
my.marker.list.v2 <- vector(mode = "list",length = length(raj.cell.types))

for (i in 1:length(my.marker.list)) {
  my.marker.list.v2[[i]] <- intersect(my.marker.list[[i]],my.detected)
}
names(my.marker.list.v2) <- names(my.marker.list)

# check marker lists
lapply(my.marker.list.v2, length)


### Generate anno table
# initialize
anno           <- data.frame(cbind(rep(names(my.marker.list.v2)[1],length(my.marker.list.v2[[1]])), my.marker.list.v2[[1]], rep(2,length(my.marker.list.v2[[1]]))))
colnames(anno) <- c("Type", "Marker", "Weight")

for (i in 2:length(my.marker.list.v2)) {
  anno <- rbind(anno,
                cbind("Type" = rep(names(my.marker.list.v2)[i],length(my.marker.list.v2[[i]])),
                      "Marker" = my.marker.list.v2[[i]],
                      "Weight" = rep(2,length(my.marker.list.v2[[i]]))))
  
}

############ Section 2 - Pre-processing the data for scSorter ############

# get highly variable genes and filter out genes with non-zero expression in less than 5% of total cells.
topgenes        <- head(VariableFeatures(killi.brain.clean), 2500)
expr            <- GetAssayData(killi.brain.clean)
topgene_filter  <- rowSums(expr[topgenes, ]!=0) > ncol(expr)*.05
topgenes        <- topgenes[topgene_filter]

# At last, we subset the preprocessed expression data.
# Keep all genes that can be markers, even if not for cell types in that tissue
picked_genes    <- unique(c(anno$Marker, topgenes))
print(length(picked_genes)) # anno  ### 1146
expr.sctype     <- expr[rownames(expr) %in% picked_genes, ]

# Now, we are ready to run scSorter.

############ Section 3 - Running scSorter ############
# run scSorter
rts.brain.z  <- scSorter(expr.sctype , anno)
save(rts.brain.z, file = paste0(Sys.Date(),"_scSorter_out_Brain_Raj2018.RData"))

############ Section 4 - gate ScSorter calls back to Seurat  ############
# get calls from scSorter calls
killi.brain.clean@meta.data$ScSorter_Raj2018 <- "other" # initialize
killi.brain.clean@meta.data[colnames(expr.sctype ), ]$ScSorter_Raj2018 <- rts.brain.z$Pred_Type

pdf(paste(Sys.Date(),"Killifish_Brain_UMAP_color_by_ScSorter_Raj2018_call_RASTER.pdf", sep = "_"), height = 5, width = 15)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "ScSorter_Raj2018", raster = T, raster.dpi = c(600,600))
dev.off()

# save annotated object
save(killi.brain.clean, file = paste0(Sys.Date(),"_Seurat_object_with_ScSorter_PanglaoDB_Raj2018_Raw_output.RData"))
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###

###############################################################################################################################################################


###############################################################################################################################################################
#### 3. Try scMagic
# https://github.com/TianLab-Bioinfo/scMAGIC

# ###%%%%%%%%%%%%%%###  0. Before running scMAGIC, please firstly run the following codes  ###%%%%%%%%%%%%%%###
library(reticulate)
py_config()  # your python environment
print(py_module_available('numpy')) # whether the "numpy" has been installed
np      <- import("numpy")
np.exp2 <- np$exp2
np.max  <- np$max
np.sum  <- np$sum

# For Windows users, please set "method_HVGene = 'SciBet_R'"
#### Also for MacOSX !!!!
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###


# ###%%%%%%%%%%%%%%### A. Mousify gene expression matrix for comparison to mouse brain atlases  ###%%%%%%%%%%%%%%###
###### "mousify" our dataset

# Load Mouse/Killifish homology table
# (best mouse hit to killifish to get conversion)
mouse.homol    <- read.csv("../Mouse_alignment/2022-10-11_Mouse_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-3_Minimal_HOMOLOGY_TABLE_REV.txt", sep = "\t", header = T)
mouse.homol.cl <- unique(mouse.homol[,3:4])

# split Seurat by cohort to do memory save
brain.list     <- SplitObject(killi.brain.clean, split.by = "Batch")

# Create object list to receive mousified counts
mm.cts.list        <- vector(length = length(brain.list), mode = "list")
names(mm.cts.list) <- names(brain.list)

# loop over cohorts to get mousified counts
for (i in 1:length(brain.list)){
  # Extract killifish data counts
  killi.cts          <- as.data.frame(brain.list[[i]][["RNA"]]@counts)
  
  # merge killi/mouse based on BLAST results
  killi.cts.ann <- merge(killi.cts, unique(mouse.homol[,c("Nfur_Symbol","Mmu_Symbol")]), by.x = "row.names", by.y = "Nfur_Symbol")
  
  # summarize based on mouse (sum paralogs if needed)
  my.killi.cts.2           <- aggregate(killi.cts.ann[,!(colnames(killi.cts.ann) %in% c("Row.names","Mmu_Symbol"))], by = list(as.factor(killi.cts.ann$Mmu_Symbol)), FUN = 'sum')
  rownames(my.killi.cts.2) <- my.killi.cts.2$Group.1
  mm.cts.list[[i]]         <- my.killi.cts.2[,!(colnames(my.killi.cts.2) %in% "Group.1")]
  
}

# save mousified count matrix list
save(mm.cts.list, file = paste0(Sys.Date(),"_Mousified_Killifish_scCountMatrices_per_Cohort.RData"))

# remove lines with too many zeros from mousified data
for (i in 1:length(mm.cts.list)) {
  my.null <- which(rowSums(mm.cts.list[[i]]) < 100)
  mm.cts.list[[i]]         <- mm.cts.list[[i]][-my.null,]
}
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###




# ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%### B. Filter and prep Zeisel training data  ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
# ### Load training [https://github.com/pcahan1/singleCellNet#cs_train ; 2023-07-28]
# Downloaded data from SCNET webiste
load("../Atlas_Info/stList_Zeisel_SCNET.rda")
# stList_zeisel

table(stList_zeisel$sampTab$Class)
# Astrocytes      Ependymal         Immune        Neurons         Oligos PeripheralGlia       Vascular 
# 200            150            100           1050            100            144            200 

table(stList_zeisel$sampTab$newAnn3)
# Astrocytes                          Cerebellum neurons       Cholinergic and monoaminergic neurons                    Choroid epithelial cells 
# 50                                          50                                          50                                          50 
# Dentate gyrus granule neurons        Dentate gyrus radial glia-like cells    Di- and mesencephalon excitatory neurons    Di- and mesencephalon inhibitory neurons 
# 50                                          50                                          50                                          50 
# Enteric glia                             Enteric neurons                             Ependymal cells                   Glutamatergic neuroblasts 
# 50                                          50                                          50                                          50 
# Hindbrain neurons                                   Microglia               Non-glutamatergic neuroblasts                 Olfactory ensheathing cells 
# 50                                          50                                          50                                          50 
# Olfactory inhibitory neurons             Oligodendrocyte precursor cells                            Oligodendrocytes                         Peptidergic neurons 
# 50                                          50                                          50                                          50 
# Pericytes    Peripheral sensory neurofilament neurons  Peripheral sensory non-peptidergic neurons      Peripheral sensory peptidergic neurons 
# 50                                          50                                          50                                          50 
# Perivascular macrophages                              Satellite glia                               Schwann cells              Spinal cord excitatory neurons 
# 50                                          50                                          44                                          50 
# Spinal cord inhibitory neurons       Subcommissural organ hypendymal cells  Subventricular zone radial glia-like cells             Sympathetic cholinergic neurons 
# 50                                          50                                          50                                          50 
# Sympathetic noradrenergic neurons       Telencephalon inhibitory interneurons Telencephalon projecting excitatory neurons Telencephalon projecting inhibitory neurons 
# 50                                          50                                          50                                          50 
# Vascular and leptomeningeal cells                  Vascular endothelial cells                Vascular smooth muscle cells 
# 50                                          50                                          50 

# format reference dataset
ref.mtx     <- stList_zeisel$expDat
ref.labels  <- stList_zeisel$sampTab$newAnn3

# Create object list to receive scMAGIC results
killibrain.scMAGIC        <- vector(mode = "list", length = length(mm.cts.list))
names(killibrain.scMAGIC) <- names(mm.cts.list)

# Loop over cohorts to run scMAGIc
for (i in 1:length(mm.cts.list)){
  # run scMAGIC
  killibrain.scMAGIC[[i]] <- scMAGIC(exp_sc_mat         =   mm.cts.list[[i]],
                                     exp_ref_mat        =   ref.mtx         ,
                                     exp_ref_label      =   ref.labels      , 
                                     atlas              =   NULL            , 
                                     method_HVGene      =   'SciBet_R'      , 
                                     num_threads        =   1               ,
                                     method_findmarker  =   'Seurat'        ,
                                     cluster_num_pc     =   19              , # from PCA analysis
                                     min_cell           =   10              ,
                                     method1            =   'spearman'      ,
                                     percent_high_exp   =    0.8            )
}

save(killibrain.scMAGIC, file = paste0(Sys.Date(),"_scMAGIC_ZeiselRef_level3_RawOutput.RData"))

#######
# gate back to killi brain seurat object
killi.brain.clean@meta.data$scMAGIC_Zeisel <- NA

# loop over cohorts
for (i in 1:length(killibrain.scMAGIC)){
  
  # grab cell ids
  my.cell.ids <- rownames(killibrain.scMAGIC[[i]])
  my.cell.ids <- my.cell.ids[-grep('NA',my.cell.ids)] # remove NAs
  
  # populate meta-data
  killi.brain.clean@meta.data[my.cell.ids,]$scMAGIC_Zeisel <- as.vector(killibrain.scMAGIC[[i]][my.cell.ids,])
  
}

# Plot Raw output
pdf(paste(Sys.Date(),"Killifish_Brain_UMAP_color_by_scMAGIC_Zeisel_call_RASTER.pdf", sep = "_"), height = 5, width = 15)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "scMAGIC_Zeisel", raster = T, raster.dpi = c(600,600))
dev.off()

# Save object with raw annotations
save(killi.brain.clean, file = paste0(Sys.Date(),"_Seurat_object_with_scSorterPanglaoDB_Raj2018_scMAGICZeisel_Raw_output.RData"))
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###


# ###%%%%%%%%%%%%%%### C. Filter and prep Ximerakis/Rubin aging Brain training data  ###%%%%%%%%%%%%%%###
# ### Load reference data [https://singlecell.broadinstitute.org/single_cell/study/SCP263/aging-mouse-brain ; 2023-08-04]

#### expression matrix
ximerakis.expr           <- read.csv("../Atlas_Info/SCP263_Ximerakis_expression_Aging_mouse_brain_portal_data_updated.txt", header = T, sep = "\t")
rownames(ximerakis.expr) <- ximerakis.expr$GENE

ximerakis.expr <- ximerakis.expr[,-1] # remove gene column

#### meta data
ximerakis.meta     <- read.csv("../Atlas_Info/SCP263_Ximerakis_meta_Aging_mouse_brain_portal_data.txt"              , header = T, sep = "\t")
ximerakis.cellkey  <- read.csv("../Atlas_Info/SCP263_Ximerakis_cellType_abbreviation_key.txt"                       , header = F, sep = "\t")

# clean up comment line
ximerakis.meta <- ximerakis.meta[-1,]

table(ximerakis.meta$cell_type)
# ABC           ARP           ASC           CPC            DC            EC           EPC         Hb_VC        HypEPC           MAC            MG 
# 307           184          6747            89            55          2413           274            81            12           377          3910 
# MNC         NendC NEUR_immature   NEUR_mature          NEUT           NRP           NSC           OEG           OLG           OPC            PC 
# 77           394           162          5135            29            82           166           892         12384          2187           735 
# TNC          VLMC          VSMC 
# 29           105           243 

table(ximerakis.meta$cell_class)
# Astrocyte_lineage         Ependymal_cells            Immune_cells        Neuronal_lineage Oligodendrocyte_lineage       Vasculature_cells 
# 7097                     404                    4448                    5773                   15463                    3884 

# Use plyr to extract full names of cell types for ease of reading
ximerakis.meta$cell_type_full <- plyr::mapvalues(ximerakis.meta$cell_type, from = ximerakis.cellkey$V1, to = ximerakis.cellkey$V2)
View(ximerakis.meta)

# check correspondence
sum(ximerakis.meta$NAME == colnames(ximerakis.expr)) # 37069
length(ximerakis.meta$NAME) # 37069


# format reference dataset
ref.mtx.x     <- ximerakis.expr
ref.labels.x  <- ximerakis.meta$cell_type_full

# Create object list to receive scMAGIC results
killibrain.scMAGIC.x        <- vector(mode = "list", length = length(mm.cts.list))
names(killibrain.scMAGIC.x) <- names(mm.cts.list)

# Loop over cohorts to run scMAGIc
for (i in 1:length(mm.cts.list)){
  # run scMAGIC
  killibrain.scMAGIC.x[[i]] <- scMAGIC(exp_sc_mat         =   mm.cts.list[[i]],
                                       exp_ref_mat        =   ref.mtx.x         ,
                                       exp_ref_label      =   ref.labels.x      ,
                                       atlas              =   NULL            ,
                                       method_HVGene      =   'SciBet_R'      ,
                                       num_threads        =   1               ,
                                       method_findmarker  =   'Seurat'        ,
                                       cluster_num_pc     =   19              ,
                                       min_cell           =   10              ,
                                       method1            =   'spearman'      ,
                                       percent_high_exp   =    0.8            )
}

save(killibrain.scMAGIC.x, file = paste0(Sys.Date(),"_scMAGIC_Ximerakis_RawOutput.RData"))


#######
# gate back to killi brain seurat object
killi.brain.clean@meta.data$scMAGIC_Ximerakis <- NA

# loop over cohorts
for (i in 1:length(killibrain.scMAGIC.x)){
  
  # grab cell ids
  my.cell.ids <- rownames(killibrain.scMAGIC.x[[i]])
  my.cell.ids <- my.cell.ids[-grep('NA',my.cell.ids)] # remove NAs
  
  # populate meta-data
  killi.brain.clean@meta.data[my.cell.ids,]$scMAGIC_Ximerakis <- as.vector(killibrain.scMAGIC.x[[i]][my.cell.ids,])
  
}

# Plot Raw output
pdf(paste(Sys.Date(),"Killifish_Brain_UMAP_color_by_scMAGIC_Ximerakis_call_RASTER.pdf", sep = "_"), height = 5, width = 15)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "scMAGIC_Ximerakis", raster = T, raster.dpi = c(600,600))
dev.off()


# Save object with raw annotations
save(killi.brain.clean, file = paste0(Sys.Date(),"_Seurat_object_with_scSorterPanglaoDB_Raj2018_scMAGICZeiselXimerakis_Raw_output.RData"))
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###
################################################################################################################################################################


################################################################################################################################################################
##### 4. Calculate cluster markers and plot potential marker gene expression


############### RNA SNN 1.2
# find markers for every cluster compared to all remaining cells, report only the positive ones
killi.brain.clean <- SetIdent(object = killi.brain.clean, value = 'RNA_snn_res.1.2')

killi.markers_1.2 <- FindAllMarkers(killi.brain.clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

# get top 5 and top 10
killi.markers_1.2 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)  -> top5_1.2
killi.markers_1.2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10_1.2

# plot heatmap
my.killi.heat.2_1.2 <- DoHeatmap(killi.brain.clean, features = top5_1.2$gene, group.by = 'RNA_snn_res.1.2',   size = 3) + scale_fill_gradientn(colors = c("blue", "white", "red"))  + theme(axis.text.y = element_text(size = 1))

png(paste0(Sys.Date(),"Top5_marker_heatmap_RNA_RNA_res_1.2.png"), width = 40, height = 20, units = "cm", res = 300)
my.killi.heat.2_1.2
dev.off()

# Dotplot of top5 findmarkers markers by cluster
pdf(paste(Sys.Date(),"Killifish_Brain_Dotplot_TOP5_markers_by_RNA_snn_res.1.2_CLEAN.pdf", sep = "_"), height = 12, width = 18)
DotPlot(killi.brain.clean,
        features = unique(top5_1.2$gene),
        group.by = "RNA_snn_res.1.2") + RotatedAxis()
dev.off()


# write markers to file
write.table(data.frame(top5_1.2), file = paste(Sys.Date(),"Seurat_top5_markers_RNA_snn_res.1.2_Clustering.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(data.frame(top10_1.2), file = paste(Sys.Date(),"Seurat_top10_markers_RNA_snn_res.1.2_Clustering.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(killi.markers_1.2[killi.markers_1.2$p_val_adj < 0.05,], file = paste(Sys.Date(),"Seurat_ALL_FDR5_markers_RNA_snn_res.1.2_Clustering.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
save(killi.markers_1.2, file = paste(Sys.Date(),"Seurat_markers_RNA_snn_res.1.2.RData", sep = "_"))



##### Annotate Markers
# read in killifish gene annotation table
killi.annot     <- read.csv('../Mouse_alignment/GCF_001465895.1_Nfu_20140520_feature_table.txt', header = T, sep = "\t")
killi.annot.flt <- unique(killi.annot[,c("symbol", "name")])

# make function to aggregate annotation
single_strg <- function (my.vec) {
  out <- paste(my.vec, collapse = ", ")
}

nfur.annot               <- aggregate(data.frame(killi.annot.flt[,2]), by = list("Nfur_gene" = killi.annot.flt$symbol), FUN = "single_strg")
colnames(nfur.annot)[2]  <- "Description"
nfur.annot$Description   <- gsub("^, ","", nfur.annot$Description)

# merge description
annot.killi.markers_1.2 <- merge(killi.markers_1.2, nfur.annot, by.x = "gene", by.y = "Nfur_gene")

# reorder
annot.killi.markers_1.2.v2 <- annot.killi.markers_1.2[,c("cluster","avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj", "gene","Description" )]
annot.killi.markers_1.2.v2 <- annot.killi.markers_1.2.v2[order(annot.killi.markers_1.2.v2$cluster, decreasing = F),]

# get top 5 and top 10
annot.killi.markers_1.2.v2 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)  -> ann.top5_1.2
annot.killi.markers_1.2.v2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> ann.top10_1.2

# write annotated markers to file
write.table( data.frame(ann.top5_1.2 ), file = paste(Sys.Date(),"Seurat_top5_markers_RNA_snn_res.1.2_Clustering_ANNOT.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table( data.frame(ann.top10_1.2), file = paste(Sys.Date(),"Seurat_top10_markers_RNA_snn_res.1.2_Clustering_ANNOT.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(annot.killi.markers_1.2.v2[annot.killi.markers_1.2.v2$p_val_adj < 0.05,], file = paste(Sys.Date(),"Seurat_ALL_FDR5_markers_RNA_snn_res.1.2_Clustering_ANNOT.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
save(annot.killi.markers_1.2.v2, file = paste(Sys.Date(),"Seurat_markers_RNA_snn_res.1.2_ANNOT.RData", sep = "_"))
###############################################################################################################################################################



################################################################################################################################################################
##### 5. Parse cell type predictions to assign an identity to each cell cluster


################################################################################
# Parse results from each prediction by SNN cluster group (to apply majority and decrease noise)
# we are using the 1.2 clustering resolution to improve granularity of assignment
my.cluster.annot                 <- data.frame(matrix(0, length(unique(killi.brain.clean@meta.data$RNA_snn_res.1.2)), 17) )
rownames(my.cluster.annot)       <- paste0("Cluster_",0:(length(unique(killi.brain.clean@meta.data$RNA_snn_res.1.2))-1))
colnames(my.cluster.annot)       <- c("RNA_snn_res.1.2",
                                      "Top_ScSorter_PanglaoDB"    , "Top_ScSorter_Perc_PanglaoDB"    ,
                                      "Second_ScSorter_PanglaoDB" , "Second_ScSorter_Perc_PanglaoDB" ,
                                      "Top_ScSorter_Raj2018"      , "Top_ScSorter_Perc_Raj2018"      ,
                                      "Second_ScSorter_Raj2018"   , "Second_ScSorter_Perc_Raj2018"   ,
                                      "Top_scMAGICZeisel"         , "Top_scMAGICZeisel_Perc"         ,
                                      "Second_scMAGICZeisel"      , "Second_scMAGICZeisel_Perc"      ,
                                      "Top_scMAGICXimerakis"      , "Top_scMAGICXimerakis_Perc"      ,
                                      "Second_scMAGICXimerakis"   , "Second_scMAGICXimerakis_Perc"   )

my.cluster.annot$RNA_snn_res.1.2 <- 0:(length(unique(killi.brain.clean@meta.data$RNA_snn_res.1.2))-1) # initialize


# function to parse predictions
get_top2_preds <- function (clust.preds, my.colname){
  
  # tabulate by cell type called by algorithm
  tab.res <- table(meta.snn.sub[,my.colname])
  
  # Sort to get top 2
  tab.res.sort <- sort(tab.res, decreasing = T)
  
  # parse information of top 2 most frequent predictions in the cluster
  a  <- names(tab.res.sort)[1]
  b  <- 100*round(tab.res.sort[1]/sum(tab.res.sort), digits = 4)
  
  c  <- names(tab.res.sort)[2]
  d  <- 100*round(tab.res.sort[2]/sum(tab.res.sort), digits = 4)
  
  return(c(a,b,c,d))
}

####
for (i in 1:nrow(my.cluster.annot)) {
  
  my.snn.clust <- my.cluster.annot$RNA_snn_res.1.2[i]
  
  # subset cluster from metadata dataframe to extract info
  meta.snn.sub <- killi.brain.clean@meta.data[killi.brain.clean@meta.data$RNA_snn_res.1.2 == my.snn.clust,]
  
  # parse information for ScSorter PanglaoDB
  scSorterPang.parse <- get_top2_preds(meta.snn.sub, "ScSorter_PanglaoDB")
  my.cluster.annot[i,]$Top_ScSorter_PanglaoDB              <- scSorterPang.parse[1]
  my.cluster.annot[i,]$Top_ScSorter_Perc_PanglaoDB         <- scSorterPang.parse[2]
  my.cluster.annot[i,]$Second_ScSorter_PanglaoDB           <- scSorterPang.parse[3]
  my.cluster.annot[i,]$Second_ScSorter_Perc_PanglaoDB      <- scSorterPang.parse[4]
  
  # parse information for ScSorter Raj2018
  scSorterRaj.parse <- get_top2_preds(meta.snn.sub, "ScSorter_Raj2018")
  my.cluster.annot[i,]$Top_ScSorter_Raj2018              <- scSorterRaj.parse[1]
  my.cluster.annot[i,]$Top_ScSorter_Perc_Raj2018         <- scSorterRaj.parse[2]
  my.cluster.annot[i,]$Second_ScSorter_Raj2018           <- scSorterRaj.parse[3]
  my.cluster.annot[i,]$Second_ScSorter_Perc_Raj2018      <- scSorterRaj.parse[4]
  
  # parse information for scMAGICZeisel
  scMAGICZeisel.parse <- get_top2_preds(meta.snn.sub, "scMAGIC_Zeisel")
  my.cluster.annot[i,]$Top_scMAGICZeisel                <- scMAGICZeisel.parse[1]
  my.cluster.annot[i,]$Top_scMAGICZeisel_Perc           <- scMAGICZeisel.parse[2]
  my.cluster.annot[i,]$Second_scMAGICZeisel             <- scMAGICZeisel.parse[3]
  my.cluster.annot[i,]$Second_scMAGICZeisel_Perc        <- scMAGICZeisel.parse[4]
  
  # parse information for scMAGICXimerakis
  scMAGICXimerakis.parse <- get_top2_preds(meta.snn.sub, "scMAGIC_Ximerakis")
  my.cluster.annot[i,]$Top_scMAGICXimerakis            <- scMAGICXimerakis.parse[1]
  my.cluster.annot[i,]$Top_scMAGICXimerakis_Perc       <- scMAGICXimerakis.parse[2]
  my.cluster.annot[i,]$Second_scMAGICXimerakis         <- scMAGICXimerakis.parse[3]
  my.cluster.annot[i,]$Second_scMAGICXimerakis_Perc    <- scMAGICXimerakis.parse[4]
  
}

write.table(my.cluster.annot, file = paste(Sys.Date(),"Parsed_Cell_Annotation_Results_by_RNA_snn_res.1.2_Clusters.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = T)

##### add additional marker plotting to help annotate
pdf(paste(Sys.Date(),"Killifish_Brain_Dotplot_PROLIFERATION_markers_by_RNA_snn_res.1.2.pdf", sep = "_"), height = 8, width = 3.5)
DotPlot(killi.brain.clean,
        features = c("pcna","mki67"),                   # Proliferation
        group.by = "RNA_snn_res.1.2") + RotatedAxis()
dev.off()
# 16,17, 25, 31, 32 have proliferation markers

# res 1.2
pdf(paste(Sys.Date(),"Killifish_Brain_Dotplot_Known_Cell_type_markers_by_RNA_snn_res.1.2_REVISED.pdf", sep = "_"), height = 8, width = 14)
DotPlot(killi.brain.clean,
        features = c("olig1","olig2","mpz","LOC107394899",                    # Oligodendrocyte/OPC                      ### LOC107394899 / mbpa myelin basic protein a; 
                     "marco","csf1r", "ptprc",                                # mph/microglia
                     "s100b", "slc1a2",                                       # astrocyte/radial glia
                     "rbfox3", "map2", "ncam1",                               # mature neuron
                     "fat2", "neurod1", "LOC107376653" ,                      # Granule excitatory neuron / glutamatergic ### LOC107376653 (Slc17a7)
                     "pvalb",                                                 # pvalb interneurons
                     "lhx1","ca8","itpr1",                                    # Purkinje cells
                     "gad1", "gad2",  "LOC107384443", "LOC107391088",         # GABAergic neurons                         ### LOC107391088 (scl6a1), LOC107384443 (gad1l)
                     "LOC107386535",  "clu",                                  # ependymal cells                           ### LOC107386535 (epd)
                     "sox5",  "sox2","fgfr2","gli3", "LOC107385497", "efna2", # NSPCs                                     ### gli3, fgfr2, efna2, LOC107385497 (msi1)
                     "LOC107378372", "LOC107378374",                          # Erythrocytes (hemoglobin)
                     "tpm1", "acta2"                                          # Smooth muscle
        ),
        group.by = "RNA_snn_res.1.2") + RotatedAxis()
dev.off()


# gli3
# LOC107385497	RNA-binding protein Musashi homolog 1
# fgfr2	fibroblast growth factor receptor 2, transcript variant X2, fibroblast growth factor receptor 2, transcript variant X3, fibroblast growth factor receptor 2, transcript variant X1, fibroblast growth factor receptor 2, transcript variant X4, fibroblast growth factor receptor 2 isoform X2, fibroblast growth factor receptor 2 isoform X1, fibroblast growth factor receptor 2 isoform X3, fibroblast growth factor receptor 2 isoform X4
# efna2	ephrin-A2

pdf(paste(Sys.Date(),"Killifish_Brain_Dotplot_Purkinje_markers_by_RNA_snn_res.1.2.pdf", sep = "_"), height = 8, width = 5)
DotPlot(killi.brain.clean,
        features = c("pvalb","lhx1","LOC107384050","ca8","itpr1","aldoc","LOC107388709"),
        group.by = "RNA_snn_res.1.2") + RotatedAxis()
dev.off()


################################################################################################################################################################

################################################################################################################################################################
##### 6. Assign cell type predictions to each cell cluster

# Used the predictions combined with marker gene expression to annotate in Excel

# Import predictions
my.annot <- read_xlsx("2023-08-21_Parsed_Cell_Annotation_Results_by_RNA_snn_res.1.2_Clusters_v2.xlsx")

# get calls from scSorter calls
killi.brain.clean@meta.data$Cell_Identity <- NA # initialize

for (i in 0:35) {
  killi.brain.clean@meta.data[killi.brain.clean@meta.data$RNA_snn_res.1.2 == i,]$Cell_Identity <- my.annot$`Final Call`[my.annot$`RNA snn res.1.2` == i]
}


table(killi.brain.clean@meta.data$Cell_Identity)
#       Astrocytes_Radial_Glia              Ependymal_cells                 Erythrocytes            GABAergic_neurons   Granule_Excitatory_Neurons 
#                        11270                         3321                         1745                         8418                        38801 
#                    Microglia               Neurons_misc_1               Neurons_misc_2               Neurons_misc_3               Neurons_misc_4 
#                         3466                          931                        14208                         6240                        86962 
#                        NSPCs             Oligodendrocytes                         OPCs               Purkinje_cells              PV_interneurons 
#                         4868                        14271                         4093                          334                         9158 
# Vascular_smooth_muscle_cells 
#                         1853 
length(unique(killi.brain.clean@meta.data$Cell_Identity)) # 16

# create your own color palette based on `seedcolors`
set.seed(123123) # stabilize
P16 = createPalette(16+3,  c("#ff0000", "#00ff00", "#0000ff"))
swatch(P16)

write.table(cbind("col" = P16[-c(1:3)], "cell_type" = sort(unique(killi.brain.clean@meta.data$Cell_Identity))), file = paste0(Sys.Date(),"_color_palette_annotation.txt"), sep = "\t", row.names = F)

# Plot Raw output
pdf(paste(Sys.Date(),"Killifish_Brain_UMAP_color_by_deNovo_annotation_RASTER.pdf", sep = "_"), height = 4, width = 7)
DimPlot(killi.brain.clean, reduction = "umap", group.by = "Cell_Identity", raster = T,
        raster.dpi = c(350,350), cols = as.vector(P16[-c(1:3)]))
dev.off()

# Save object with raw annotations
save(killi.brain.clean, file = paste0(Sys.Date(),"_Seurat_object_with_Manual_annotation.RData"))
################################################################################################################################################################


################################################################################################################################################################
##### 7. Marker plotting and QC

# Load object
load('2023-08-23_Seurat_object_with_Manual_annotation.RData')

# recalculate marker genes on the manual annotations
killi.brain.clean <- SetIdent(object = killi.brain.clean, value = 'Cell_Identity')

killi.markers.celltype <- FindAllMarkers(killi.brain.clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

# get top 5 and top 10
killi.markers.celltype %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)  -> top5_annot
killi.markers.celltype %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10_annot

# make function to aggregate annotation (keep only 1st description)
keep1 <- function (my.vec) {
  null.desc <- my.vec == "" # to keep the first non null string
  out <- my.vec[!null.desc][1]
}

nfur.annot.2               <- aggregate(data.frame(killi.annot.flt[,2]), by = list("Nfur_gene" = killi.annot.flt$symbol), FUN = "keep1")
colnames(nfur.annot.2)[2]  <- "Description"

# merge description
annot.killi.markers_celltype <- merge(killi.markers.celltype, nfur.annot.2, by.x = "gene", by.y = "Nfur_gene")

# reorder
annot.killi.markers_celltype.v2 <- annot.killi.markers_celltype[,c("cluster","avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj", "gene","Description" )]
annot.killi.markers_celltype.v2 <- annot.killi.markers_celltype.v2[order(annot.killi.markers_celltype.v2$cluster, decreasing = F),]

# get top 5 and top 10
annot.killi.markers_celltype.v2 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)  -> ann.top5_celltype
annot.killi.markers_celltype.v2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> ann.top10_celltype

# write annotated markers to file
write.table( data.frame(ann.top5_celltype ), file = paste(Sys.Date(),"Seurat_top5_markers_RNA_snn_res.celltype_Clustering_ANNOT.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table( data.frame(ann.top10_celltype), file = paste(Sys.Date(),"Seurat_top10_markers_RNA_snn_res.celltype_Clustering_ANNOT.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(annot.killi.markers_celltype.v2[annot.killi.markers_celltype.v2$p_val_adj < 0.05,], file = paste(Sys.Date(),"Seurat_ALL_FDR5_markers_RNA_snn_res.celltype_Clustering_ANNOT.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)

# reorder cell types for ease of plotting (non neuronal/neuronal)
killi.brain.clean@meta.data$Cell_Identity <- factor(x = killi.brain.clean@meta.data$Cell_Identity, 
                                                    levels = c("Astrocytes_Radial_Glia"      , 
                                                               "Ependymal_cells"             ,
                                                               "Erythrocytes"                ,
                                                               "Microglia"                   ,
                                                               "NSPCs"                       ,
                                                               "Oligodendrocytes"            ,
                                                               "OPCs"                        ,
                                                               "Vascular_smooth_muscle_cells",
                                                               "GABAergic_neurons"           ,
                                                               "Granule_Excitatory_Neurons"  ,       
                                                               "Neurons_misc_1"              ,
                                                               "Neurons_misc_2"              ,
                                                               "Neurons_misc_3"              ,
                                                               "Neurons_misc_4"              ,       
                                                               "Purkinje_cells"              ,
                                                               "PV_interneurons"              ))

###### Plot marker genes - mix of previously known and from the marker calculation above
pdf(paste(Sys.Date(),"Killifish_Brain_Dotplot_KnownAndPredicted_Cell_type_markers_by_ManualAnnotation.pdf", sep = "_"), height = 8, width = 16)
DotPlot(killi.brain.clean,
        features = c("s100b", "slc1a2","kcnj10",                              # astrocyte/radial glia
                     "LOC107386535", "LOC107383970","clu",                    # ependymal cells               ### LOC107386535 (ependymin-2-like), LOC107383970	(serotransferrin-like)
                     "LOC107378372", "LOC107378374","LOC107378381",           # Erythrocytes (hemoglobin)
                     "LOC107387973", "LOC107379395", "csf1r",                 # mph/microglia                 ### LOC107387973 (itgam/cd11b), LOC107379395	(apoeb)
                     "gli3","efna2", "LOC107385497",                          # NSPCs                         ### gli3, efna2, LOC107385497 (msi1)
                     "mpz","LOC107394899", "LOC107386530",                    # Oligodendrocyte               ### LOC107394899 (mbpa), LOC107386530 (plp)
                     "olig1","olig2", "sema5a","sox5",                        # OPCs
                     "tpm1", "LOC107375895", "ahnak",                         # Vascular smooth muscle        ### LOC107375895 (kcnq5)
                     
                     "LOC107384443","gad2", "LOC107391088",                   # GABAergic neurons             ### LOC107384443 (gad1l), LOC107391088 (scl6a1; GABA transporter)
                     "fat2", "LOC107392205", "LOC107376653",                  # Granule excitatory neuron     ### LOC107392205 (qka), LOC107376653 (Slc17a7/VGLUT1)
                     "LOC107392873","LOC107397057","LOC107387921",            # Neurons_misc_1                ### LOC107392873 (fcho1),LOC107397057 (gjc1),LOC107387921 (ngfr)
                     "LOC107385736","LOC107385201", "cdh24",                  # Neurons_misc_2                ### LOC107385736 (grm2), LOC107385201 (camkv)
                     "elavl4", "tmem163","LOC107380370",                      # Neurons_misc_3                ### LOC107380370 (elnl)
                     "nrp1","meis2","LOC107390816",                           # Neurons_misc_4                ### LOC107375999 (ebf3) / LOC107390816 cpne5a
                     "itpr1","mpped1","LOC107379976",                         # Purkinje cells                ### LOC107379976 (aldoca)
                     "pvalb", "pde2a",                                        # PV interneurons
                     "rbfox3", "map2", "ncam1"                                # mature neuron markers
        ),
        group.by = "Cell_Identity") + RotatedAxis()
dev.off()

#### for heatmap, needs more balance: subset 200 cells per cell type to plot (otherwise, readability is terrible)
# https://satijalab.org/seurat/v3.0/multimodal_vignette.html
killi.small.plot <- subset(killi.brain.clean, downsample = 500)

# Get normalized RNA values
killi.small.plot <- NormalizeData(killi.small.plot)
killi.small.plot <- FindVariableFeatures(killi.small.plot, selection.method = "vst", nfeatures = 5000)
killi.small.plot <- ScaleData(killi.small.plot, features = rownames(killi.small.plot))

pdf(paste(Sys.Date(),"Killifish_Brain_Heatmap_KnownAndPredicted_Cell_type_markers_by_ManualAnnotation_500cell_DS.pdf", sep = "_"), width = 17, height = 11)
DoHeatmap(killi.small.plot, features = c("s100b", "slc1a2","kcnj10",                              # astrocyte/radial glia
                                         "LOC107386535", "LOC107383970","clu",                    # ependymal cells               ### LOC107386535 (ependymin-2-like), LOC107383970	(serotransferrin-like)
                                         "LOC107378372", "LOC107378374","LOC107378381",           # Erythrocytes (hemoglobin)
                                         "LOC107387973", "LOC107379395", "csf1r",                 # mph/microglia                 ### LOC107387973 (itgam/cd11b), LOC107379395	(apoeb)
                                         "gli3","efna2", "LOC107385497",                          # NSPCs                         ### gli3, efna2, LOC107385497 (msi1)
                                         "mpz","LOC107394899", "LOC107386530",                    # Oligodendrocyte               ### LOC107394899 (mbpa), LOC107386530 (plp)
                                         "olig1","olig2", "sema5a","sox5",                        # OPCs
                                         "tpm1", "LOC107375895", "ahnak",                         # Vascular smooth muscle        ### LOC107375895 (kcnq5)
                                         
                                         "LOC107384443","gad2", "LOC107391088",                   # GABAergic neurons             ### LOC107384443 (gad1l), LOC107391088 (scl6a1; GABA transporter)
                                         "fat2", "LOC107392205", "LOC107376653",                  # Granule excitatory neuron     ### LOC107392205 (qka), LOC107376653 (Slc17a7/VGLUT1)
                                         "LOC107392873","LOC107397057","LOC107387921",            # Neurons_misc_1                ### LOC107392873 (fcho1),LOC107397057 (gjc1),LOC107387921 (ngfr)
                                         "LOC107385736","LOC107385201", "cdh24",                  # Neurons_misc_2                ### LOC107385736 (grm2), LOC107385201 (camkv)
                                         "elavl4", "tmem163","LOC107380370",                      # Neurons_misc_3                ### LOC107380370 (elnl)
                                         "nrp1","meis2","LOC107390816",                           # Neurons_misc_4                ### LOC107375999 (ebf3)
                                         "itpr1","mpped1","LOC107379976",                         # Purkinje cells                ### LOC107379976 (aldoca)
                                         "pvalb", "pde2a",                                        # PV interneurons
                                         "rbfox3", "map2", "ncam1"                                # mature neuron markers
), group.by = "Cell_Identity") + theme(axis.text.y = element_text(size = 12)) + scale_fill_gradientn(colors = c("darkblue", "white", "red"))
dev.off()

# Save object with FINAL annotations
save(killi.brain.clean, file = paste0(Sys.Date(),"_Seurat_object_with_Manual_annotation_FINAL.RData"))
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###


###############################################################################################################################################################
sink(file = paste(Sys.Date(),"_Killi_Brain_Atlas_ANNOT_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()


