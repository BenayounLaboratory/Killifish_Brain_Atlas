setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/snATAC_Brain_Aging_Meta/Subsetted_Bams')
options(stringsAsFactors = F)

#### Packages
library(Seurat)        # single cell general package
library(Signac)        # scATAC processing
library(bitops)        # 

# 2024-04-24
# Process scATAC to get lists of barcodes per cell type/per replicate
# to generate sample specific cell type specific bam files
# https://github.com/10XGenomics/subset-bam


#########################################################################################
# 0. Load annotated object
load("../Signac/2023-08-23_Killi_brain_ATAC_2Cohorts_FiltNorm_Singlets_ANNOTATEDfromRNA.RData")

brain.atac.singlets
# An object of class Seurat 
# 148827 features across 27171 samples within 2 assays 
# Active assay: RNA (22135 features, 0 variable features)
# 1 other assay present: ATAC
# 2 dimensional reductions calculated: lsi, umap

# extract barcode information
brain.atac.singlets@meta.data$barcode <- unlist(lapply(strsplit(rownames(brain.atac.singlets@meta.data),"_"),'[[',2))
brain.atac.singlets@meta.data <- brain.atac.singlets@meta.data[,-grep("prediction.score",colnames(brain.atac.singlets@meta.data))]

# extract meta data data frame
meta.data <- brain.atac.singlets@meta.data

# get sample names
my.samples   <- sort(unique(brain.atac.singlets@meta.data$sample))
# [1] "OF1" "OF2" "OM1" "OM2" "YF1" "YF2" "YM1" "YM2"

# get QC cell types
my.celltypes <- sort(unique(brain.atac.singlets@meta.data$predicted.id))
# [1] "Astrocytes_Radial_Glia"       "Ependymal_cells"              "Erythrocytes"                 "GABAergic_neurons"            "Granule_Excitatory_Neurons"  
# [6] "Microglia"                    "Neurons_misc_1"               "Neurons_misc_2"               "Neurons_misc_3"               "Neurons_misc_4"              
# [11] "NSPCs"                        "Oligodendrocytes"             "OPCs"                         "Purkinje_cells"               "PV_interneurons"             
# [16] "Vascular_smooth_muscle_cells"

# nb. of cells per cluster-sample
cell.per.samp.tab <- t(table(brain.atac.singlets@meta.data$predicted.id, brain.atac.singlets@meta.data$sample))

# cell types with at least 10 cells from every each sex/cohort sample
celltype.qc <- colnames(cell.per.samp.tab)[apply(cell.per.samp.tab >= 10, 2, sum) == 8]

# REMOVE "MISC" CELL TYPES
celltype.qc <- celltype.qc[-grep("misc",celltype.qc)]
celltype.qc
# keep:
# [1] "Astrocytes_Radial_Glia"     "Ependymal_cells"            "GABAergic_neurons"          "Granule_Excitatory_Neurons" "Microglia"                 
# [6] "NSPCs"                      "Oligodendrocytes"           "OPCs"                       "PV_interneurons"    

for (i in 1:length(my.samples)) {
  
  cur.samp <- my.samples[i]
  
  for (j in 1:length(celltype.qc)) {
    cur.celltype <- celltype.qc[j]
    
    meta.data.cur <- meta.data[bitAnd(meta.data$sample == cur.samp, meta.data$predicted.id == cur.celltype)>0,]
    
    write.table(meta.data.cur$barcode, file = paste(Sys.Date(),cur.celltype,cur.samp,"barcode_list.txt", sep = "_"), quote = F, col.names = F, row.names = F)
  }
  
}


###################################################################################################################


#######################
sink(file = paste(Sys.Date(),"_Killi_brain_ATAC_extract_barcodes_session_Info.txt", sep =""))
sessionInfo()
sink()
