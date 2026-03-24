setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Preprocessing/DecontX/Cohort_1/')
options(stringsAsFactors = F)
options (future.globals.maxSize = 32000 * 1024^2)

# General use packages
library('Seurat')
library(bitops)
library(sctransform)

# removal of ambient RNA with DecontX
# https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html#seurat
library('celda')

# Doublet identification packages
# https://github.com/chris-mcginnis-ucsf/DoubletFinder
library(DoubletFinder)

# cxds_bcds_hybrid related packages
# https://www.bioconductor.org/packages/release/bioc/vignettes/scds/inst/doc/scds.html
library(scds)
library(scater)
library(bitops)


################################################################################################
# 2022-12-23
# Start homogeneous processing of all snRNAseq brain dataset
#      - use soupx to diminish the impact of ambient RNA (incompatible with phantom purge)
#      - clean doublets, poor QC cells, etc in each library
#      - use intersection of "generous" run of scds and doubletfinder (~15% to account for 10x doublets and doublets from isolation)
# for groups: use Y/M/O/G instead of weeks, since even though we initially designed to use 5/15w, Ari used 6/16w
#
# 2023-01-05
# After running up to integration
# Cell types look jumbled.... May need more aggressive doublet filtering
#
# 2023-06-30
# marker genes aren't super tight to clusters across cohorts and after integration
# try alternate method to clean up ambient RNA, DecontX, to see if improvements are observed
#
# 2023-07-02
# Increasing pre-filtering stringency (ncoutRNA, etc)
# make predicted doublet number 1.15 fold what 10x says, since we know nuclei are more susceptible
#
# 2023-07-03
# Add decontX filter to remove any cell with predicted contamination level ≥ 0.25
################################################################################################

################################################################################################
#### 0. Assumed doublet information/to calculate %age for prediction
# Targeted Cell Recovery  # of Cells Loaded	Barcodes Detected	Singlets Multiplets	Multiplet Rate
# 3,000	                         4,950           ~3,000        ~2,900	   ~80	    ~2.4%
# 4,000	                         6,600           ~3,900	       ~3,800	   ~140     ~3.2%
# 5,000	                         8,250           ~4,800        ~4,600	   ~210     ~4.0%
# 6,000                          9,900           ~5,700	       ~5,400	   ~300	    ~4.8%
# 7,000                          11,550	         ~6,600	       ~6,200	   ~400	    ~5.6%
# 8,000                          13,200	         ~7,500	       ~7,000	   ~510	    ~6.4%
# 9,000                          14,850	         ~8,400	       ~7,700	   ~640	    ~7.2%
#10,000                          16,500	         ~9,200	       ~8,400	   ~780	    ~8.0%
#12,000                          19,800	         ~10,900	     ~9,800	   ~1,100   ~9.6%

pred.10x.dblt <- data.frame( "cell_number" = c(3000,4000,5000,6000,7000,8000,9000, 10000, 12000),
                             "dblt_rate"   = c(2.4 ,3.2 ,4.0 ,4.8 ,5.6 ,6.4 ,7.2 , 8.0  , 9.6))

pred_dblt_lm <- lm(dblt_rate ~ cell_number, data = pred.10x.dblt)

pdf(paste0(Sys.Date(), "_10x_cel_number_vs_doublet_rate.pdf"))
plot(dblt_rate ~ cell_number, data = pred.10x.dblt)
abline(pred_dblt_lm, col = "red", lty = "dashed")
dev.off()

# predict(pred_dblt_lm, data.frame("cell_number" = 8500))
################################################################################################


################################################################################################
#### 1a. Read data from CellRanger, perform background removal

# Calculate and clean the contribution of ambient RNA with DecontX
# read 10x libraries cell ranger gene barcode matrices for DecontX
cts.ZMZ_Y_F_1     <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-05-25_snRNAseq_brain_aging_set_1/CellRanger/ZMZ_5w_F_1_SIGAG5_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.ZMZ_Y_M_1     <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-05-25_snRNAseq_brain_aging_set_1/CellRanger/ZMZ_5w_M_1_SIGAG1_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.ZMZ_M_F_1     <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-05-25_snRNAseq_brain_aging_set_1/CellRanger/ZMZ_10w_F_1_SIGAG6_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.ZMZ_M_M_1     <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-05-25_snRNAseq_brain_aging_set_1/CellRanger/ZMZ_10w_M_1_SIGAG2_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.ZMZ_O_F_1     <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-05-25_snRNAseq_brain_aging_set_1/CellRanger/ZMZ_15w_F_1_SIGAG7_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.ZMZ_O_M_1     <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-05-25_snRNAseq_brain_aging_set_1/CellRanger/ZMZ_15w_M_1_SIGAG3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.ZMZ_G_F_1     <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-05-25_snRNAseq_brain_aging_set_1/CellRanger/ZMZ_26w_F_1_SIGAG8_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.ZMZ_G_M_1     <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-05-25_snRNAseq_brain_aging_set_1/CellRanger/ZMZ_26w_M_1_SIGAG4_FishTEDB_NR/outs/filtered_feature_bc_matrix/')

cts_raw.ZMZ_Y_F_1 <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-05-25_snRNAseq_brain_aging_set_1/CellRanger/ZMZ_5w_F_1_SIGAG5_FishTEDB_NR/outs/raw_feature_bc_matrix/')
cts_raw.ZMZ_Y_M_1 <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-05-25_snRNAseq_brain_aging_set_1/CellRanger/ZMZ_5w_M_1_SIGAG1_FishTEDB_NR/outs/raw_feature_bc_matrix/')
cts_raw.ZMZ_M_F_1 <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-05-25_snRNAseq_brain_aging_set_1/CellRanger/ZMZ_10w_F_1_SIGAG6_FishTEDB_NR/outs/raw_feature_bc_matrix/')
cts_raw.ZMZ_M_M_1 <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-05-25_snRNAseq_brain_aging_set_1/CellRanger/ZMZ_10w_M_1_SIGAG2_FishTEDB_NR/outs/raw_feature_bc_matrix/')
cts_raw.ZMZ_O_F_1 <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-05-25_snRNAseq_brain_aging_set_1/CellRanger/ZMZ_15w_F_1_SIGAG7_FishTEDB_NR/outs/raw_feature_bc_matrix/')
cts_raw.ZMZ_O_M_1 <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-05-25_snRNAseq_brain_aging_set_1/CellRanger/ZMZ_15w_M_1_SIGAG3_FishTEDB_NR/outs/raw_feature_bc_matrix/')
cts_raw.ZMZ_G_F_1 <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-05-25_snRNAseq_brain_aging_set_1/CellRanger/ZMZ_26w_F_1_SIGAG8_FishTEDB_NR/outs/raw_feature_bc_matrix/')
cts_raw.ZMZ_G_M_1 <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-05-25_snRNAseq_brain_aging_set_1/CellRanger/ZMZ_26w_M_1_SIGAG4_FishTEDB_NR/outs/raw_feature_bc_matrix/')

# Create SingleCellExperiment objects
sce.ZMZ_Y_F_1     <- SingleCellExperiment(list(counts = cts.ZMZ_Y_F_1    ))
sce.ZMZ_Y_M_1     <- SingleCellExperiment(list(counts = cts.ZMZ_Y_M_1    ))
sce.ZMZ_M_F_1     <- SingleCellExperiment(list(counts = cts.ZMZ_M_F_1    ))
sce.ZMZ_M_M_1     <- SingleCellExperiment(list(counts = cts.ZMZ_M_M_1    ))
sce.ZMZ_O_F_1     <- SingleCellExperiment(list(counts = cts.ZMZ_O_F_1    ))
sce.ZMZ_O_M_1     <- SingleCellExperiment(list(counts = cts.ZMZ_O_M_1    ))
sce.ZMZ_G_F_1     <- SingleCellExperiment(list(counts = cts.ZMZ_G_F_1    ))
sce.ZMZ_G_M_1     <- SingleCellExperiment(list(counts = cts.ZMZ_G_M_1    ))

sce_raw.ZMZ_Y_F_1 <- SingleCellExperiment(list(counts = cts_raw.ZMZ_Y_F_1    ))
sce_raw.ZMZ_Y_M_1 <- SingleCellExperiment(list(counts = cts_raw.ZMZ_Y_M_1    ))
sce_raw.ZMZ_M_F_1 <- SingleCellExperiment(list(counts = cts_raw.ZMZ_M_F_1    ))
sce_raw.ZMZ_M_M_1 <- SingleCellExperiment(list(counts = cts_raw.ZMZ_M_M_1    ))
sce_raw.ZMZ_O_F_1 <- SingleCellExperiment(list(counts = cts_raw.ZMZ_O_F_1    ))
sce_raw.ZMZ_O_M_1 <- SingleCellExperiment(list(counts = cts_raw.ZMZ_O_M_1    ))
sce_raw.ZMZ_G_F_1 <- SingleCellExperiment(list(counts = cts_raw.ZMZ_G_F_1    ))
sce_raw.ZMZ_G_M_1 <- SingleCellExperiment(list(counts = cts_raw.ZMZ_G_M_1    ))

# Run decontX
sce.ZMZ_Y_F_1 <- decontX(sce.ZMZ_Y_F_1, background = sce_raw.ZMZ_Y_F_1 )
sce.ZMZ_Y_M_1 <- decontX(sce.ZMZ_Y_M_1, background = sce_raw.ZMZ_Y_M_1 )
sce.ZMZ_M_F_1 <- decontX(sce.ZMZ_M_F_1, background = sce_raw.ZMZ_M_F_1 )
sce.ZMZ_M_M_1 <- decontX(sce.ZMZ_M_M_1, background = sce_raw.ZMZ_M_M_1 )
sce.ZMZ_O_F_1 <- decontX(sce.ZMZ_O_F_1, background = sce_raw.ZMZ_O_F_1 )
sce.ZMZ_O_M_1 <- decontX(sce.ZMZ_O_M_1, background = sce_raw.ZMZ_O_M_1 )
sce.ZMZ_G_F_1 <- decontX(sce.ZMZ_G_F_1, background = sce_raw.ZMZ_G_F_1 )
sce.ZMZ_G_M_1 <- decontX(sce.ZMZ_G_M_1, background = sce_raw.ZMZ_G_M_1 )

# get seurat objects
seurat.ZMZ_Y_F_1 <- CreateSeuratObject( round(decontXcounts(sce.ZMZ_Y_F_1)) )
seurat.ZMZ_Y_M_1 <- CreateSeuratObject( round(decontXcounts(sce.ZMZ_Y_M_1)) )
seurat.ZMZ_M_F_1 <- CreateSeuratObject( round(decontXcounts(sce.ZMZ_M_F_1)) )
seurat.ZMZ_M_M_1 <- CreateSeuratObject( round(decontXcounts(sce.ZMZ_M_M_1)) )
seurat.ZMZ_O_F_1 <- CreateSeuratObject( round(decontXcounts(sce.ZMZ_O_F_1)) )
seurat.ZMZ_O_M_1 <- CreateSeuratObject( round(decontXcounts(sce.ZMZ_O_M_1)) )
seurat.ZMZ_G_F_1 <- CreateSeuratObject( round(decontXcounts(sce.ZMZ_G_F_1)) )
seurat.ZMZ_G_M_1 <- CreateSeuratObject( round(decontXcounts(sce.ZMZ_G_M_1)) )

seurat.ZMZ_Y_F_1 <- AddMetaData(object = seurat.ZMZ_Y_F_1, colData(sce.ZMZ_Y_F_1)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.ZMZ_Y_M_1 <- AddMetaData(object = seurat.ZMZ_Y_M_1, colData(sce.ZMZ_Y_M_1)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.ZMZ_M_F_1 <- AddMetaData(object = seurat.ZMZ_M_F_1, colData(sce.ZMZ_M_F_1)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.ZMZ_M_M_1 <- AddMetaData(object = seurat.ZMZ_M_M_1, colData(sce.ZMZ_M_M_1)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.ZMZ_O_F_1 <- AddMetaData(object = seurat.ZMZ_O_F_1, colData(sce.ZMZ_O_F_1)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.ZMZ_O_M_1 <- AddMetaData(object = seurat.ZMZ_O_M_1, colData(sce.ZMZ_O_M_1)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.ZMZ_G_F_1 <- AddMetaData(object = seurat.ZMZ_G_F_1, colData(sce.ZMZ_G_F_1)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.ZMZ_G_M_1 <- AddMetaData(object = seurat.ZMZ_G_M_1, colData(sce.ZMZ_G_M_1)$decontX_contamination  , col.name = "decontX_contamination"  )

# Merge objects for the cohort
brain.c1.combined <- merge(seurat.ZMZ_Y_F_1,
                           y =  c(seurat.ZMZ_Y_M_1,
                                  seurat.ZMZ_M_F_1,
                                  seurat.ZMZ_M_M_1,
                                  seurat.ZMZ_O_F_1,
                                  seurat.ZMZ_O_M_1,
                                  seurat.ZMZ_G_F_1,
                                  seurat.ZMZ_G_M_1),
                           add.cell.ids = c("ZMZ_Y_F_1"  ,
                                            "ZMZ_Y_M_1"  ,
                                            "ZMZ_M_F_1" ,
                                            "ZMZ_M_M_1" ,
                                            "ZMZ_O_F_1" ,
                                            "ZMZ_O_M_1" ,
                                            "ZMZ_G_F_1" ,
                                            "ZMZ_G_M_1" ),
                           project = "10x_Killi_Brain_Aging")
brain.c1.combined
# An object of class Seurat 
# 27834 features across 47860 samples within 1 assay 
# Active assay: RNA (27834 features, 0 variable features)

#clean memory
rm(cts.ZMZ_Y_F_1, cts.ZMZ_Y_M_1, cts.ZMZ_M_F_1, cts.ZMZ_M_M_1, cts.ZMZ_O_F_1, cts.ZMZ_O_M_1, cts.ZMZ_G_F_1, cts.ZMZ_G_M_1, cts_raw.ZMZ_Y_F_1, cts_raw.ZMZ_Y_M_1, cts_raw.ZMZ_M_F_1, cts_raw.ZMZ_M_M_1, cts_raw.ZMZ_O_F_1, cts_raw.ZMZ_O_M_1, cts_raw.ZMZ_G_F_1, cts_raw.ZMZ_G_M_1, sce.ZMZ_Y_F_1, sce.ZMZ_Y_M_1, sce.ZMZ_M_F_1, sce.ZMZ_M_M_1, sce.ZMZ_O_F_1, sce.ZMZ_O_M_1, sce.ZMZ_G_F_1, sce.ZMZ_G_M_1, sce_raw.ZMZ_Y_F_1, sce_raw.ZMZ_Y_M_1, sce_raw.ZMZ_M_F_1, sce_raw.ZMZ_M_M_1, sce_raw.ZMZ_O_F_1, sce_raw.ZMZ_O_M_1, sce_raw.ZMZ_G_F_1, sce_raw.ZMZ_G_M_1, seurat.ZMZ_Y_F_1, seurat.ZMZ_Y_M_1, seurat.ZMZ_M_F_1, seurat.ZMZ_M_M_1, seurat.ZMZ_O_F_1, seurat.ZMZ_O_M_1, seurat.ZMZ_G_F_1, seurat.ZMZ_G_M_1)
################################################################################################


################################################################################################
#### 1b. Add key metadata to Seurat object

# create Group label
my.ZMZ_Y_F_1   <- grep("ZMZ_Y_F_1"  , colnames(brain.c1.combined@assays$RNA))
my.ZMZ_Y_M_1   <- grep("ZMZ_Y_M_1"  , colnames(brain.c1.combined@assays$RNA))
my.ZMZ_M_F_1   <- grep("ZMZ_M_F_1"  , colnames(brain.c1.combined@assays$RNA))
my.ZMZ_M_M_1   <- grep("ZMZ_M_M_1"  , colnames(brain.c1.combined@assays$RNA))
my.ZMZ_O_F_1   <- grep("ZMZ_O_F_1"  , colnames(brain.c1.combined@assays$RNA))
my.ZMZ_O_M_1   <- grep("ZMZ_O_M_1"  , colnames(brain.c1.combined@assays$RNA))
my.ZMZ_G_F_1   <- grep("ZMZ_G_F_1"  , colnames(brain.c1.combined@assays$RNA))
my.ZMZ_G_M_1   <- grep("ZMZ_G_M_1"  , colnames(brain.c1.combined@assays$RNA))

#####
Group <- rep("NA", length(colnames(brain.c1.combined@assays$RNA)))
Group[ my.ZMZ_Y_F_1  ]   <- "ZMZ_Y_F"
Group[ my.ZMZ_Y_M_1  ]   <- "ZMZ_Y_M"
Group[ my.ZMZ_M_F_1  ]   <- "ZMZ_M_F"
Group[ my.ZMZ_M_M_1  ]   <- "ZMZ_M_M"
Group[ my.ZMZ_O_F_1  ]   <- "ZMZ_O_F"
Group[ my.ZMZ_O_M_1  ]   <- "ZMZ_O_M"
Group[ my.ZMZ_G_F_1  ]   <- "ZMZ_G_F"
Group[ my.ZMZ_G_M_1  ]   <- "ZMZ_G_M"
Group <- data.frame(Group)
rownames(Group) <- colnames(brain.c1.combined@assays$RNA)

#####
Sex <- rep("NA", length(colnames(brain.c1.combined@assays$RNA)))
Sex[ my.ZMZ_Y_F_1 ]   <- "F"  
Sex[ my.ZMZ_Y_M_1 ]   <- "M" 
Sex[ my.ZMZ_M_F_1 ]   <- "F" 
Sex[ my.ZMZ_M_M_1 ]   <- "M"
Sex[ my.ZMZ_O_F_1 ]   <- "F" 
Sex[ my.ZMZ_O_M_1 ]   <- "M"  
Sex[ my.ZMZ_G_F_1 ]   <- "F"   
Sex[ my.ZMZ_G_M_1 ]   <- "M" 
Sex <- data.frame(Sex)
rownames(Sex) <- colnames(brain.c1.combined@assays$RNA)

##### Y, M, O, G
Age.gp <- rep("NA", length(colnames(brain.c1.combined@assays$RNA)))
Age.gp[ my.ZMZ_Y_F_1 ]   <- "Y"
Age.gp[ my.ZMZ_Y_M_1 ]   <- "Y"
Age.gp[ my.ZMZ_M_F_1 ]   <- "M"
Age.gp[ my.ZMZ_M_M_1 ]   <- "M"
Age.gp[ my.ZMZ_O_F_1 ]   <- "O"
Age.gp[ my.ZMZ_O_M_1 ]   <- "O"
Age.gp[ my.ZMZ_G_F_1 ]   <- "G"
Age.gp[ my.ZMZ_G_M_1 ]   <- "G"
Age.gp <- data.frame(Age.gp)
rownames(Age.gp) <- colnames(brain.c1.combined@assays$RNA)

##### (average of the pool) 
Age.w <- rep("NA", length(colnames(brain.c1.combined@assays$RNA)))
Age.w[ my.ZMZ_Y_F_1 ]   <-  5.9   #
Age.w[ my.ZMZ_Y_M_1 ]   <-  5.9   #
Age.w[ my.ZMZ_M_F_1 ]   <- 10.5   #
Age.w[ my.ZMZ_M_M_1 ]   <- 10.5   #
Age.w[ my.ZMZ_O_F_1 ]   <- 16.0   #
Age.w[ my.ZMZ_O_M_1 ]   <- 16.3   #
Age.w[ my.ZMZ_G_F_1 ]   <- 26.1   #
Age.w[ my.ZMZ_G_M_1 ]   <- 26.1   #
Age.w <- data.frame(Age.w)
rownames(Age.w) <- colnames(brain.c1.combined@assays$RNA)

#####
Strain <- rep("NA", length(colnames(brain.c1.combined@assays$RNA)))
Strain[ my.ZMZ_Y_F_1 ]   <- "ZMZ"
Strain[ my.ZMZ_Y_M_1 ]   <- "ZMZ"
Strain[ my.ZMZ_M_F_1 ]   <- "ZMZ"
Strain[ my.ZMZ_M_M_1 ]   <- "ZMZ"
Strain[ my.ZMZ_O_F_1 ]   <- "ZMZ"
Strain[ my.ZMZ_O_M_1 ]   <- "ZMZ"
Strain[ my.ZMZ_G_F_1 ]   <- "ZMZ"
Strain[ my.ZMZ_G_M_1 ]   <- "ZMZ"
Strain <- data.frame(Strain)
rownames(Strain) <- colnames(brain.c1.combined@assays$RNA)


#####
Batch <- rep("Set_1", length(colnames(brain.c1.combined@assays$RNA)))
Batch <- data.frame(Batch)
rownames(Batch) <- colnames(brain.c1.combined@assays$RNA)

# update Seurat with metadata
brain.c1.combined <- AddMetaData(object = brain.c1.combined, metadata = as.vector(Group)    , col.name = "Group"       )
brain.c1.combined <- AddMetaData(object = brain.c1.combined, metadata = as.vector(Sex)      , col.name = "Sex"         )
brain.c1.combined <- AddMetaData(object = brain.c1.combined, metadata = as.vector(Age.w)    , col.name = "Age_weeks"   )
brain.c1.combined <- AddMetaData(object = brain.c1.combined, metadata = as.vector(Age.gp)   , col.name = "Age_Group"   )
brain.c1.combined <- AddMetaData(object = brain.c1.combined, metadata = as.vector(Strain)   , col.name = "Strain"      )
brain.c1.combined <- AddMetaData(object = brain.c1.combined, metadata = as.vector(Batch)    , col.name = "Batch"       )
################################################################################################


################################################################################################
#### 2. Basic QC and filtering with Seurat and decontX

### No filtering on genes at this stage - only after all cohorts merged for fairness
brain.c1.combined <- SetIdent(brain.c1.combined, value = "Group")

# DecontX contamination levels for filtration
pdf(paste(Sys.Date(),"Killifish_Brain_Cohort1_violinPlots_QC_DecontX.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = brain.c1.combined, features = c("decontX_contamination"), pt.size = 0)
dev.off()

# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
brain.c1.combined[["percent.mito"]] <- PercentageFeatureSet(brain.c1.combined, pattern = "^MT-")
head(brain.c1.combined@meta.data)

pdf(paste(Sys.Date(),"Killifish_Brain_Cohort1_violinPlots_QC_gene_UMI_mito.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = brain.c1.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(brain.c1.combined, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(brain.c1.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste(Sys.Date(),"Killifish_Brain_Cohort1_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
plot1 + plot2
dev.off()

# filter dead/low Q cells
brain.c1.combined <- subset(brain.c1.combined, subset = nFeature_RNA > 250 & nFeature_RNA < 5000 & percent.mito < 10 & nCount_RNA < 25000 & decontX_contamination < 0.25 )
brain.c1.combined
# An object of class Seurat 
# 27834 features across 41842 samples within 1 assay 
# Active assay: RNA (27834 features, 0 variable features)

### Check data after cell filtering
head(brain.c1.combined@meta.data)

table(brain.c1.combined@meta.data$Group)
# ZMZ_G_F ZMZ_G_M ZMZ_M_F ZMZ_M_M ZMZ_O_F ZMZ_O_M ZMZ_Y_F ZMZ_Y_M
#     4271    4887    5681    5018    4765    5599    5201    6420 

#### Normalize the data for doublet analysis, etc
# global-scaling normalization method 'LogNormalize' normalizes gene expression measurements for each cell 
# by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
brain.c1.combined <- NormalizeData(object = brain.c1.combined, normalization.method = "LogNormalize",  scale.factor = 10000)
################################################################################################

################################################################################################
#### 2b. Cell cycle prediction and storage
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "../../cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")

# make into mouse gene names
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

cc.genes.mouse <- firstup(tolower(cc.genes))

# We can segregate this list into markers of G2/M phase and markers of S
# phase
s.genes   <- cc.genes.mouse[1:43]
g2m.genes <- cc.genes.mouse[44:97]

# Read parsed Blast results for killi 2015 genome to Mouse
my.nf_mm.homo <- read.csv('../../../Mouse_alignment/2022-10-11_Mouse_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-3_Minimal_HOMOLOGY_TABLE_REV.txt', sep = "\t", header = T)

# get the cell cycle gene homologs
s.genes.k   <- my.nf_mm.homo$Nfur_Symbol[my.nf_mm.homo$Mmu_Symbol %in% s.genes  ]
g2m.genes.k <- my.nf_mm.homo$Nfur_Symbol[my.nf_mm.homo$Mmu_Symbol %in% g2m.genes]

# Assign Cell-Cycle Scores
brain.c1.combined <- CellCycleScoring(object = brain.c1.combined, s.features = s.genes.k, g2m.features = g2m.genes.k, set.ident = TRUE)

# write predictions to file
write.table(brain.c1.combined@meta.data, file = paste0(Sys.Date(),"_Killi_Brain_Cohort1_CellCycle_predictions.txt"), sep = "\t", quote = F)
################################################################################################

################################################################################################
#### 3. Find and remove doublets using doublet finder & scds workflow

# https://github.com/chris-mcginnis-ucsf/DoubletFinder
brain.c1.combined <- SCTransform(object = brain.c1.combined, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mito", "Phase"))
save(brain.c1.combined, file = paste0(Sys.Date(),"_Killi_Brain_Cohort1_Seurat_object_postSCT.RData"))

# Run first pass analysis just for doublet identification (not final clustering)
brain.c1.combined <- RunPCA(brain.c1.combined, npcs = 30)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_Killi_Brain_Cohort1_ElbowPlot.pdf"))
ElbowPlot(brain.c1.combined, ndims = 30)
dev.off()

# run dimensionality reduction
# Keep all PCs here, we'll do the clean clustering analysis on the merged object across all cohorts
brain.c1.combined <- RunUMAP(brain.c1.combined, dims = 1:30)
brain.c1.combined <- FindNeighbors(brain.c1.combined, dims = 1:30)
brain.c1.combined <- FindClusters(object = brain.c1.combined)

#### need to split by 10x sample to make sure to identify real doublets
# will run on one object at a time
cohort.1.list <- SplitObject(brain.c1.combined, split.by = "Group")

## Assume doublet rate based on 10x information (add a 15% fudge factor due to nuclei being more sticky)
pred.dblt.rate <- 1.15 * predict(pred_dblt_lm, data.frame("cell_number" = unlist(lapply(cohort.1.list, ncol))))/100


## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& A. Run DoubletFinder
# loop over samples
for (i in 1:length(cohort.1.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list_killi <- paramSweep_v3(cohort.1.list[[i]], PCs = 1:30, sct = TRUE, num.cores	 = 4)
  sweep.stats_killi    <- summarizeSweep(sweep.res.list_killi, GT = FALSE)
  bcmvn_killi          <- find.pK(sweep.stats_killi)
  
  # need some R gymnastics since the Pk is stored as a factor for some reason
  # to get the pK number, need to first convert to character and THEN to numeric
  # numeric first yield row number
  pk.killi <- as.numeric(as.character(bcmvn_killi[as.numeric(bcmvn_killi$pK[bcmvn_killi$BCmetric == max(bcmvn_killi$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(cohort.1.list[[i]]@meta.data$seurat_clusters)             ## ex: annotations
  nExp_poi       <- round((pred.dblt.rate[i]) *length(cohort.1.list[[i]]@meta.data$Group))      ## Assume doublets based on nuclei isolation protocol performance
  
  ## Run DoubletFinder with varying classification stringencies
  cohort.1.list[[i]] <- doubletFinder_v3(cohort.1.list[[i]], PCs = 1:30, pN = 0.25, pK = pk.killi, nExp = nExp_poi,     reuse.pANN = FALSE, sct = T)
  
  # get classification name
  my.DF.res.col <- colnames(cohort.1.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(cohort.1.list[[i]]@meta.data))]
  
  # rename column to enable subsetting
  colnames(cohort.1.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(cohort.1.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# run UMAP plots
for (i in 1:length(cohort.1.list)) {
  
  pdf(paste(Sys.Date(),"Killifish_Tissue",names(cohort.1.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(cohort.1.list[[i]], reduction = "umap", group.by = "DoubletFinder"), raster = T)
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
killi.singlets.annot.c1 <- merge(cohort.1.list[[1]],
                                 y = c(cohort.1.list[[2]],
                                       cohort.1.list[[3]],
                                       cohort.1.list[[4]],
                                       cohort.1.list[[5]],
                                       cohort.1.list[[6]],
                                       cohort.1.list[[7]],
                                       cohort.1.list[[8]]),
                                 project = "Killi_Brain_Cohort1")
killi.singlets.annot.c1


# remove pANN columns that are 10xGenomics library lane specific
killi.singlets.annot.c1@meta.data <- killi.singlets.annot.c1@meta.data[,-grep("pANN",colnames(killi.singlets.annot.c1@meta.data))]


## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& B. Run scds:single cell doublet scoring (hybrid method)
# cxds is based on co-expression of gene pairs and works with absence/presence calls only, 
# bcds uses the full count information and a binary classification approach using artificially generated doublets. 
# cxds_bcds_hybrid combines both approaches

# create scds working object - convert list to SingleCellExperiment
cohort.1.list.scds        <- lapply(cohort.1.list, as.SingleCellExperiment)

# loop over sample
for (i in 1:length(cohort.1.list.scds)) {
  
  # Annotate doublets using co-expression based doublet scoring:
  cohort.1.list.scds[[i]] <- cxds_bcds_hybrid(cohort.1.list.scds[[i]])
  
  # predicted doublet rate
  n.db <- round((pred.dblt.rate[i])*ncol(cohort.1.list.scds[[i]]))                         ## Assume doublets based on nuclei isolation protocol performance
  
  # sort prediction, get top n.db cells
  srt.db.score <- sort(cohort.1.list.scds[[i]]$hybrid_score, index.return = T, decreasing = T)
  cohort.1.list.scds[[i]]$scds <- "Singlet"
  cohort.1.list.scds[[i]]$scds[srt.db.score$ix[1:n.db]] <- "Doublet"
  
}

# run UMAP plots
for (i in 1:length(cohort.1.list.scds)) {
  
  p <- plotReducedDim(cohort.1.list.scds[[i]], dimred = "UMAP", colour_by = "scds")
  
  pdf(paste(Sys.Date(),"Killifish_Tissue",names(cohort.1.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to doubletFinder annotated Seurat object
killi.singlets.annot.c1@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(cohort.1.list.scds)) {
  
  # for each object compare and move doublet annotations over
  killi.singlets.annot.c1@meta.data[colnames(cohort.1.list.scds[[i]]), ]$scds_hybrid <- cohort.1.list.scds[[i]]$scds
  
}

## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& C. Merge and summarize doublet findings

table(killi.singlets.annot.c1@meta.data$DoubletFinder, killi.singlets.annot.c1@meta.data$scds_hybrid)
#            Doublet Singlet
#    Doublet     205    1837
#    Singlet    1837   37963

# Union (more conservative)
killi.singlets.annot.c1@meta.data$DoubletCall <- ifelse( bitOr(killi.singlets.annot.c1@meta.data$DoubletFinder == "Doublet", killi.singlets.annot.c1@meta.data$scds_hybrid == "Doublet") > 0, 
                                                         "Doublet", "Singlet")
table(killi.singlets.annot.c1@meta.data$DoubletCall)
# Doublet Singlet 
#     3879   37963 

# re-run dimensionality reduction for plotting purposes
killi.singlets.annot.c1 <- SCTransform(object = killi.singlets.annot.c1, vars.to.regress =  c("nFeature_RNA", "nCount_RNA", "percent.mito"))
killi.singlets.annot.c1 <- RunPCA(killi.singlets.annot.c1, npcs = 30)
killi.singlets.annot.c1 <- RunUMAP(killi.singlets.annot.c1, dims = 1:30)

pdf(paste0(Sys.Date(),"_Killi_Brain_Cohort1_UMAP_Singlets_labelled_UNION.pdf"), width = 6, height = 5)
DimPlot(killi.singlets.annot.c1, reduction = "umap", group.by = "DoubletCall", raster = T)
dev.off()

pdf(paste0(Sys.Date(),"_Killi_Brain_Cohort1_Key_Marker_Gene_expression_plots.pdf"), width = 3.5, height = 3)
FeaturePlot(killi.singlets.annot.c1, features = "olig2", raster = T)
FeaturePlot(killi.singlets.annot.c1, features = "olig1", raster = T)
FeaturePlot(killi.singlets.annot.c1, features = "mpz"  , raster = T)
FeaturePlot(killi.singlets.annot.c1, features = "csf1r", raster = T)
dev.off()

# DimPlot(killi.singlets.annot.c1, reduction = "umap", group.by = "DoubletFinder")
# DimPlot(killi.singlets.annot.c1, reduction = "umap", group.by = "scds_hybrid")

# save annotated object
save(killi.singlets.annot.c1, file = paste0(Sys.Date(),"_Killifish_Brain_Cohort1_Seurat_object_with_AnnotatedDoublets.RData"))


### extract/subset only singlets
# save data for singlets df
killi.singlets.c1   <- subset(killi.singlets.annot.c1, subset = DoubletCall %in% "Singlet")  # only keep singlets
killi.singlets.c1
# An object of class Seurat 
# 52640 features across 37963 samples within 2 assays 
# Active assay: SCT (24806 features, 3000 variable features)
#  1 other assay present: RNA
#  2 dimensional reductions calculated: pca, umap

pdf(paste0(Sys.Date(),"_Killi_Brain_Cohort1_UMAP_Singlets_ONLY_UNION.pdf"), width = 6, height = 5)
DimPlot(killi.singlets.c1, reduction = "umap", group.by = "DoubletCall", raster = T)
dev.off()

table(killi.singlets.c1@meta.data$Group)
# ZMZ_G_F ZMZ_G_M ZMZ_M_F ZMZ_M_M ZMZ_O_F ZMZ_O_M ZMZ_Y_F ZMZ_Y_M 
#      3947    4472    5108    4567    4358    5036    4736    5739 


# save filtered/annotated object
save(killi.singlets.c1, file = paste0(Sys.Date(),"_Killifish_Brain_Cohort1_Seurat_object_SINGLETS_ONLY.RData"))
################################################################################################


#######################
sink(file = paste(Sys.Date(),"_Killifish_Brain_Cohort1_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()




