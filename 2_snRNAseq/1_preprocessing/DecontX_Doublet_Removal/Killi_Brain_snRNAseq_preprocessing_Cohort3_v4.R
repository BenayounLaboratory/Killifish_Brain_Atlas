setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Preprocessing/DecontX/Cohort_3/')
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
# make predicted doublet number 1.1 fold what 10x says, since we know nuclei are more susceptible
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
cts.ZMZ_Y_F_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/ZMZ_5w_F_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.ZMZ_Y_M_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/ZMZ_5w_M_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.ZMZ_M_F_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/ZMZ_10w_F_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.ZMZ_M_M_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/ZMZ_10w_M_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.ZMZ_O_F_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/ZMZ_15w_F_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.ZMZ_O_M_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/ZMZ_15w_M_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.ZMZ_G_F_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/ZMZ_26w_F_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.ZMZ_G_M_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/ZMZ_26w_M_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.GRZ_Y_F_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/GRZ_5w_F_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.GRZ_Y_M_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/GRZ_5w_M_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.GRZ_M_F_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/GRZ_10w_F_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.GRZ_M_M_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/GRZ_10w_M_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.GRZ_O_F_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/GRZ_15w_F_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.GRZ_O_M_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/GRZ_15w_M_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')

cts_raw.ZMZ_Y_F_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/ZMZ_5w_F_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts_raw.ZMZ_Y_M_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/ZMZ_5w_M_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts_raw.ZMZ_M_F_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/ZMZ_10w_F_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts_raw.ZMZ_M_M_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/ZMZ_10w_M_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts_raw.ZMZ_O_F_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/ZMZ_15w_F_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts_raw.ZMZ_O_M_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/ZMZ_15w_M_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts_raw.ZMZ_G_F_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/ZMZ_26w_F_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts_raw.ZMZ_G_M_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/ZMZ_26w_M_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts_raw.GRZ_Y_F_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/GRZ_5w_F_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts_raw.GRZ_Y_M_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/GRZ_5w_M_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts_raw.GRZ_M_F_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/GRZ_10w_F_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts_raw.GRZ_M_M_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/GRZ_10w_M_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts_raw.GRZ_O_F_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/GRZ_15w_F_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts_raw.GRZ_O_M_3  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-10-01_snRNAseq_brain_aging_set_3/CellRanger/GRZ_15w_M_3_FishTEDB_NR/outs/filtered_feature_bc_matrix/')

# Create SingleCellExperiment object
sce.ZMZ_Y_F_3  <- SingleCellExperiment(list(counts = cts.ZMZ_Y_F_3   ))
sce.ZMZ_Y_M_3  <- SingleCellExperiment(list(counts = cts.ZMZ_Y_M_3   ))
sce.ZMZ_M_F_3  <- SingleCellExperiment(list(counts = cts.ZMZ_M_F_3   ))
sce.ZMZ_M_M_3  <- SingleCellExperiment(list(counts = cts.ZMZ_M_M_3   ))
sce.ZMZ_O_F_3  <- SingleCellExperiment(list(counts = cts.ZMZ_O_F_3   ))
sce.ZMZ_O_M_3  <- SingleCellExperiment(list(counts = cts.ZMZ_O_M_3   ))
sce.ZMZ_G_F_3  <- SingleCellExperiment(list(counts = cts.ZMZ_G_F_3   ))
sce.ZMZ_G_M_3  <- SingleCellExperiment(list(counts = cts.ZMZ_G_M_3   ))
sce.GRZ_Y_F_3  <- SingleCellExperiment(list(counts = cts.GRZ_Y_F_3   ))
sce.GRZ_Y_M_3  <- SingleCellExperiment(list(counts = cts.GRZ_Y_M_3   ))
sce.GRZ_M_F_3  <- SingleCellExperiment(list(counts = cts.GRZ_M_F_3   ))
sce.GRZ_M_M_3  <- SingleCellExperiment(list(counts = cts.GRZ_M_M_3   ))
sce.GRZ_O_F_3  <- SingleCellExperiment(list(counts = cts.GRZ_O_F_3   ))
sce.GRZ_O_M_3  <- SingleCellExperiment(list(counts = cts.GRZ_O_M_3   ))

sce_raw.ZMZ_Y_F_3  <- SingleCellExperiment(list(counts = cts_raw.ZMZ_Y_F_3   ))
sce_raw.ZMZ_Y_M_3  <- SingleCellExperiment(list(counts = cts_raw.ZMZ_Y_M_3   ))
sce_raw.ZMZ_M_F_3  <- SingleCellExperiment(list(counts = cts_raw.ZMZ_M_F_3   ))
sce_raw.ZMZ_M_M_3  <- SingleCellExperiment(list(counts = cts_raw.ZMZ_M_M_3   ))
sce_raw.ZMZ_O_F_3  <- SingleCellExperiment(list(counts = cts_raw.ZMZ_O_F_3   ))
sce_raw.ZMZ_O_M_3  <- SingleCellExperiment(list(counts = cts_raw.ZMZ_O_M_3   ))
sce_raw.ZMZ_G_F_3  <- SingleCellExperiment(list(counts = cts_raw.ZMZ_G_F_3   ))
sce_raw.ZMZ_G_M_3  <- SingleCellExperiment(list(counts = cts_raw.ZMZ_G_M_3   ))
sce_raw.GRZ_Y_F_3  <- SingleCellExperiment(list(counts = cts_raw.GRZ_Y_F_3   ))
sce_raw.GRZ_Y_M_3  <- SingleCellExperiment(list(counts = cts_raw.GRZ_Y_M_3   ))
sce_raw.GRZ_M_F_3  <- SingleCellExperiment(list(counts = cts_raw.GRZ_M_F_3   ))
sce_raw.GRZ_M_M_3  <- SingleCellExperiment(list(counts = cts_raw.GRZ_M_M_3   ))
sce_raw.GRZ_O_F_3  <- SingleCellExperiment(list(counts = cts_raw.GRZ_O_F_3   ))
sce_raw.GRZ_O_M_3  <- SingleCellExperiment(list(counts = cts_raw.GRZ_O_M_3   ))

# Run decontX
sce.ZMZ_Y_F_3  <- decontX(sce.ZMZ_Y_F_3 , background = sce_raw.ZMZ_Y_F_3  )
sce.ZMZ_Y_M_3  <- decontX(sce.ZMZ_Y_M_3 , background = sce_raw.ZMZ_Y_M_3  )
sce.ZMZ_M_F_3  <- decontX(sce.ZMZ_M_F_3 , background = sce_raw.ZMZ_M_F_3  )
sce.ZMZ_M_M_3  <- decontX(sce.ZMZ_M_M_3 , background = sce_raw.ZMZ_M_M_3  )
sce.ZMZ_O_F_3  <- decontX(sce.ZMZ_O_F_3 , background = sce_raw.ZMZ_O_F_3  )
sce.ZMZ_O_M_3  <- decontX(sce.ZMZ_O_M_3 , background = sce_raw.ZMZ_O_M_3  )
sce.ZMZ_G_F_3  <- decontX(sce.ZMZ_G_F_3 , background = sce_raw.ZMZ_G_F_3  )
sce.ZMZ_G_M_3  <- decontX(sce.ZMZ_G_M_3 , background = sce_raw.ZMZ_G_M_3  )
sce.GRZ_Y_F_3  <- decontX(sce.GRZ_Y_F_3 , background = sce_raw.GRZ_Y_F_3  )
sce.GRZ_Y_M_3  <- decontX(sce.GRZ_Y_M_3 , background = sce_raw.GRZ_Y_M_3  )
sce.GRZ_M_F_3  <- decontX(sce.GRZ_M_F_3 , background = sce_raw.GRZ_M_F_3  )
sce.GRZ_M_M_3  <- decontX(sce.GRZ_M_M_3 , background = sce_raw.GRZ_M_M_3  )
sce.GRZ_O_F_3  <- decontX(sce.GRZ_O_F_3 , background = sce_raw.GRZ_O_F_3  )
sce.GRZ_O_M_3  <- decontX(sce.GRZ_O_M_3 , background = sce_raw.GRZ_O_M_3  )


# get seurat objects
seurat.ZMZ_Y_F_3  <- CreateSeuratObject( round(decontXcounts( sce.ZMZ_Y_F_3 ) ) )
seurat.ZMZ_Y_M_3  <- CreateSeuratObject( round(decontXcounts( sce.ZMZ_Y_M_3 ) ) )
seurat.ZMZ_M_F_3  <- CreateSeuratObject( round(decontXcounts( sce.ZMZ_M_F_3 ) ) )
seurat.ZMZ_M_M_3  <- CreateSeuratObject( round(decontXcounts( sce.ZMZ_M_M_3 ) ) )
seurat.ZMZ_O_F_3  <- CreateSeuratObject( round(decontXcounts( sce.ZMZ_O_F_3 ) ) )
seurat.ZMZ_O_M_3  <- CreateSeuratObject( round(decontXcounts( sce.ZMZ_O_M_3 ) ) )
seurat.ZMZ_G_F_3  <- CreateSeuratObject( round(decontXcounts( sce.ZMZ_G_F_3 ) ) )
seurat.ZMZ_G_M_3  <- CreateSeuratObject( round(decontXcounts( sce.ZMZ_G_M_3 ) ) )
seurat.GRZ_Y_F_3  <- CreateSeuratObject( round(decontXcounts( sce.GRZ_Y_F_3 ) ) )
seurat.GRZ_Y_M_3  <- CreateSeuratObject( round(decontXcounts( sce.GRZ_Y_M_3 ) ) )
seurat.GRZ_M_F_3  <- CreateSeuratObject( round(decontXcounts( sce.GRZ_M_F_3 ) ) )
seurat.GRZ_M_M_3  <- CreateSeuratObject( round(decontXcounts( sce.GRZ_M_M_3 ) ) )
seurat.GRZ_O_F_3  <- CreateSeuratObject( round(decontXcounts( sce.GRZ_O_F_3 ) ) )
seurat.GRZ_O_M_3  <- CreateSeuratObject( round(decontXcounts( sce.GRZ_O_M_3 ) ) )

seurat.ZMZ_Y_F_3  <- AddMetaData(object = seurat.ZMZ_Y_F_3, colData(sce.ZMZ_Y_F_3)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.ZMZ_Y_M_3  <- AddMetaData(object = seurat.ZMZ_Y_M_3, colData(sce.ZMZ_Y_M_3)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.ZMZ_M_F_3  <- AddMetaData(object = seurat.ZMZ_M_F_3, colData(sce.ZMZ_M_F_3)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.ZMZ_M_M_3  <- AddMetaData(object = seurat.ZMZ_M_M_3, colData(sce.ZMZ_M_M_3)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.ZMZ_O_F_3  <- AddMetaData(object = seurat.ZMZ_O_F_3, colData(sce.ZMZ_O_F_3)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.ZMZ_O_M_3  <- AddMetaData(object = seurat.ZMZ_O_M_3, colData(sce.ZMZ_O_M_3)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.ZMZ_G_F_3  <- AddMetaData(object = seurat.ZMZ_G_F_3, colData(sce.ZMZ_G_F_3)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.ZMZ_G_M_3  <- AddMetaData(object = seurat.ZMZ_G_M_3, colData(sce.ZMZ_G_M_3)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.GRZ_Y_F_3  <- AddMetaData(object = seurat.GRZ_Y_F_3, colData(sce.GRZ_Y_F_3)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.GRZ_Y_M_3  <- AddMetaData(object = seurat.GRZ_Y_M_3, colData(sce.GRZ_Y_M_3)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.GRZ_M_F_3  <- AddMetaData(object = seurat.GRZ_M_F_3, colData(sce.GRZ_M_F_3)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.GRZ_M_M_3  <- AddMetaData(object = seurat.GRZ_M_M_3, colData(sce.GRZ_M_M_3)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.GRZ_O_F_3  <- AddMetaData(object = seurat.GRZ_O_F_3, colData(sce.GRZ_O_F_3)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.GRZ_O_M_3  <- AddMetaData(object = seurat.GRZ_O_M_3, colData(sce.GRZ_O_M_3)$decontX_contamination  , col.name = "decontX_contamination"  )


# Merge objects for the cohort
# Merge objects for the cohort
brain.c3.combined <- merge(seurat.ZMZ_Y_F_3,
                           y =  c(seurat.ZMZ_Y_M_3,
                                  seurat.ZMZ_M_F_3,
                                  seurat.ZMZ_M_M_3,
                                  seurat.ZMZ_O_F_3,
                                  seurat.ZMZ_O_M_3,
                                  seurat.ZMZ_G_F_3,
                                  seurat.ZMZ_G_M_3,
                                  
                                  seurat.GRZ_Y_F_3,
                                  seurat.GRZ_Y_M_3,
                                  seurat.GRZ_M_F_3,
                                  seurat.GRZ_M_M_3,
                                  seurat.GRZ_O_F_3,
                                  seurat.GRZ_O_M_3),
                           add.cell.ids = c("ZMZ_Y_F_3"  ,
                                            "ZMZ_Y_M_3"  ,
                                            "ZMZ_M_F_3" ,
                                            "ZMZ_M_M_3" ,
                                            "ZMZ_O_F_3" ,
                                            "ZMZ_O_M_3" ,
                                            "ZMZ_G_F_3" ,
                                            "ZMZ_G_M_3" ,
                                            
                                            "GRZ_Y_F_3"  ,
                                            "GRZ_Y_M_3"  ,
                                            "GRZ_M_F_3" ,
                                            "GRZ_M_M_3" ,
                                            "GRZ_O_F_3" ,
                                            "GRZ_O_M_3"),
                           project = "10x_Killi_Brain_Aging")
brain.c3.combined
# An object of class Seurat
# 27834 features across 74454 samples within 1 assay
# Active assay: RNA (27834 features, 0 variable features)

# clean memory
rm(cts.ZMZ_Y_F_3,cts.ZMZ_Y_M_3,cts.ZMZ_M_F_3,cts.ZMZ_M_M_3,cts.ZMZ_O_F_3,cts.ZMZ_O_M_3,cts.ZMZ_G_F_3,cts.ZMZ_G_M_3,cts.GRZ_Y_F_3,cts.GRZ_Y_M_3,cts.GRZ_M_F_3,cts.GRZ_M_M_3,cts.GRZ_O_F_3,cts.GRZ_O_M_3,cts_raw.ZMZ_Y_F_3,cts_raw.ZMZ_Y_M_3,cts_raw.ZMZ_M_F_3,cts_raw.ZMZ_M_M_3,cts_raw.ZMZ_O_F_3,cts_raw.ZMZ_O_M_3,cts_raw.ZMZ_G_F_3,cts_raw.ZMZ_G_M_3,cts_raw.GRZ_Y_F_3,cts_raw.GRZ_Y_M_3,cts_raw.GRZ_M_F_3,cts_raw.GRZ_M_M_3,cts_raw.GRZ_O_F_3,cts_raw.GRZ_O_M_3,sce.ZMZ_Y_F_3,sce.ZMZ_Y_M_3,sce.ZMZ_M_F_3,sce.ZMZ_M_M_3,sce.ZMZ_O_F_3,sce.ZMZ_O_M_3,sce.ZMZ_G_F_3,sce.ZMZ_G_M_3,sce.GRZ_Y_F_3,sce.GRZ_Y_M_3,sce.GRZ_M_F_3,sce.GRZ_M_M_3,sce.GRZ_O_F_3,sce.GRZ_O_M_3,sce_raw.ZMZ_Y_F_3,sce_raw.ZMZ_Y_M_3,sce_raw.ZMZ_M_F_3,sce_raw.ZMZ_M_M_3,sce_raw.ZMZ_O_F_3,sce_raw.ZMZ_O_M_3,sce_raw.ZMZ_G_F_3,sce_raw.ZMZ_G_M_3,sce_raw.GRZ_Y_F_3,sce_raw.GRZ_Y_M_3,sce_raw.GRZ_M_F_3,sce_raw.GRZ_M_M_3,sce_raw.GRZ_O_F_3,sce_raw.GRZ_O_M_3,seurat.ZMZ_Y_F_3,seurat.ZMZ_Y_M_3,seurat.ZMZ_M_F_3,seurat.ZMZ_M_M_3,seurat.ZMZ_O_F_3,seurat.ZMZ_O_M_3,seurat.ZMZ_G_F_3,seurat.ZMZ_G_M_3,seurat.GRZ_Y_F_3,seurat.GRZ_Y_M_3,seurat.GRZ_M_F_3,seurat.GRZ_M_M_3,seurat.GRZ_O_F_3,seurat.GRZ_O_M_3)
################################################################################################


################################################################################################
#### 1b. Add key metadata to Seurat object

# create Group label
my.ZMZ_Y_F_3   <- grep("ZMZ_Y_F_3"  , colnames(brain.c3.combined@assays$RNA))
my.ZMZ_Y_M_3   <- grep("ZMZ_Y_M_3"  , colnames(brain.c3.combined@assays$RNA))
my.ZMZ_M_F_3   <- grep("ZMZ_M_F_3"  , colnames(brain.c3.combined@assays$RNA))
my.ZMZ_M_M_3   <- grep("ZMZ_M_M_3"  , colnames(brain.c3.combined@assays$RNA))
my.ZMZ_O_F_3   <- grep("ZMZ_O_F_3"  , colnames(brain.c3.combined@assays$RNA))
my.ZMZ_O_M_3   <- grep("ZMZ_O_M_3"  , colnames(brain.c3.combined@assays$RNA))
my.ZMZ_G_F_3   <- grep("ZMZ_G_F_3"  , colnames(brain.c3.combined@assays$RNA))
my.ZMZ_G_M_3   <- grep("ZMZ_G_M_3"  , colnames(brain.c3.combined@assays$RNA))

my.GRZ_Y_F_3   <- grep("GRZ_Y_F_3"  , colnames(brain.c3.combined@assays$RNA))
my.GRZ_Y_M_3   <- grep("GRZ_Y_M_3"  , colnames(brain.c3.combined@assays$RNA))
my.GRZ_M_F_3   <- grep("GRZ_M_F_3"  , colnames(brain.c3.combined@assays$RNA))
my.GRZ_M_M_3   <- grep("GRZ_M_M_3"  , colnames(brain.c3.combined@assays$RNA))
my.GRZ_O_F_3   <- grep("GRZ_O_F_3"  , colnames(brain.c3.combined@assays$RNA))
my.GRZ_O_M_3   <- grep("GRZ_O_M_3"  , colnames(brain.c3.combined@assays$RNA))

##### (even though we initially labeled 5/15w, Ari used 6/16w - correct in meta data attribution)
Group <- rep("NA", length(colnames(brain.c3.combined@assays$RNA)))
Group[ my.ZMZ_Y_F_3  ]   <- "ZMZ_Y_F"
Group[ my.ZMZ_Y_M_3  ]   <- "ZMZ_Y_M"
Group[ my.ZMZ_M_F_3  ]   <- "ZMZ_M_F"
Group[ my.ZMZ_M_M_3  ]   <- "ZMZ_M_M"
Group[ my.ZMZ_O_F_3  ]   <- "ZMZ_O_F"
Group[ my.ZMZ_O_M_3  ]   <- "ZMZ_O_M"
Group[ my.ZMZ_G_F_3  ]   <- "ZMZ_G_F"
Group[ my.ZMZ_G_M_3  ]   <- "ZMZ_G_M"
Group[ my.GRZ_Y_F_3  ]   <- "GRZ_Y_F"
Group[ my.GRZ_Y_M_3  ]   <- "GRZ_Y_M"
Group[ my.GRZ_M_F_3  ]   <- "GRZ_M_F"
Group[ my.GRZ_M_M_3  ]   <- "GRZ_M_M"
Group[ my.GRZ_O_F_3  ]   <- "GRZ_O_F"
Group[ my.GRZ_O_M_3  ]   <- "GRZ_O_M"
Group <- data.frame(Group)
rownames(Group) <- colnames(brain.c3.combined@assays$RNA)

#####
Sex <- rep("NA", length(colnames(brain.c3.combined@assays$RNA)))
Sex[ my.ZMZ_Y_F_3 ]   <- "F"  
Sex[ my.ZMZ_Y_M_3 ]   <- "M" 
Sex[ my.ZMZ_M_F_3 ]   <- "F" 
Sex[ my.ZMZ_M_M_3 ]   <- "M"
Sex[ my.ZMZ_O_F_3 ]   <- "F" 
Sex[ my.ZMZ_O_M_3 ]   <- "M"  
Sex[ my.ZMZ_G_F_3 ]   <- "F"   
Sex[ my.ZMZ_G_M_3 ]   <- "M" 
Sex[ my.GRZ_Y_F_3 ]   <- "F"  
Sex[ my.GRZ_Y_M_3 ]   <- "M" 
Sex[ my.GRZ_M_F_3 ]   <- "F" 
Sex[ my.GRZ_M_M_3 ]   <- "M"
Sex[ my.GRZ_O_F_3 ]   <- "F" 
Sex[ my.GRZ_O_M_3 ]   <- "M"  
Sex <- data.frame(Sex)
rownames(Sex) <- colnames(brain.c3.combined@assays$RNA)

##### Y, M, O, G
Age.gp <- rep("NA", length(colnames(brain.c3.combined@assays$RNA)))
Age.gp[ my.ZMZ_Y_F_3 ]   <- "Y"
Age.gp[ my.ZMZ_Y_M_3 ]   <- "Y"
Age.gp[ my.ZMZ_M_F_3 ]   <- "M"
Age.gp[ my.ZMZ_M_M_3 ]   <- "M"
Age.gp[ my.ZMZ_O_F_3 ]   <- "O"
Age.gp[ my.ZMZ_O_M_3 ]   <- "O"
Age.gp[ my.ZMZ_G_F_3 ]   <- "G"
Age.gp[ my.ZMZ_G_M_3 ]   <- "G"
Age.gp[ my.GRZ_Y_F_3 ]   <- "Y"
Age.gp[ my.GRZ_Y_M_3 ]   <- "Y"
Age.gp[ my.GRZ_M_F_3 ]   <- "M"
Age.gp[ my.GRZ_M_M_3 ]   <- "M"
Age.gp[ my.GRZ_O_F_3 ]   <- "O"
Age.gp[ my.GRZ_O_M_3 ]   <- "O"
Age.gp <- data.frame(Age.gp)
rownames(Age.gp) <- colnames(brain.c3.combined@assays$RNA)

##### (average of the pool) 
Age.w <- rep("NA", length(colnames(brain.c3.combined@assays$RNA)))
Age.w[ my.ZMZ_Y_F_3 ]   <-   5.9  # 
Age.w[ my.ZMZ_Y_M_3 ]   <-   5.9  # 
Age.w[ my.ZMZ_M_F_3 ]   <-  10.4  # 
Age.w[ my.ZMZ_M_M_3 ]   <-  10.6  # 
Age.w[ my.ZMZ_O_F_3 ]   <-  15.9  # 
Age.w[ my.ZMZ_O_M_3 ]   <-  16.3  # 
Age.w[ my.ZMZ_G_F_3 ]   <-  26.2  # 
Age.w[ my.ZMZ_G_M_3 ]   <-  26.2  # 
Age.w[ my.GRZ_Y_F_3 ]   <-   6.5  # 
Age.w[ my.GRZ_Y_M_3 ]   <-   6.5  # 
Age.w[ my.GRZ_M_F_3 ]   <-  10.5  # 
Age.w[ my.GRZ_M_M_3 ]   <-  10.1  # 
Age.w[ my.GRZ_O_F_3 ]   <-  16.4  # 
Age.w[ my.GRZ_O_M_3 ]   <-  16.4  # 
Age.w <- data.frame(Age.w)
rownames(Age.w) <- colnames(brain.c3.combined@assays$RNA)

#####
Strain <- rep("NA", length(colnames(brain.c3.combined@assays$RNA)))
Strain[ my.ZMZ_Y_F_3 ]   <- "ZMZ"
Strain[ my.ZMZ_Y_M_3 ]   <- "ZMZ"
Strain[ my.ZMZ_M_F_3 ]   <- "ZMZ"
Strain[ my.ZMZ_M_M_3 ]   <- "ZMZ"
Strain[ my.ZMZ_O_F_3 ]   <- "ZMZ"
Strain[ my.ZMZ_O_M_3 ]   <- "ZMZ"
Strain[ my.ZMZ_G_F_3 ]   <- "ZMZ"
Strain[ my.ZMZ_G_M_3 ]   <- "ZMZ"
Strain[ my.GRZ_Y_F_3 ]   <- "GRZ"
Strain[ my.GRZ_Y_M_3 ]   <- "GRZ"
Strain[ my.GRZ_M_F_3 ]   <- "GRZ"
Strain[ my.GRZ_M_M_3 ]   <- "GRZ"
Strain[ my.GRZ_O_F_3 ]   <- "GRZ"
Strain[ my.GRZ_O_M_3 ]   <- "GRZ"
Strain <- data.frame(Strain)
rownames(Strain) <- colnames(brain.c3.combined@assays$RNA)


#####
Batch <- rep("Set_3", length(colnames(brain.c3.combined@assays$RNA)))
Batch <- data.frame(Batch)
rownames(Batch) <- colnames(brain.c3.combined@assays$RNA)

# update Seurat with metadata
brain.c3.combined <- AddMetaData(object = brain.c3.combined, metadata = as.vector(Group)    , col.name = "Group"       )
brain.c3.combined <- AddMetaData(object = brain.c3.combined, metadata = as.vector(Sex)      , col.name = "Sex"         )
brain.c3.combined <- AddMetaData(object = brain.c3.combined, metadata = as.vector(Age.w)    , col.name = "Age_weeks"   )
brain.c3.combined <- AddMetaData(object = brain.c3.combined, metadata = as.vector(Age.gp)   , col.name = "Age_Group"   )
brain.c3.combined <- AddMetaData(object = brain.c3.combined, metadata = as.vector(Strain)   , col.name = "Strain"      )
brain.c3.combined <- AddMetaData(object = brain.c3.combined, metadata = as.vector(Batch)    , col.name = "Batch"       )
################################################################################################


################################################################################################
#### 2. Basic QC and filtering with Seurat

### No filtering on genes at this stage - only after all cohorts merged for fairness
brain.c3.combined <- SetIdent(brain.c3.combined, value = "Group")

# DecontX contamination levels for filtration
pdf(paste(Sys.Date(),"Killifish_Brain_Cohort3_violinPlots_QC_DecontX.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = brain.c3.combined, features = c("decontX_contamination"), pt.size = 0)
dev.off()

# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
brain.c3.combined[["percent.mito"]] <- PercentageFeatureSet(brain.c3.combined, pattern = "^MT-")

pdf(paste(Sys.Date(),"Killifish_Brain_Cohort3_violinPlots_QC_gene_UMI_mito.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = brain.c3.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(brain.c3.combined, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(brain.c3.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste(Sys.Date(),"Killifish_Brain_Cohort3_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
plot1 + plot2
dev.off()

# filter dead/low Q cells
brain.c3.combined <- subset(brain.c3.combined, subset = nFeature_RNA > 250 & nFeature_RNA < 5000 & percent.mito < 10 & nCount_RNA < 25000 & decontX_contamination < 0.25 )
brain.c3.combined
# An object of class Seurat
# 27834 features across 73780 samples within 1 assay
# Active assay: RNA (27834 features, 0 variable features)

### Check data after cell filtering
head(brain.c3.combined@meta.data)


table(brain.c3.combined@meta.data$Group)
# GRZ_M_F GRZ_M_M GRZ_O_F GRZ_O_M GRZ_Y_F GRZ_Y_M ZMZ_G_F ZMZ_G_M ZMZ_M_F ZMZ_M_M ZMZ_O_F ZMZ_O_M ZMZ_Y_F ZMZ_Y_M
#    6720    6067    5793    6197    6730    7471    1542    1787    5701    8406   4868    2209    4483    5806


#### Normalize the data for doublet analysis, etc
# global-scaling normalization method 'LogNormalize' normalizes gene expression measurements for each cell 
# by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
brain.c3.combined <- NormalizeData(object = brain.c3.combined, normalization.method = "LogNormalize",  scale.factor = 10000)
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
brain.c3.combined <- CellCycleScoring(object = brain.c3.combined, s.features = s.genes.k, g2m.features = g2m.genes.k, set.ident = TRUE)

# write predictions to file
write.table(brain.c3.combined@meta.data, file = paste0(Sys.Date(),"_Killi_Brain_Cohort3_CellCycle_predictions.txt"), sep = "\t", quote = F)
################################################################################################


################################################################################################
#### 3. Find and remove doublets using doublet finder & scds workflow

# pre normalize the data for doubletFinder
brain.c3.combined <- SCTransform(object = brain.c3.combined, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mito", "Phase"))
save(brain.c3.combined, file = paste0(Sys.Date(),"_Killi_Brain_Cohort3_Seurat_object_postSCT.RData"))

# Run first pass analysis just for doublet identification (not final clustering)
brain.c3.combined <- RunPCA(brain.c3.combined, npcs = 30)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_Killi_Brain_Cohort3_ElbowPlot.pdf"))
ElbowPlot(brain.c3.combined, ndims = 30)
dev.off()

# run dimensionality reduction
# Keep all PCs here, we'll do the clean clustering analysis on the merged object across all cohorts
brain.c3.combined <- RunUMAP(brain.c3.combined, dims = 1:30)
brain.c3.combined <- FindNeighbors(brain.c3.combined, dims = 1:30)
brain.c3.combined <- FindClusters(object = brain.c3.combined)

#### need to split by 10x sample to make sure to identify real doublets
# will run on one object at a time
cohort.3.list <- SplitObject(brain.c3.combined, split.by = "Group")

## Assume doublet rate based on 10x information (add a 15% fudge factor due to nuclei being more sticky)
pred.dblt.rate <- 1.15 * predict(pred_dblt_lm, data.frame("cell_number" = unlist(lapply(cohort.3.list, ncol))))/100

## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& A. Run DoubletFinder
# loop over samples
for (i in 1:length(cohort.3.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list_killi <- paramSweep_v3(cohort.3.list[[i]], PCs = 1:30, sct = TRUE, num.cores	 = 4)
  sweep.stats_killi    <- summarizeSweep(sweep.res.list_killi, GT = FALSE)
  bcmvn_killi          <- find.pK(sweep.stats_killi)
  
  # need some R gymnastics since the Pk is stored as a factor for some reason
  # to get the pK number, need to first convert to character and THEN to numeric
  # numeric first yield row number
  pk.killi <- as.numeric(as.character(bcmvn_killi[as.numeric(bcmvn_killi$pK[bcmvn_killi$BCmetric == max(bcmvn_killi$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(cohort.3.list[[i]]@meta.data$seurat_clusters)             ## ex: annotations
  nExp_poi       <- round((pred.dblt.rate[i]) *length(cohort.3.list[[i]]@meta.data$Group))      ## Assume doublets based on nuclei isolation protocol performance
  
  ## Run DoubletFinder with varying classification stringencies
  cohort.3.list[[i]] <- doubletFinder_v3(cohort.3.list[[i]], PCs = 1:30, pN = 0.25, pK = pk.killi, nExp = nExp_poi,     reuse.pANN = FALSE, sct = T)
  
  # get classification name
  my.DF.res.col <- colnames(cohort.3.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(cohort.3.list[[i]]@meta.data))]
  
  # rename column to enable subsetting
  colnames(cohort.3.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(cohort.3.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# run UMAP plots
for (i in 1:length(cohort.3.list)) {
  
  pdf(paste(Sys.Date(),"Killifish_Tissue",names(cohort.3.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(cohort.3.list[[i]], reduction = "umap", group.by = "DoubletFinder"), raster = T)
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
killi.singlets.annot.c3 <- merge(cohort.3.list[[1]],
                                 y = c(cohort.3.list[[ 2]],
                                       cohort.3.list[[ 3]],
                                       cohort.3.list[[ 4]],
                                       cohort.3.list[[ 5]],
                                       cohort.3.list[[ 6]],
                                       cohort.3.list[[ 7]],
                                       cohort.3.list[[ 8]],
                                       cohort.3.list[[ 9]],
                                       cohort.3.list[[10]],
                                       cohort.3.list[[11]],
                                       cohort.3.list[[12]],
                                       cohort.3.list[[13]],
                                       cohort.3.list[[14]]
                                       ),
                                 project = "Killi_Brain_Cohort3")
killi.singlets.annot.c3


# remove pANN columns that are 10xGenomics library lane specific
killi.singlets.annot.c3@meta.data <- killi.singlets.annot.c3@meta.data[,-grep("pANN",colnames(killi.singlets.annot.c3@meta.data))]


## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& B. Run scds:single cell doublet scoring (hybrid method)
# cxds is based on co-expression of gene pairs and works with absence/presence calls only, 
# bcds uses the full count information and a binary classification approach using artificially generated doublets. 
# cxds_bcds_hybrid combines both approaches

# create scds working object - convert list to SingleCellExperiment
cohort.3.list.scds        <- lapply(cohort.3.list, as.SingleCellExperiment)

# loop over sample
for (i in 1:length(cohort.3.list.scds)) {
  
  # Annotate doublets using co-expression based doublet scoring:
  cohort.3.list.scds[[i]] <- cxds_bcds_hybrid(cohort.3.list.scds[[i]])
  
  # predicted doublet rate
  n.db <- round((pred.dblt.rate[i])*ncol(cohort.3.list.scds[[i]]))                         ## Assume doublets based on nuclei isolation protocol performance
  
  # sort prediction, get top n.db cells
  srt.db.score <- sort(cohort.3.list.scds[[i]]$hybrid_score, index.return = T, decreasing = T)
  cohort.3.list.scds[[i]]$scds <- "Singlet"
  cohort.3.list.scds[[i]]$scds[srt.db.score$ix[1:n.db]] <- "Doublet"
  
}

# run UMAP plots
for (i in 1:length(cohort.3.list.scds)) {
  
  p <- plotReducedDim(cohort.3.list.scds[[i]], dimred = "UMAP", colour_by = "scds")
  
  pdf(paste(Sys.Date(),"Killifish_",names(cohort.3.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to doubletFinder annotated Seurat object
killi.singlets.annot.c3@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(cohort.3.list.scds)) {
  
  # for each object compare and move doublet annotations over
  killi.singlets.annot.c3@meta.data[colnames(cohort.3.list.scds[[i]]), ]$scds_hybrid <- cohort.3.list.scds[[i]]$scds
  
}

# free some memory
rm(cohort.3.list, cohort.3.list.scds, brain.c3.combined)

## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& C. Merge and summarize doublet findings

table(killi.singlets.annot.c3@meta.data$DoubletFinder, killi.singlets.annot.c3@meta.data$scds_hybrid)
#           Doublet Singlet
#   Doublet 
#   Singlet 

# Union (more conservative)
killi.singlets.annot.c3@meta.data$DoubletCall <- ifelse( bitOr(killi.singlets.annot.c3@meta.data$DoubletFinder == "Doublet", killi.singlets.annot.c3@meta.data$scds_hybrid == "Doublet") > 0, 
                                                           "Doublet", "Singlet")
table(killi.singlets.annot.c3@meta.data$DoubletCall)
#   Doublet Singlet
#    

# re-run dimensionality reduction for plotting purposes
killi.singlets.annot.c3 <- SCTransform(object = killi.singlets.annot.c3, vars.to.regress =  c("nFeature_RNA", "nCount_RNA", "percent.mito"))
killi.singlets.annot.c3 <- RunPCA(killi.singlets.annot.c3, npcs = 30)
killi.singlets.annot.c3 <- RunUMAP(killi.singlets.annot.c3, dims = 1:30)

pdf(paste0(Sys.Date(),"_Killi_Brain_Cohort3_UMAP_Singlets_labelled_UNION.pdf"), width = 6, height = 5)
DimPlot(killi.singlets.annot.c3, reduction = "umap", group.by = "DoubletCall", raster = T)
dev.off()

pdf(paste0(Sys.Date(),"_Killi_Brain_Cohort3_Key_Marker_Gene_expression_plots.pdf"), width = 3.5, height = 3)
FeaturePlot(killi.singlets.annot.c3, features = "olig2", raster = T)
FeaturePlot(killi.singlets.annot.c3, features = "olig1", raster = T)
FeaturePlot(killi.singlets.annot.c3, features = "mpz"  , raster = T)
FeaturePlot(killi.singlets.annot.c3, features = "csf1r", raster = T)
dev.off()

# DimPlot(killi.singlets.annot.c3, reduction = "umap", group.by = "DoubletFinder")
# DimPlot(killi.singlets.annot.c3, reduction = "umap", group.by = "scds_hybrid")

# save annotated object
save(killi.singlets.annot.c3, file = paste0(Sys.Date(),"_Killifish_Brain_Cohort3_Seurat_object_with_AnnotatedDoublets.RData"))


### extract/subset only singlets
# save data for singlets df
killi.singlets.c3   <- subset(killi.singlets.annot.c3, subset = DoubletCall %in% "Singlet")  # only keep singlets
killi.singlets.c3
# An object of class Seurat
# 52641 features across 66352 samples within 2 assays
# Active assay: SCT (24807 features, 3000 variable features)
#  1 other assay present: RNA
#  2 dimensional reductions calculated: pca, umap

pdf(paste0(Sys.Date(),"_Killi_Brain_Cohort3_UMAP_Singlets_ONLY_UNION.pdf"), width = 6, height = 5)
DimPlot(killi.singlets.c3, reduction = "umap", group.by = "DoubletCall", raster = T)
dev.off()

table(killi.singlets.c3@meta.data$Group)
# GRZ_M_F GRZ_M_M GRZ_O_F GRZ_O_M GRZ_Y_F GRZ_Y_M ZMZ_G_F ZMZ_G_M ZMZ_M_F ZMZ_M_M  ZMZ_O_F ZMZ_O_M ZMZ_Y_F ZMZ_Y_M
#    5990    5451    5225    5551    5988    6596    1498    1734    5143    7205    4468    2121    4147    5235


# save filtered/annotated object
save(killi.singlets.c3, file = paste0(Sys.Date(),"_Killifish_Brain_Cohort3_Seurat_object_SINGLETS_ONLY.RData"))
################################################################################################


#######################
sink(file = paste(Sys.Date(),"_Killifish_Brain_Cohort3_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()


