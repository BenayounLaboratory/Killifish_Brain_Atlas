setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Preprocessing/DecontX/Cohort_4/')
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
# 2022-12-24
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
cts.GRZ_Y_F_4     <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-11-22_snRNAseq_brain_aging_set_4/CellRanger/GRZ_5w_F_4_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.GRZ_Y_M_4     <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-11-22_snRNAseq_brain_aging_set_4/CellRanger/GRZ_5w_M_4_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.GRZ_M_F_4     <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-11-22_snRNAseq_brain_aging_set_4/CellRanger/GRZ_10w_F_4_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.GRZ_M_M_4     <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-11-22_snRNAseq_brain_aging_set_4/CellRanger/GRZ_10w_M_4_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.GRZ_O_F_4     <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-11-22_snRNAseq_brain_aging_set_4/CellRanger/GRZ_15w_F_4_FishTEDB_NR/outs/filtered_feature_bc_matrix/')
cts.GRZ_O_M_4     <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-11-22_snRNAseq_brain_aging_set_4/CellRanger/GRZ_15w_M_4_FishTEDB_NR/outs/filtered_feature_bc_matrix/')

cts_raw.GRZ_Y_F_4  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-11-22_snRNAseq_brain_aging_set_4/CellRanger/GRZ_5w_F_4_FishTEDB_NR/outs/raw_feature_bc_matrix/')
cts_raw.GRZ_Y_M_4  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-11-22_snRNAseq_brain_aging_set_4/CellRanger/GRZ_5w_M_4_FishTEDB_NR/outs/raw_feature_bc_matrix/')
cts_raw.GRZ_M_F_4  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-11-22_snRNAseq_brain_aging_set_4/CellRanger/GRZ_10w_F_4_FishTEDB_NR/outs/raw_feature_bc_matrix/')
cts_raw.GRZ_M_M_4  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-11-22_snRNAseq_brain_aging_set_4/CellRanger/GRZ_10w_M_4_FishTEDB_NR/outs/raw_feature_bc_matrix/')
cts_raw.GRZ_O_F_4  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-11-22_snRNAseq_brain_aging_set_4/CellRanger/GRZ_15w_F_4_FishTEDB_NR/outs/raw_feature_bc_matrix/')
cts_raw.GRZ_O_M_4  <- Read10X('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-11-22_snRNAseq_brain_aging_set_4/CellRanger/GRZ_15w_M_4_FishTEDB_NR/outs/raw_feature_bc_matrix/')


# Create SingleCellExperiment objects
sce.GRZ_Y_F_4     <- SingleCellExperiment(list(counts = cts.GRZ_Y_F_4 ))
sce.GRZ_Y_M_4     <- SingleCellExperiment(list(counts = cts.GRZ_Y_M_4 ))
sce.GRZ_M_F_4     <- SingleCellExperiment(list(counts = cts.GRZ_M_F_4 ))
sce.GRZ_M_M_4     <- SingleCellExperiment(list(counts = cts.GRZ_M_M_4 ))
sce.GRZ_O_F_4     <- SingleCellExperiment(list(counts = cts.GRZ_O_F_4 ))
sce.GRZ_O_M_4     <- SingleCellExperiment(list(counts = cts.GRZ_O_M_4 ))

sce_raw.GRZ_Y_F_4  <- SingleCellExperiment(list(counts = cts_raw.GRZ_Y_F_4 ))
sce_raw.GRZ_Y_M_4  <- SingleCellExperiment(list(counts = cts_raw.GRZ_Y_M_4 ))
sce_raw.GRZ_M_F_4  <- SingleCellExperiment(list(counts = cts_raw.GRZ_M_F_4 ))
sce_raw.GRZ_M_M_4  <- SingleCellExperiment(list(counts = cts_raw.GRZ_M_M_4 ))
sce_raw.GRZ_O_F_4  <- SingleCellExperiment(list(counts = cts_raw.GRZ_O_F_4 ))
sce_raw.GRZ_O_M_4  <- SingleCellExperiment(list(counts = cts_raw.GRZ_O_M_4 ))


# Run decontX
sce.GRZ_Y_F_4     <- decontX(sce.GRZ_Y_F_4, background = sce_raw.GRZ_Y_F_4 )
sce.GRZ_Y_M_4     <- decontX(sce.GRZ_Y_M_4, background = sce_raw.GRZ_Y_M_4 )
sce.GRZ_M_F_4     <- decontX(sce.GRZ_M_F_4, background = sce_raw.GRZ_M_F_4 )
sce.GRZ_M_M_4     <- decontX(sce.GRZ_M_M_4, background = sce_raw.GRZ_M_M_4 )
sce.GRZ_O_F_4     <- decontX(sce.GRZ_O_F_4, background = sce_raw.GRZ_O_F_4 )
sce.GRZ_O_M_4     <- decontX(sce.GRZ_O_M_4, background = sce_raw.GRZ_O_M_4 )


# get seurat objects
seurat.GRZ_Y_F_4  <- CreateSeuratObject( round(decontXcounts(sce.GRZ_Y_F_4)) )
seurat.GRZ_Y_M_4  <- CreateSeuratObject( round(decontXcounts(sce.GRZ_Y_M_4)) )
seurat.GRZ_M_F_4  <- CreateSeuratObject( round(decontXcounts(sce.GRZ_M_F_4)) )
seurat.GRZ_M_M_4  <- CreateSeuratObject( round(decontXcounts(sce.GRZ_M_M_4)) )
seurat.GRZ_O_F_4  <- CreateSeuratObject( round(decontXcounts(sce.GRZ_O_F_4)) )
seurat.GRZ_O_M_4  <- CreateSeuratObject( round(decontXcounts(sce.GRZ_O_M_4)) )

seurat.GRZ_Y_F_4  <- AddMetaData(object = seurat.GRZ_Y_F_4 , colData(sce.GRZ_Y_F_4)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.GRZ_Y_M_4  <- AddMetaData(object = seurat.GRZ_Y_M_4 , colData(sce.GRZ_Y_M_4)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.GRZ_M_F_4  <- AddMetaData(object = seurat.GRZ_M_F_4 , colData(sce.GRZ_M_F_4)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.GRZ_M_M_4  <- AddMetaData(object = seurat.GRZ_M_M_4 , colData(sce.GRZ_M_M_4)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.GRZ_O_F_4  <- AddMetaData(object = seurat.GRZ_O_F_4 , colData(sce.GRZ_O_F_4)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.GRZ_O_M_4  <- AddMetaData(object = seurat.GRZ_O_M_4 , colData(sce.GRZ_O_M_4)$decontX_contamination  , col.name = "decontX_contamination"  )



# Merge objects for the cohort
brain.c4.combined <- merge(seurat.GRZ_Y_F_4,
                           y =  c(seurat.GRZ_Y_M_4,
                                  seurat.GRZ_M_F_4,
                                  seurat.GRZ_M_M_4,
                                  seurat.GRZ_O_F_4,
                                  seurat.GRZ_O_M_4),
                           add.cell.ids = c("GRZ_Y_F_4"  ,
                                            "GRZ_Y_M_4"  ,
                                            "GRZ_M_F_4" ,
                                            "GRZ_M_M_4" ,
                                            "GRZ_O_F_4" ,
                                            "GRZ_O_M_4"),
                           project = "10x_Killi_Brain_Aging")
brain.c4.combined
# An object of class Seurat
# 27834 features across 29142 samples within 1 assay
# Active assay: RNA (27834 features, 0 variable features)

# clean memory
rm(cts.GRZ_Y_F_4     ,cts.GRZ_Y_M_4     ,cts.GRZ_M_F_4     ,cts.GRZ_M_M_4     ,cts.GRZ_O_F_4     ,cts.GRZ_O_M_4     ,cts_raw.GRZ_Y_F_4 ,cts_raw.GRZ_Y_M_4 ,cts_raw.GRZ_M_F_4 ,cts_raw.GRZ_M_M_4 ,cts_raw.GRZ_O_F_4 ,cts_raw.GRZ_O_M_4 ,sce.GRZ_Y_F_4    ,sce.GRZ_Y_M_4    ,sce.GRZ_M_F_4    ,sce.GRZ_M_M_4    ,sce.GRZ_O_F_4    ,sce.GRZ_O_M_4    ,sce_raw.GRZ_Y_F_4,sce_raw.GRZ_Y_M_4,sce_raw.GRZ_M_F_4,sce_raw.GRZ_M_M_4,sce_raw.GRZ_O_F_4,sce_raw.GRZ_O_M_4,seurat.GRZ_Y_F_4,seurat.GRZ_Y_M_4,seurat.GRZ_M_F_4,seurat.GRZ_M_M_4,seurat.GRZ_O_F_4,seurat.GRZ_O_M_4)
################################################################################################


################################################################################################
#### 1b. Add key metadata to Seurat object

# create Group label
my.GRZ_Y_F_4   <- grep("GRZ_Y_F_4"  , colnames(brain.c4.combined@assays$RNA))
my.GRZ_Y_M_4   <- grep("GRZ_Y_M_4"  , colnames(brain.c4.combined@assays$RNA))
my.GRZ_M_F_4   <- grep("GRZ_M_F_4"  , colnames(brain.c4.combined@assays$RNA))
my.GRZ_M_M_4   <- grep("GRZ_M_M_4"  , colnames(brain.c4.combined@assays$RNA))
my.GRZ_O_F_4   <- grep("GRZ_O_F_4"  , colnames(brain.c4.combined@assays$RNA))
my.GRZ_O_M_4   <- grep("GRZ_O_M_4"  , colnames(brain.c4.combined@assays$RNA))

##### (even though we initially labeled 5/15w, Ari used 6/16w - correct in meta data attribution)
Group <- rep("NA", length(colnames(brain.c4.combined@assays$RNA)))
Group[ my.GRZ_Y_F_4  ]   <- "GRZ_Y_F"
Group[ my.GRZ_Y_M_4  ]   <- "GRZ_Y_M"
Group[ my.GRZ_M_F_4  ]   <- "GRZ_M_F"
Group[ my.GRZ_M_M_4  ]   <- "GRZ_M_M"
Group[ my.GRZ_O_F_4  ]   <- "GRZ_O_F"
Group[ my.GRZ_O_M_4  ]   <- "GRZ_O_M"
Group <- data.frame(Group)
rownames(Group) <- colnames(brain.c4.combined@assays$RNA)

#####
Sex <- rep("NA", length(colnames(brain.c4.combined@assays$RNA)))
Sex[ my.GRZ_Y_F_4 ]   <- "F"  
Sex[ my.GRZ_Y_M_4 ]   <- "M" 
Sex[ my.GRZ_M_F_4 ]   <- "F" 
Sex[ my.GRZ_M_M_4 ]   <- "M"
Sex[ my.GRZ_O_F_4 ]   <- "F" 
Sex[ my.GRZ_O_M_4 ]   <- "M"  
Sex <- data.frame(Sex)
rownames(Sex) <- colnames(brain.c4.combined@assays$RNA)

##### Y, M, O, G
Age.gp <- rep("NA", length(colnames(brain.c4.combined@assays$RNA)))
Age.gp[ my.GRZ_Y_F_4 ]   <- "Y"
Age.gp[ my.GRZ_Y_M_4 ]   <- "Y"
Age.gp[ my.GRZ_M_F_4 ]   <- "M"
Age.gp[ my.GRZ_M_M_4 ]   <- "M"
Age.gp[ my.GRZ_O_F_4 ]   <- "O"
Age.gp[ my.GRZ_O_M_4 ]   <- "O"
Age.gp <- data.frame(Age.gp)
rownames(Age.gp) <- colnames(brain.c4.combined@assays$RNA)

##### (average of the pool) 
Age.w <- rep("NA", length(colnames(brain.c4.combined@assays$RNA)))
Age.w[ my.GRZ_Y_F_4 ]   <-   6.1  # 
Age.w[ my.GRZ_Y_M_4 ]   <-   6.1  # 
Age.w[ my.GRZ_M_F_4 ]   <-  10.6  # 
Age.w[ my.GRZ_M_M_4 ]   <-   9.4  # 
Age.w[ my.GRZ_O_F_4 ]   <-  15.7  # 
Age.w[ my.GRZ_O_M_4 ]   <-  15.4  # 
Age.w <- data.frame(Age.w)
rownames(Age.w) <- colnames(brain.c4.combined@assays$RNA)

#####
Strain <- rep("NA", length(colnames(brain.c4.combined@assays$RNA)))
Strain[ my.GRZ_Y_F_4 ]   <- "GRZ"
Strain[ my.GRZ_Y_M_4 ]   <- "GRZ"
Strain[ my.GRZ_M_F_4 ]   <- "GRZ"
Strain[ my.GRZ_M_M_4 ]   <- "GRZ"
Strain[ my.GRZ_O_F_4 ]   <- "GRZ"
Strain[ my.GRZ_O_M_4 ]   <- "GRZ"
Strain <- data.frame(Strain)
rownames(Strain) <- colnames(brain.c4.combined@assays$RNA)


#####
Batch <- rep("Set_4", length(colnames(brain.c4.combined@assays$RNA)))
Batch <- data.frame(Batch)
rownames(Batch) <- colnames(brain.c4.combined@assays$RNA)

# update Seurat with metadata
brain.c4.combined <- AddMetaData(object = brain.c4.combined, metadata = as.vector(Group)    , col.name = "Group"       )
brain.c4.combined <- AddMetaData(object = brain.c4.combined, metadata = as.vector(Sex)      , col.name = "Sex"         )
brain.c4.combined <- AddMetaData(object = brain.c4.combined, metadata = as.vector(Age.w)    , col.name = "Age_weeks"   )
brain.c4.combined <- AddMetaData(object = brain.c4.combined, metadata = as.vector(Age.gp)   , col.name = "Age_Group"   )
brain.c4.combined <- AddMetaData(object = brain.c4.combined, metadata = as.vector(Strain)   , col.name = "Strain"      )
brain.c4.combined <- AddMetaData(object = brain.c4.combined, metadata = as.vector(Batch)    , col.name = "Batch"       )
################################################################################################


################################################################################################
#### 2. Basic QC and filtering with Seurat

### No filtering on genes at this stage - only after all cohorts merged for fairness
brain.c4.combined <- SetIdent(brain.c4.combined, value = "Group")

# DecontX contamination levels for filtration
pdf(paste(Sys.Date(),"Killifish_Brain_Cohort4_violinPlots_QC_DecontX.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = brain.c4.combined, features = c("decontX_contamination"), pt.size = 0)
dev.off()

# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
brain.c4.combined[["percent.mito"]] <- PercentageFeatureSet(brain.c4.combined, pattern = "^MT-")

pdf(paste(Sys.Date(),"Killifish_Brain_Cohort4_violinPlots_QC_gene_UMI_mito.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = brain.c4.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(brain.c4.combined, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(brain.c4.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste(Sys.Date(),"Killifish_Brain_Cohort4_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
plot1 + plot2
dev.off()

# filter dead/low Q cells
brain.c4.combined <- subset(brain.c4.combined, subset = nFeature_RNA > 250 & nFeature_RNA < 5000 & percent.mito < 10 & nCount_RNA < 25000 & decontX_contamination < 0.25 )
brain.c4.combined
# An object of class Seurat
# 27834 features across 24874 samples within 1 assay
# Active assay: RNA (27834 features, 0 variable features)


### Check data after cell filtering
head(brain.c4.combined@meta.data)

table(brain.c4.combined@meta.data$Group)
# GRZ_M_F GRZ_M_M GRZ_O_F GRZ_O_M GRZ_Y_F GRZ_Y_M
#    4643    4241    3924    4491    4043    3532

#### Normalize the data for doublet analysis, etc
# global-scaling normalization method 'LogNormalize' normalizes gene expression measurements for each cell 
# by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
brain.c4.combined <- NormalizeData(object = brain.c4.combined, normalization.method = "LogNormalize",  scale.factor = 10000)
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
brain.c4.combined <- CellCycleScoring(object = brain.c4.combined, s.features = s.genes.k, g2m.features = g2m.genes.k, set.ident = TRUE)

# write predictions to file
write.table(brain.c4.combined@meta.data, file = paste0(Sys.Date(),"_Killi_Brain_Cohort4_CellCycle_predictions.txt"), sep = "\t", quote = F)
################################################################################################


################################################################################################
#### 3. Find and remove doublets using doublet finder & scds workflow

# pre normalize the data for doubletFinder
brain.c4.combined <- SCTransform(object = brain.c4.combined, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mito", "Phase"))
save(brain.c4.combined, file = paste0(Sys.Date(),"_Killi_Brain_Cohort4_Seurat_object_postSCT.RData"))

# Run first pass analysis just for doublet identification (not final clustering)
brain.c4.combined <- RunPCA(brain.c4.combined, npcs = 30)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_Killi_Brain_Cohort4_ElbowPlot.pdf"))
ElbowPlot(brain.c4.combined, ndims = 30)
dev.off()

# run dimensionality reduction
# Keep all PCs here, we'll do the clean clustering analysis on the merged object across all cohorts
brain.c4.combined <- RunUMAP(brain.c4.combined, dims = 1:30)
brain.c4.combined <- FindNeighbors(brain.c4.combined, dims = 1:30)
brain.c4.combined <- FindClusters(object = brain.c4.combined)

#### need to split by 10x sample to make sure to identify real doublets
# will run on one object at a time
cohort.4.list <- SplitObject(brain.c4.combined, split.by = "Group")

## Assume doublet rate based on 10x information (add a 15% fudge factor due to nuclei being more sticky)
pred.dblt.rate <- 1.15 * predict(pred_dblt_lm, data.frame("cell_number" = unlist(lapply(cohort.4.list, ncol))))/100

## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& A. Run DoubletFinder
# loop over samples
for (i in 1:length(cohort.4.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list_killi <- paramSweep_v3(cohort.4.list[[i]], PCs = 1:30, sct = TRUE, num.cores	 = 4)
  sweep.stats_killi    <- summarizeSweep(sweep.res.list_killi, GT = FALSE)
  bcmvn_killi          <- find.pK(sweep.stats_killi)
  
  # need some R gymnastics since the Pk is stored as a factor for some reason
  # to get the pK number, need to first convert to character and THEN to numeric
  # numeric first yield row number
  pk.killi <- as.numeric(as.character(bcmvn_killi[as.numeric(bcmvn_killi$pK[bcmvn_killi$BCmetric == max(bcmvn_killi$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(cohort.4.list[[i]]@meta.data$seurat_clusters)             ## ex: annotations
  nExp_poi       <- round((pred.dblt.rate[i]) *length(cohort.4.list[[i]]@meta.data$Group))      ## Assume doublets based on nuclei isolation protocol performance
  
  ## Run DoubletFinder with varying classification stringencies
  cohort.4.list[[i]] <- doubletFinder_v3(cohort.4.list[[i]], PCs = 1:30, pN = 0.25, pK = pk.killi, nExp = nExp_poi,     reuse.pANN = FALSE, sct = T)
  
  # get classification name
  my.DF.res.col <- colnames(cohort.4.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(cohort.4.list[[i]]@meta.data))]
  
  # rename column to enable subsetting
  colnames(cohort.4.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(cohort.4.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# run UMAP plots
for (i in 1:length(cohort.4.list)) {
  pdf(paste(Sys.Date(),"Killifish_Tissue",names(cohort.4.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(cohort.4.list[[i]], reduction = "umap", group.by = "DoubletFinder"), raster = T)
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
killi.singlets.annot.c4 <- merge(cohort.4.list[[1]],
                                 y = c(cohort.4.list[[ 2]],
                                       cohort.4.list[[ 3]],
                                       cohort.4.list[[ 4]],
                                       cohort.4.list[[ 5]],
                                       cohort.4.list[[ 6]]),
                                 project = "Killi_Brain_Cohort4")
killi.singlets.annot.c4

 
# remove pANN columns that are 10xGenomics library lane specific
killi.singlets.annot.c4@meta.data <- killi.singlets.annot.c4@meta.data[,-grep("pANN",colnames(killi.singlets.annot.c4@meta.data))]


## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& B. Run scds:single cell doublet scoring (hybrid method)
# cxds is based on co-expression of gene pairs and works with absence/presence calls only, 
# bcds uses the full count information and a binary classification approach using artificially generated doublets. 
# cxds_bcds_hybrid combines both approaches

# create scds working object - convert list to SingleCellExperiment
cohort.4.list.scds        <- lapply(cohort.4.list, as.SingleCellExperiment)

# loop over sample
for (i in 1:length(cohort.4.list.scds)) {
  
  # Annotate doublets using co-expression based doublet scoring:
  cohort.4.list.scds[[i]] <- cxds_bcds_hybrid(cohort.4.list.scds[[i]])
  
  # predicted doublet rate
  n.db <- round((pred.dblt.rate[i])*ncol(cohort.4.list.scds[[i]]))                         ## Assume doublets based on nuclei isolation protocol performance
  
  # sort prediction, get top n.db cells
  srt.db.score <- sort(cohort.4.list.scds[[i]]$hybrid_score, index.return = T, decreasing = T)
  cohort.4.list.scds[[i]]$scds <- "Singlet"
  cohort.4.list.scds[[i]]$scds[srt.db.score$ix[1:n.db]] <- "Doublet"
  
}

# run UMAP plots
for (i in 1:length(cohort.4.list.scds)) {
  
  p <- plotReducedDim(cohort.4.list.scds[[i]], dimred = "UMAP", colour_by = "scds")
  
  pdf(paste(Sys.Date(),"Killifish_",names(cohort.4.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to doubletFinder annotated Seurat object
killi.singlets.annot.c4@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(cohort.4.list.scds)) {
  
  # for each object compare and move doublet annotations over
  killi.singlets.annot.c4@meta.data[colnames(cohort.4.list.scds[[i]]), ]$scds_hybrid <- cohort.4.list.scds[[i]]$scds
  
}

# free some memory
rm(cohort.4.list, cohort.4.list.scds, brain.c4.combined)

## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& C. Merge and summarize doublet findings

table(killi.singlets.annot.c4@meta.data$DoubletFinder, killi.singlets.annot.c4@meta.data$scds_hybrid)
#             Doublet Singlet
#     Doublet     199     757
#     Singlet     757   23161
  
  
# Union (more conservative)
killi.singlets.annot.c4@meta.data$DoubletCall <- ifelse( bitOr(killi.singlets.annot.c4@meta.data$DoubletFinder == "Doublet", killi.singlets.annot.c4@meta.data$scds_hybrid == "Doublet") > 0, 
                                                           "Doublet", "Singlet")
table(killi.singlets.annot.c4@meta.data$DoubletCall)
# Doublet Singlet 
#    1713   23161


# re-run dimensionality reduction for plotting purposes
killi.singlets.annot.c4 <- SCTransform(object = killi.singlets.annot.c4, vars.to.regress =  c("nFeature_RNA", "nCount_RNA", "percent.mito"))
killi.singlets.annot.c4 <- RunPCA(killi.singlets.annot.c4, npcs = 30)
killi.singlets.annot.c4 <- RunUMAP(killi.singlets.annot.c4, dims = 1:30)

pdf(paste0(Sys.Date(),"_Killi_Brain_Cohort4_UMAP_Singlets_labelled_UNION.pdf"), width = 6, height = 5)
DimPlot(killi.singlets.annot.c4, reduction = "umap", group.by = "DoubletCall", raster = T)
dev.off()

pdf(paste0(Sys.Date(),"_Killi_Brain_Cohort4_Key_Marker_Gene_expression_plots.pdf"), width = 3.5, height = 3)
FeaturePlot(killi.singlets.annot.c4, features = "olig2", raster = T)
FeaturePlot(killi.singlets.annot.c4, features = "olig1", raster = T)
FeaturePlot(killi.singlets.annot.c4, features = "mpz"  , raster = T)
FeaturePlot(killi.singlets.annot.c4, features = "csf1r", raster = T)
dev.off()

# DimPlot(killi.singlets.annot.c4, reduction = "umap", group.by = "DoubletFinder")
# DimPlot(killi.singlets.annot.c4, reduction = "umap", group.by = "scds_hybrid")

# save annotated object
save(killi.singlets.annot.c4, file = paste0(Sys.Date(),"_Killifish_Brain_Cohort4_Seurat_object_with_AnnotatedDoublets.RData"))


### extract/subset only singlets
# save data for singlets df
killi.singlets.c4   <- subset(killi.singlets.annot.c4, subset = DoubletCall %in% "Singlet")  # only keep singlets
killi.singlets.c4
# An object of class Seurat
# 51665 features across 23161 samples within 2 assays
# Active assay: SCT (23831 features, 3000 variable features)
#  1 other assay present: RNA
#  2 dimensional reductions calculated: pca, umap

pdf(paste0(Sys.Date(),"_Killi_Brain_Cohort4_UMAP_Singlets_ONLY_UNION.pdf"), width = 6, height = 5)
DimPlot(killi.singlets.c4, reduction = "umap", group.by = "DoubletCall", raster = T)
dev.off()

table(killi.singlets.c4@meta.data$Group)
# GRZ_M_F GRZ_M_M GRZ_O_F GRZ_O_M GRZ_Y_F GRZ_Y_M
#    4289    3950    3665    4142    3786    3329

# save filtered/annotated object
save(killi.singlets.c4, file = paste0(Sys.Date(),"_Killifish_Brain_Cohort4_Seurat_object_SINGLETS_ONLY.RData"))
################################################################################################


#######################
sink(file = paste(Sys.Date(),"_Killifish_Brain_Cohort4_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()




