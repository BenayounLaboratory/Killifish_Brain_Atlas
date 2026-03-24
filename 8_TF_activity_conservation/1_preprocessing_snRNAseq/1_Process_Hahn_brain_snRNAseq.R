setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Mouse_Datasets_for_Comparison/Processed_Objects/Hahn_GSE212576')
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
# 2024-03-25
# Process cellranger data from Hahn et al, 2023 using our pipeline since processed object was not availble
# cellranger data was provided in GEO
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
################################################################################################


# ################################################################################################
# #### 1a. Read data from CellRanger, perform background removal
# 
# # Calculate and clean the contribution of ambient RNA with DecontX
# # read 10x libraries cell ranger gene barcode matrices for DecontX
# # https://github.com/satijalab/seurat/issues/4096
# cts.HIP_Y_1     <- ReadMtx(mtx = './GSM6537880/GSM6537880_HIP_Y_1_matrix.mtx.gz', cells = './GSM6537880/GSM6537880_HIP_Y_1_barcodes.tsv.gz', features = './GSM6537880/GSM6537880_HIP_Y_1_features.tsv.gz')
# cts.HIP_Y_2     <- ReadMtx(mtx = './GSM6537881/GSM6537881_HIP_Y_2_matrix.mtx.gz', cells = './GSM6537881/GSM6537881_HIP_Y_2_barcodes.tsv.gz', features = './GSM6537881/GSM6537881_HIP_Y_2_features.tsv.gz')
# cts.HIP_O_1     <- ReadMtx(mtx = './GSM6537882/GSM6537882_HIP_O_1_matrix.mtx.gz', cells = './GSM6537882/GSM6537882_HIP_O_1_barcodes.tsv.gz', features = './GSM6537882/GSM6537882_HIP_O_1_features.tsv.gz')
# cts.HIP_O_2     <- ReadMtx(mtx = './GSM6537883/GSM6537883_HIP_O_2_matrix.mtx.gz', cells = './GSM6537883/GSM6537883_HIP_O_2_barcodes.tsv.gz', features = './GSM6537883/GSM6537883_HIP_O_2_features.tsv.gz')
# cts.CP_Y_1      <- ReadMtx(mtx = './GSM6537884/GSM6537884_CP_Y_1_matrix.mtx.gz', cells = './GSM6537884/GSM6537884_CP_Y_1_barcodes.tsv.gz', features = './GSM6537884/GSM6537884_CP_Y_1_features.tsv.gz')
# cts.CP_Y_2      <- ReadMtx(mtx = './GSM6537885/GSM6537885_CP_Y_2_matrix.mtx.gz', cells = './GSM6537885/GSM6537885_CP_Y_2_barcodes.tsv.gz', features = './GSM6537885/GSM6537885_CP_Y_2_features.tsv.gz')
# cts.CP_O_1      <- ReadMtx(mtx = './GSM6537886/GSM6537886_CP_O_1_matrix.mtx.gz', cells = './GSM6537886/GSM6537886_CP_O_1_barcodes.tsv.gz', features = './GSM6537886/GSM6537886_CP_O_1_features.tsv.gz')
# cts.CP_O_2      <- ReadMtx(mtx = './GSM6537887/GSM6537887_CP_O_2_matrix.mtx.gz', cells = './GSM6537887/GSM6537887_CP_O_2_barcodes.tsv.gz', features = './GSM6537887/GSM6537887_CP_O_2_features.tsv.gz')
# 
# # Create SingleCellExperiment objects
# sce.HIP_Y_1      <- SingleCellExperiment(list(counts = cts.HIP_Y_1      ))
# sce.HIP_Y_2      <- SingleCellExperiment(list(counts = cts.HIP_Y_2      ))
# sce.HIP_O_1      <- SingleCellExperiment(list(counts = cts.HIP_O_1      ))
# sce.HIP_O_2      <- SingleCellExperiment(list(counts = cts.HIP_O_2      ))
# sce.CP_Y_1       <- SingleCellExperiment(list(counts = cts.CP_Y_1       ))
# sce.CP_Y_2       <- SingleCellExperiment(list(counts = cts.CP_Y_2       ))
# sce.CP_O_1       <- SingleCellExperiment(list(counts = cts.CP_O_1       ))
# sce.CP_O_2       <- SingleCellExperiment(list(counts = cts.CP_O_2       ))
# 
# # Run decontX (no available background)
# sce.HIP_Y_1   <- decontX(sce.HIP_Y_1  )
# sce.HIP_Y_2   <- decontX(sce.HIP_Y_2  )
# sce.HIP_O_1   <- decontX(sce.HIP_O_1  )
# sce.HIP_O_2   <- decontX(sce.HIP_O_2  )
# sce.CP_Y_1    <- decontX(sce.CP_Y_1   )
# sce.CP_Y_2    <- decontX(sce.CP_Y_2   )
# sce.CP_O_1    <- decontX(sce.CP_O_1   )
# sce.CP_O_2    <- decontX(sce.CP_O_2   )
# 
# # get seurat objects
# seurat.HIP_Y_1   <- CreateSeuratObject( round(decontXcounts(sce.HIP_Y_1  )) )
# seurat.HIP_Y_2   <- CreateSeuratObject( round(decontXcounts(sce.HIP_Y_2  )) )
# seurat.HIP_O_1   <- CreateSeuratObject( round(decontXcounts(sce.HIP_O_1  )) )
# seurat.HIP_O_2   <- CreateSeuratObject( round(decontXcounts(sce.HIP_O_2  )) )
# seurat.CP_Y_1    <- CreateSeuratObject( round(decontXcounts(sce.CP_Y_1   )) )
# seurat.CP_Y_2    <- CreateSeuratObject( round(decontXcounts(sce.CP_Y_2   )) )
# seurat.CP_O_1    <- CreateSeuratObject( round(decontXcounts(sce.CP_O_1   )) )
# seurat.CP_O_2    <- CreateSeuratObject( round(decontXcounts(sce.CP_O_2   )) )
# 
# seurat.HIP_Y_1   <- AddMetaData(object = seurat.HIP_Y_1  , colData(sce.HIP_Y_1  )$decontX_contamination  , col.name = "decontX_contamination"  )
# seurat.HIP_Y_2   <- AddMetaData(object = seurat.HIP_Y_2  , colData(sce.HIP_Y_2  )$decontX_contamination  , col.name = "decontX_contamination"  )
# seurat.HIP_O_1   <- AddMetaData(object = seurat.HIP_O_1  , colData(sce.HIP_O_1  )$decontX_contamination  , col.name = "decontX_contamination"  )
# seurat.HIP_O_2   <- AddMetaData(object = seurat.HIP_O_2  , colData(sce.HIP_O_2  )$decontX_contamination  , col.name = "decontX_contamination"  )
# seurat.CP_Y_1    <- AddMetaData(object = seurat.CP_Y_1   , colData(sce.CP_Y_1   )$decontX_contamination  , col.name = "decontX_contamination"  )
# seurat.CP_Y_2    <- AddMetaData(object = seurat.CP_Y_2   , colData(sce.CP_Y_2   )$decontX_contamination  , col.name = "decontX_contamination"  )
# seurat.CP_O_1    <- AddMetaData(object = seurat.CP_O_1   , colData(sce.CP_O_1   )$decontX_contamination  , col.name = "decontX_contamination"  )
# seurat.CP_O_2    <- AddMetaData(object = seurat.CP_O_2   , colData(sce.CP_O_2   )$decontX_contamination  , col.name = "decontX_contamination"  )
# 
# # Merge objects for the cohort
# hahn.combined <- merge(seurat.HIP_Y_1,
#                            y =  c(seurat.HIP_Y_2 ,
#                                   seurat.HIP_O_1 ,
#                                   seurat.HIP_O_2 ,
#                                   seurat.CP_Y_1  ,
#                                   seurat.CP_Y_2  ,
#                                   seurat.CP_O_1  ,
#                                   seurat.CP_O_2  ),
#                            add.cell.ids = c("HIP_Y_1"  ,
#                                             "HIP_Y_2"  ,
#                                             "HIP_O_1"  ,
#                                             "HIP_O_2"  ,
#                                             "CP_Y_1"   ,
#                                             "CP_Y_2"   ,
#                                             "CP_O_1"   ,
#                                             "CP_O_2"   ),
#                            project = "Hahn_Brain_Aging")
# hahn.combined
# # An object of class Seurat 
# # 32285 features across 109826 samples within 1 assay 
# # Active assay: RNA (32285 features, 0 variable features)
# 
# #clean memory
# rm(   cts.HIP_Y_1 ,   cts.HIP_Y_2 ,   cts.HIP_O_1 ,   cts.HIP_O_2 ,   cts.CP_Y_1  ,   cts.CP_Y_2  ,   cts.CP_O_1  ,   cts.CP_O_2  ,   sce.HIP_Y_1  ,   sce.HIP_Y_2  ,   sce.HIP_O_1  ,   sce.HIP_O_2  ,   sce.CP_Y_1   ,   sce.CP_Y_2   ,   sce.CP_O_1   ,   sce.CP_O_2   ,   seurat.HIP_Y_1 ,   seurat.HIP_Y_2 ,   seurat.HIP_O_1 ,   seurat.HIP_O_2 ,   seurat.CP_Y_1  ,   seurat.CP_Y_2  ,   seurat.CP_O_1  ,   seurat.CP_O_2   )
# ################################################################################################
# 
# 
# ################################################################################################
# #### 1b. Add key metadata to Seurat object
# 
# # create Sample label
# my.HIP_Y_1   <- grep("HIP_Y_1"  , colnames(hahn.combined@assays$RNA))
# my.HIP_Y_2   <- grep("HIP_Y_2"  , colnames(hahn.combined@assays$RNA))
# my.HIP_O_1   <- grep("HIP_O_1"  , colnames(hahn.combined@assays$RNA))
# my.HIP_O_2   <- grep("HIP_O_2"  , colnames(hahn.combined@assays$RNA))
# my.CP_Y_1    <- grep("CP_Y_1"  , colnames(hahn.combined@assays$RNA))
# my.CP_Y_2    <- grep("CP_Y_2"  , colnames(hahn.combined@assays$RNA))
# my.CP_O_1    <- grep("CP_O_1"  , colnames(hahn.combined@assays$RNA))
# my.CP_O_2    <- grep("CP_O_2"  , colnames(hahn.combined@assays$RNA))
# 
# #####
# SampleID <- rep("NA", length(colnames(hahn.combined@assays$RNA)))
# SampleID[ my.HIP_Y_1   ]   <- "HIP_Y_1"
# SampleID[ my.HIP_Y_2   ]   <- "HIP_Y_2"
# SampleID[ my.HIP_O_1   ]   <- "HIP_O_1"
# SampleID[ my.HIP_O_2   ]   <- "HIP_O_2"
# SampleID[ my.CP_Y_1    ]   <- "CP_Y_1"
# SampleID[ my.CP_Y_2    ]   <- "CP_Y_2"
# SampleID[ my.CP_O_1    ]   <- "CP_O_1"
# SampleID[ my.CP_O_2    ]   <- "CP_O_2"
# SampleID <- data.frame(SampleID)
# rownames(SampleID) <- colnames(hahn.combined@assays$RNA)
# 
# #####
# Sex <- rep("NA", length(colnames(hahn.combined@assays$RNA)))
# Sex[ my.HIP_Y_1  ]   <- "MF"  
# Sex[ my.HIP_Y_2  ]   <- "MF" 
# Sex[ my.HIP_O_1  ]   <- "MF" 
# Sex[ my.HIP_O_2  ]   <- "MF"
# Sex[ my.CP_Y_1   ]   <- "MF" 
# Sex[ my.CP_Y_2   ]   <- "MF"  
# Sex[ my.CP_O_1   ]   <- "MF"   
# Sex[ my.CP_O_2   ]   <- "MF" 
# Sex <- data.frame(Sex)
# rownames(Sex) <- colnames(hahn.combined@assays$RNA)
# 
# ##### Y, O
# Age.gp <- rep("NA", length(colnames(hahn.combined@assays$RNA)))
# Age.gp[ my.HIP_Y_1  ]   <- "Y"
# Age.gp[ my.HIP_Y_2  ]   <- "Y"
# Age.gp[ my.HIP_O_1  ]   <- "O"
# Age.gp[ my.HIP_O_2  ]   <- "O"
# Age.gp[ my.CP_Y_1   ]   <- "Y"
# Age.gp[ my.CP_Y_2   ]   <- "Y"
# Age.gp[ my.CP_O_1   ]   <- "O"
# Age.gp[ my.CP_O_2   ]   <- "O"
# Age.gp <- data.frame(Age.gp)
# rownames(Age.gp) <- colnames(hahn.combined@assays$RNA)
# 
# ##### (average of the pool) 
# Age.m <- rep("NA", length(colnames(hahn.combined@assays$RNA)))
# Age.m[ my.HIP_Y_1   ]   <-  3   #
# Age.m[ my.HIP_Y_2   ]   <-  3   #
# Age.m[ my.HIP_O_1   ]   <- 21   #
# Age.m[ my.HIP_O_2   ]   <- 21   #
# Age.m[ my.CP_Y_1    ]   <-  3    #
# Age.m[ my.CP_Y_2    ]   <-  3    #
# Age.m[ my.CP_O_1    ]   <- 21    #
# Age.m[ my.CP_O_2    ]   <- 21    #
# Age.m <- data.frame(Age.m)
# rownames(Age.m) <- colnames(hahn.combined@assays$RNA)
# 
# #####
# Tissue <- rep("NA", length(colnames(hahn.combined@assays$RNA)))
# Tissue[ my.HIP_Y_1  ]   <- "Hippocampus"
# Tissue[ my.HIP_Y_2  ]   <- "Hippocampus"
# Tissue[ my.HIP_O_1  ]   <- "Hippocampus"
# Tissue[ my.HIP_O_2  ]   <- "Hippocampus"
# Tissue[ my.CP_Y_1   ]   <- "Caudate_Putamen"
# Tissue[ my.CP_Y_2   ]   <- "Caudate_Putamen"
# Tissue[ my.CP_O_1   ]   <- "Caudate_Putamen"
# Tissue[ my.CP_O_2   ]   <- "Caudate_Putamen"
# Tissue <- data.frame(Tissue)
# rownames(Tissue) <- colnames(hahn.combined@assays$RNA)
# 
# 
# #####
# 
# # update Seurat with metadata
# hahn.combined <- AddMetaData(object = hahn.combined, metadata = as.vector(SampleID) , col.name = "SampleID"       )
# hahn.combined <- AddMetaData(object = hahn.combined, metadata = as.vector(Sex)      , col.name = "Sex"         )
# hahn.combined <- AddMetaData(object = hahn.combined, metadata = as.vector(Age.m)    , col.name = "Age_months"   )
# hahn.combined <- AddMetaData(object = hahn.combined, metadata = as.vector(Age.gp)   , col.name = "Age_Group"   )
# hahn.combined <- AddMetaData(object = hahn.combined, metadata = as.vector(Tissue)   , col.name = "Tissue"      )
# ################################################################################################
# 
# 
# ################################################################################################
# #### 2. Basic QC and filtering with Seurat and decontX
# 
# ### No filtering on genes at this stage
# hahn.combined <- SetIdent(hahn.combined, value = "SampleID")
# 
# # DecontX contamination levels for filtration
# pdf(paste(Sys.Date(),"Hahn_Brain_Cohort1_violinPlots_QC_DecontX.pdf", sep = "_"), height = 5, width = 10)
# VlnPlot(object = hahn.combined, features = c("decontX_contamination"), pt.size = 0)
# dev.off()
# 
# # The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.
# # The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
# hahn.combined[["percent.mito"]] <- PercentageFeatureSet(hahn.combined, pattern = "^MT-")
# head(hahn.combined@meta.data)
# 
# pdf(paste(Sys.Date(),"Hahn_Brain_Cohort1_violinPlots_QC_gene_UMI_mito.pdf", sep = "_"), height = 5, width = 10)
# VlnPlot(object = hahn.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
# dev.off()
# 
# # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# 
# plot1 <- FeatureScatter(hahn.combined, feature1 = "nCount_RNA", feature2 = "percent.mito")
# plot2 <- FeatureScatter(hahn.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# 
# pdf(paste(Sys.Date(),"Hahn_Brain_Cohort1_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
# plot1 + plot2
# dev.off()
# 
# # filter dead/low Q cells
# hahn.combined <- subset(hahn.combined, subset = nFeature_RNA > 250 & nFeature_RNA < 5000 & percent.mito < 10 & nCount_RNA < 25000 & decontX_contamination < 0.25 )
# hahn.combined
# # An object of class Seurat 
# # 32285 features across 97865 samples within 1 assay 
# # Active assay: RNA (32285 features, 0 variable features)
# 
# ### Check data after cell filtering
# head(hahn.combined@meta.data)
# 
# table(hahn.combined@meta.data$SampleID)
# # CP_O_1  CP_O_2  CP_Y_1  CP_Y_2 HIP_O_1 HIP_O_2 HIP_Y_1 HIP_Y_2 
# #   13752   17411   11015   13730    8428    6004   10282   17243 
# 
# #### Normalize the data for doublet analysis, etc
# # global-scaling normalization method 'LogNormalize' normalizes gene expression measurements for each cell 
# # by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
# hahn.combined <- NormalizeData(object = hahn.combined, normalization.method = "LogNormalize",  scale.factor = 10000)
# ################################################################################################
# 
# ################################################################################################
# #### 2b. Cell cycle prediction and storage
# # Read in a list of cell cycle markers, from Tirosh et al, 2015
# cc.genes <- readLines(con = "../../../../../cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")
# 
# # make into mouse gene names
# firstup <- function(x) {
#   substr(x, 1, 1) <- toupper(substr(x, 1, 1))
#   x
# }
# 
# cc.genes.mouse <- firstup(tolower(cc.genes))
# 
# # We can segregate this list into markers of G2/M phase and markers of S
# # phase
# s.genes   <- cc.genes.mouse[1:43]
# g2m.genes <- cc.genes.mouse[44:97]
# 
# # Assign Cell-Cycle Scores
# hahn.combined <- CellCycleScoring(object = hahn.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# 
# # write predictions to file
# write.table(hahn.combined@meta.data, file = paste0(Sys.Date(),"_Hahn_CellCycle_predictions.txt"), sep = "\t", quote = F)
# ################################################################################################
# 
# ################################################################################################
# #### 3. Find and remove doublets using doublet finder & scds workflow
# # save(hahn.combined, file = paste0(Sys.Date(),"_Hahn_Seurat_object_preSCT.RData"))
# 
# # https://github.com/chris-mcginnis-ucsf/DoubletFinder
# hahn.combined <- SCTransform(object = hahn.combined, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mito", "Phase"))
# save(hahn.combined, file = paste0(Sys.Date(),"_Hahn_Seurat_object_postSCT.RData"))
# 

load("2024-03-25_Hahn_Seurat_object_postSCT.RData")

# Run first pass analysis just for doublet identification (not final clustering)
hahn.combined <- RunPCA(hahn.combined, npcs = 30)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_Hahn_ElbowPlot.pdf"))
ElbowPlot(hahn.combined, ndims = 30)
dev.off()

# run dimensionality reduction
# Keep all PCs here, we'll do the clean clustering analysis on the merged object across all cohorts
hahn.combined <- RunUMAP(hahn.combined, dims = 1:30)
hahn.combined <- FindNeighbors(hahn.combined, dims = 1:30)
hahn.combined <- FindClusters(object = hahn.combined)

#### need to split by 10x sample to make sure to identify real doublets
# will run on one object at a time
cohort.1.list <- SplitObject(hahn.combined, split.by = "SampleID")

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
  nExp_poi       <- round((pred.dblt.rate[i]) *length(cohort.1.list[[i]]@meta.data$SampleID))      ## Assume doublets based on nuclei isolation protocol performance
  
  ## Run DoubletFinder with varying classification stringencies
  cohort.1.list[[i]] <- doubletFinder_v3(cohort.1.list[[i]], PCs = 1:30, pN = 0.25, pK = pk.killi, nExp = nExp_poi,     reuse.pANN = FALSE, sct = T)
  
  # get classification name
  my.DF.res.col <- colnames(cohort.1.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(cohort.1.list[[i]]@meta.data))]
  
  # rename column to enable subsetting
  colnames(cohort.1.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(cohort.1.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# run UMAP plots
for (i in 1:length(cohort.1.list)) {
  
  pdf(paste(Sys.Date(),"Hahn_Tissue",names(cohort.1.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(cohort.1.list[[i]], reduction = "umap", group.by = "DoubletFinder"), raster = T)
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
hahn.singlets.annot.c1 <- merge(cohort.1.list[[1]],
                                 y = c(cohort.1.list[[2]],
                                       cohort.1.list[[3]],
                                       cohort.1.list[[4]],
                                       cohort.1.list[[5]],
                                       cohort.1.list[[6]],
                                       cohort.1.list[[7]],
                                       cohort.1.list[[8]]),
                                 project = "Hahn")
hahn.singlets.annot.c1


# remove pANN columns that are 10xGenomics library lane specific
hahn.singlets.annot.c1@meta.data <- hahn.singlets.annot.c1@meta.data[,-grep("pANN",colnames(hahn.singlets.annot.c1@meta.data))]


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
  
  pdf(paste(Sys.Date(),"Hahn_Tissue",names(cohort.1.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to doubletFinder annotated Seurat object
hahn.singlets.annot.c1@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(cohort.1.list.scds)) {
  
  # for each object compare and move doublet annotations over
  hahn.singlets.annot.c1@meta.data[colnames(cohort.1.list.scds[[i]]), ]$scds_hybrid <- cohort.1.list.scds[[i]]$scds
  
}

## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& C. Merge and summarize doublet findings

table(hahn.singlets.annot.c1@meta.data$DoubletFinder, hahn.singlets.annot.c1@meta.data$scds_hybrid)
#          Doublet Singlet
#  Doublet    1425   10647
#  Singlet   10647   75146

# Union (more conservative)
hahn.singlets.annot.c1@meta.data$DoubletCall <- ifelse( bitOr(hahn.singlets.annot.c1@meta.data$DoubletFinder == "Doublet", hahn.singlets.annot.c1@meta.data$scds_hybrid == "Doublet") > 0, 
                                                         "Doublet", "Singlet")
table(hahn.singlets.annot.c1@meta.data$DoubletCall)
# Doublet Singlet 
#     22719   75146

# re-run dimensionality reduction for plotting purposes
hahn.singlets.annot.c1 <- SCTransform(object = hahn.singlets.annot.c1, vars.to.regress =  c("nFeature_RNA", "nCount_RNA", "percent.mito"))
hahn.singlets.annot.c1 <- RunPCA(hahn.singlets.annot.c1, npcs = 30)
hahn.singlets.annot.c1 <- RunUMAP(hahn.singlets.annot.c1, dims = 1:30)

pdf(paste0(Sys.Date(),"_Hahn_UMAP_Singlets_labelled_UNION.pdf"), width = 6, height = 5)
DimPlot(hahn.singlets.annot.c1, reduction = "umap", group.by = "DoubletCall", raster = T)
dev.off()

# save annotated object
save(hahn.singlets.annot.c1, file = paste0(Sys.Date(),"_Hahn_Brain_Seurat_object_with_AnnotatedDoublets.RData"))


### extract/subset only singlets
# save data for singlets df
hahn.singlets   <- subset(hahn.singlets.annot.c1, subset = DoubletCall %in% "Singlet")  # only keep singlets
hahn.singlets
# An object of class Seurat 
# 52640 features across 37963 samples within 2 assays 
# Active assay: SCT (24806 features, 3000 variable features)
#  1 other assay present: RNA
#  2 dimensional reductions calculated: pca, umap

pdf(paste0(Sys.Date(),"_Hahn_UMAP_Singlets_ONLY_UNION.pdf"), width = 6, height = 5)
DimPlot(hahn.singlets, reduction = "umap", group.by = "DoubletCall", raster = T)
dev.off()


# save filtered/annotated object
save(hahn.singlets, file = paste0(Sys.Date(),"_Hahn_Brain_Seurat_object_SINGLETS_ONLY.RData"))
################################################################################################


#######################
sink(file = paste(Sys.Date(),"_Hahn_Brain_Data_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()




