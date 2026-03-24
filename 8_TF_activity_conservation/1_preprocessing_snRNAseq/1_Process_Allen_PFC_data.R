setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Mouse_Datasets_for_Comparison/Processed_Objects/Allen_GSE207848')
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

theme_set(theme_bw())   

# 2024-03-22
# Parse/annotate Allen Brain data


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



# ################################################################################################################################################################
# #### 1. read data and metadata, create new Seurat and SCE objects
# 
# # read in deposited cellranger object
# Y1      <- Read10X_h5("GSM6321073_PFC_4wk_1_matrix.h5")
# Y2      <- Read10X_h5("GSM6321074_PFC_4wk_2_matrix.h5")
# Y3      <- Read10X_h5("GSM6321075_PFC_4wk_3_matrix.h5")
# Y4      <- Read10X_h5("GSM6321076_PFC_4wk_4_matrix.h5")
# O1      <- Read10X_h5("GSM6321077_PFC_90wk_1_matrix.h5")
# O2      <- Read10X_h5("GSM6321078_PFC_90wk_2_matrix.h5")
# O3      <- Read10X_h5("GSM6321079_PFC_90wk_3_matrix.h5")
# O4      <- Read10X_h5("GSM6321080_PFC_90wk_4_matrix.h5")
# 
# # Create SingleCellExperiment objects
# sce.Y1     <- SingleCellExperiment(list(counts = Y1     ))
# sce.Y2     <- SingleCellExperiment(list(counts = Y2     ))
# sce.Y3     <- SingleCellExperiment(list(counts = Y3     ))
# sce.Y4     <- SingleCellExperiment(list(counts = Y4     ))
# sce.O1     <- SingleCellExperiment(list(counts = O1     ))
# sce.O2     <- SingleCellExperiment(list(counts = O2     ))
# sce.O3     <- SingleCellExperiment(list(counts = O3     ))
# sce.O4     <- SingleCellExperiment(list(counts = O4     ))
# 
# # Run decontX (no available background)
# sce.Y1   <- decontX(sce.Y1  )
# sce.Y2   <- decontX(sce.Y2  )
# sce.Y3   <- decontX(sce.Y3  )
# sce.Y4   <- decontX(sce.Y4  )
# sce.O1   <- decontX(sce.O1  )
# sce.O2   <- decontX(sce.O2  )
# sce.O3   <- decontX(sce.O3  )
# sce.O4   <- decontX(sce.O4  )
# 
# # get seurat objects
# seurat.Y1 <- CreateSeuratObject( round(decontXcounts(sce.Y1 )) )
# seurat.Y2 <- CreateSeuratObject( round(decontXcounts(sce.Y2 )) )
# seurat.Y3 <- CreateSeuratObject( round(decontXcounts(sce.Y3 )) )
# seurat.Y4 <- CreateSeuratObject( round(decontXcounts(sce.Y4 )) )
# seurat.O1 <- CreateSeuratObject( round(decontXcounts(sce.O1 )) )
# seurat.O2 <- CreateSeuratObject( round(decontXcounts(sce.O2 )) )
# seurat.O3 <- CreateSeuratObject( round(decontXcounts(sce.O3 )) )
# seurat.O4 <- CreateSeuratObject( round(decontXcounts(sce.O4 )) )
# 
# seurat.Y1   <- AddMetaData(object = seurat.Y1 , colData(sce.Y1 )$decontX_contamination  , col.name = "decontX_contamination"  )
# seurat.Y2   <- AddMetaData(object = seurat.Y2 , colData(sce.Y2 )$decontX_contamination  , col.name = "decontX_contamination"  )
# seurat.Y3   <- AddMetaData(object = seurat.Y3 , colData(sce.Y3 )$decontX_contamination  , col.name = "decontX_contamination"  )
# seurat.Y4   <- AddMetaData(object = seurat.Y4 , colData(sce.Y4 )$decontX_contamination  , col.name = "decontX_contamination"  )
# seurat.O1   <- AddMetaData(object = seurat.O1 , colData(sce.O1 )$decontX_contamination  , col.name = "decontX_contamination"  )
# seurat.O2   <- AddMetaData(object = seurat.O2 , colData(sce.O2 )$decontX_contamination  , col.name = "decontX_contamination"  )
# seurat.O3   <- AddMetaData(object = seurat.O3 , colData(sce.O3 )$decontX_contamination  , col.name = "decontX_contamination"  )
# seurat.O4   <- AddMetaData(object = seurat.O4 , colData(sce.O4 )$decontX_contamination  , col.name = "decontX_contamination"  )
# 
# # Merge objects for the cohort
# allen.combined <- merge(seurat.Y1,
#                        y =  c(seurat.Y2 ,
#                               seurat.Y3 ,
#                               seurat.Y4 ,
#                               seurat.O1 ,
#                               seurat.O2 ,
#                               seurat.O3 ,
#                               seurat.O4 ),
#                        add.cell.ids = c("Y1"  ,
#                                         "Y2"  ,
#                                         "Y3"  ,
#                                         "Y4"  ,
#                                         "O1"   ,
#                                         "O2"   ,
#                                         "O3"   ,
#                                         "O4"   ),
#                        project = "Allen_Brain_Aging")
# allen.combined
# # An object of class Seurat 
# # 32285 features across 102831 samples within 1 assay 
# # Active assay: RNA (32285 features, 0 variable features)
# 
# #clean memory
# rm( Y1,     Y2,     Y3,     Y4,     O1,     O2,     O3,     O4,     sce.Y1,     sce.Y2,     sce.Y3,     sce.Y4,     sce.O1,     sce.O2,     sce.O3,     sce.O4,     seurat.Y1,     seurat.Y2,     seurat.Y3,     seurat.Y4,     seurat.O1,     seurat.O2,     seurat.O3,     seurat.O4     )
# ################################################################################################
# 
# 
# ################################################################################################
# #### 1b. Add key metadata to Seurat object
# 
# # create Group label
# my.Y1   <- grep("Y1"  , colnames(allen.combined@assays$RNA))
# my.Y2   <- grep("Y2"  , colnames(allen.combined@assays$RNA))
# my.Y3   <- grep("Y3"  , colnames(allen.combined@assays$RNA))
# my.Y4   <- grep("Y4"  , colnames(allen.combined@assays$RNA))
# my.O1   <- grep("O1"  , colnames(allen.combined@assays$RNA))
# my.O2   <- grep("O2"  , colnames(allen.combined@assays$RNA))
# my.O3   <- grep("O3"  , colnames(allen.combined@assays$RNA))
# my.O4   <- grep("O4"  , colnames(allen.combined@assays$RNA))
# 
# #####
# SampleID <- rep("NA", length(colnames(allen.combined@assays$RNA)))
# SampleID[ my.Y1  ]   <- "Y1"
# SampleID[ my.Y2  ]   <- "Y2"
# SampleID[ my.Y3  ]   <- "Y3"
# SampleID[ my.Y4  ]   <- "Y4"
# SampleID[ my.O1  ]   <- "O1"
# SampleID[ my.O2  ]   <- "O2"
# SampleID[ my.O3  ]   <- "O3"
# SampleID[ my.O4  ]   <- "O4"
# SampleID <- data.frame(SampleID)
# rownames(SampleID) <- colnames(allen.combined@assays$RNA)
# 
# #####
# Sex <- rep("NA", length(colnames(allen.combined@assays$RNA)))
# Sex[ my.Y1 ]   <- "F"  
# Sex[ my.Y2 ]   <- "F" 
# Sex[ my.Y3 ]   <- "F" 
# Sex[ my.Y4 ]   <- "F"
# Sex[ my.O1 ]   <- "F" 
# Sex[ my.O2 ]   <- "F"  
# Sex[ my.O3 ]   <- "F"   
# Sex[ my.O4 ]   <- "F" 
# Sex <- data.frame(Sex)
# rownames(Sex) <- colnames(allen.combined@assays$RNA)
# 
# ##### Y, M, O, G
# Age.gp <- rep("NA", length(colnames(allen.combined@assays$RNA)))
# Age.gp[ my.Y1 ]   <- "Y"
# Age.gp[ my.Y2 ]   <- "Y"
# Age.gp[ my.Y3 ]   <- "Y"
# Age.gp[ my.Y4 ]   <- "Y"
# Age.gp[ my.O1 ]   <- "O"
# Age.gp[ my.O2 ]   <- "O"
# Age.gp[ my.O3 ]   <- "O"
# Age.gp[ my.O4 ]   <- "O"
# Age.gp <- data.frame(Age.gp)
# rownames(Age.gp) <- colnames(allen.combined@assays$RNA)
# 
# ##### (average of the pool) 
# Age.m <- rep("NA", length(colnames(allen.combined@assays$RNA)))
# Age.m[ my.Y1 ]   <- 1   #
# Age.m[ my.Y2 ]   <- 1   #
# Age.m[ my.Y3 ]   <- 1   #
# Age.m[ my.Y4 ]   <- 1   #
# Age.m[ my.O1 ]   <- 21   #
# Age.m[ my.O2 ]   <- 21   #
# Age.m[ my.O3 ]   <- 21   #
# Age.m[ my.O4 ]   <- 21   #
# Age.m <- data.frame(Age.m)
# rownames(Age.m) <- colnames(allen.combined@assays$RNA)
# 
# 
# # update Seurat with metadata
# allen.combined <- AddMetaData(object = allen.combined, metadata = as.vector(SampleID) , col.name = "SampleID"       )
# allen.combined <- AddMetaData(object = allen.combined, metadata = as.vector(Sex)      , col.name = "Sex"         )
# allen.combined <- AddMetaData(object = allen.combined, metadata = as.vector(Age.m)    , col.name = "Age_months"   )
# allen.combined <- AddMetaData(object = allen.combined, metadata = as.vector(Age.gp)   , col.name = "Age_Group"   )
# ################################################################################################
# 
# 
# ################################################################################################
# #### 2. Basic QC and filtering with Seurat and decontX
# 
# ### No filtering on genes at this stage - only after all cohorts merged for fairness
# allen.combined <- SetIdent(allen.combined, value = "SampleID")
# 
# # DecontX contamination levels for filtration
# pdf(paste(Sys.Date(),"Allen_violinPlots_QC_DecontX.pdf", sep = "_"), height = 5, width = 10)
# VlnPlot(object = allen.combined, features = c("decontX_contamination"), pt.size = 0)
# dev.off()
# 
# # The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.
# # The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
# allen.combined[["percent.mito"]] <- PercentageFeatureSet(allen.combined, pattern = "^MT-")
# head(allen.combined@meta.data)
# 
# pdf(paste(Sys.Date(),"Allen_violinPlots_QC_gene_UMI_mito.pdf", sep = "_"), height = 5, width = 10)
# VlnPlot(object = allen.combined, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
# dev.off()
# 
# # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# 
# plot1 <- FeatureScatter(allen.combined, feature1 = "nCount_RNA", feature2 = "percent.mito")
# plot2 <- FeatureScatter(allen.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# 
# pdf(paste(Sys.Date(),"Allen_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
# plot1 + plot2
# dev.off()
# 
# # filter dead/low Q cells
# allen.combined <- subset(allen.combined, subset = nFeature_RNA > 250 & nFeature_RNA < 5000 & percent.mito < 10 & nCount_RNA < 25000 & decontX_contamination < 0.25 )
# allen.combined
# # An object of class Seurat 
# # 32285 features across 69152 samples within 1 assay 
# # Active assay: RNA (32285 features, 0 variable features)
# 
# ### Check data after cell filtering
# head(allen.combined@meta.data)
# 
# table(allen.combined@meta.data$SampleID)
# #    O1    O2    O3    O4    Y1    Y2    Y3    Y4 
# #  9268  7952  7303  8353 10927  9364  7594  8391 
# 
# #### Normalize the data for doublet analysis, etc
# # global-scaling normalization method 'LogNormalize' normalizes gene expression measurements for each cell 
# # by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
# allen.combined <- NormalizeData(object = allen.combined, normalization.method = "LogNormalize",  scale.factor = 10000)
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
# allen.combined <- CellCycleScoring(object = allen.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# 
# # write predictions to file
# write.table(allen.combined@meta.data, file = paste0(Sys.Date(),"_Allen_CellCycle_predictions.txt"), sep = "\t", quote = F)
# ################################################################################################
# 
# 
# ################################################################################################
# #### 3. Find and remove doublets using doublet finder & scds workflow
# 
# # https://github.com/chris-mcginnis-ucsf/DoubletFinder
# allen.combined <- SCTransform(object = allen.combined, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mito", "Phase"))
# save(allen.combined, file = paste0(Sys.Date(),"_Allen_Seurat_object_postSCT.RData"))

load('2024-03-25_Allen_Seurat_object_postSCT.RData')

# Run first pass analysis just for doublet identification (not final clustering)
allen.combined <- RunPCA(allen.combined, npcs = 30)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_Allen_ElbowPlot.pdf"))
ElbowPlot(allen.combined, ndims = 30)
dev.off()

# run dimensionality reduction
# Keep all PCs here, we'll do the clean clustering analysis on the merged object across all cohorts
allen.combined <- RunUMAP(allen.combined, dims = 1:30)
allen.combined <- FindNeighbors(allen.combined, dims = 1:30)
allen.combined <- FindClusters(object = allen.combined)

#### need to split by 10x sample to make sure to identify real doublets
# will run on one object at a time
cohort.1.list <- SplitObject(allen.combined, split.by = "SampleID")

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
  
  pdf(paste(Sys.Date(),"Allen_Tissue",names(cohort.1.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(cohort.1.list[[i]], reduction = "umap", group.by = "DoubletFinder"), raster = T)
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
allen.singlets.annot.c1 <- merge(cohort.1.list[[1]],
                                 y = c(cohort.1.list[[2]],
                                       cohort.1.list[[3]],
                                       cohort.1.list[[4]],
                                       cohort.1.list[[5]],
                                       cohort.1.list[[6]],
                                       cohort.1.list[[7]],
                                       cohort.1.list[[8]]),
                                 project = "Allen")
allen.singlets.annot.c1
# An object of class Seurat 
# 58305 features across 69152 samples within 2 assays 
# Active assay: SCT (26020 features, 0 variable features)
# 1 other assay present: RNA

# remove pANN columns that are 10xGenomics library lane specific
allen.singlets.annot.c1@meta.data <- allen.singlets.annot.c1@meta.data[,-grep("pANN",colnames(allen.singlets.annot.c1@meta.data))]


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
  
  pdf(paste(Sys.Date(),"Allen_Tissue",names(cohort.1.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to doubletFinder annotated Seurat object
allen.singlets.annot.c1@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(cohort.1.list.scds)) {
  
  # for each object compare and move doublet annotations over
  allen.singlets.annot.c1@meta.data[colnames(cohort.1.list.scds[[i]]), ]$scds_hybrid <- cohort.1.list.scds[[i]]$scds
  
}

## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& C. Merge and summarize doublet findings

table(allen.singlets.annot.c1@meta.data$DoubletFinder, allen.singlets.annot.c1@meta.data$scds_hybrid)
#         Doublet Singlet
#   Doublet     500    5089
#   Singlet    5089   58474

# Union (more conservative)
allen.singlets.annot.c1@meta.data$DoubletCall <- ifelse( bitOr(allen.singlets.annot.c1@meta.data$DoubletFinder == "Doublet", allen.singlets.annot.c1@meta.data$scds_hybrid == "Doublet") > 0, 
                                                         "Doublet", "Singlet")
table(allen.singlets.annot.c1@meta.data$DoubletCall)
# Doublet Singlet 
#     10678   58474

# re-run dimensionality reduction for plotting purposes
allen.singlets.annot.c1 <- SCTransform(object = allen.singlets.annot.c1, vars.to.regress =  c("nFeature_RNA", "nCount_RNA", "percent.mito"))
allen.singlets.annot.c1 <- RunPCA(allen.singlets.annot.c1, npcs = 30)
allen.singlets.annot.c1 <- RunUMAP(allen.singlets.annot.c1, dims = 1:30)

pdf(paste0(Sys.Date(),"_Allen_UMAP_Singlets_labelled_UNION.pdf"), width = 6, height = 5)
DimPlot(allen.singlets.annot.c1, reduction = "umap", group.by = "DoubletCall", raster = T)
dev.off()

# save annotated object
save(allen.singlets.annot.c1, file = paste0(Sys.Date(),"_Allen_Seurat_object_with_AnnotatedDoublets.RData"))


### extract/subset only singlets
# save data for singlets df
allen.singlets   <- subset(allen.singlets.annot.c1, subset = DoubletCall %in% "Singlet")  # only keep singlets
allen.singlets
# An object of class Seurat 
# 58305 features across 58474 samples within 2 assays 
# Active assay: SCT (26020 features, 3000 variable features)
# 1 other assay present: RNA

pdf(paste0(Sys.Date(),"_Allen_UMAP_Singlets_ONLY_UNION.pdf"), width = 6, height = 5)
DimPlot(allen.singlets, reduction = "umap", group.by = "DoubletCall", raster = T)
dev.off()

table(allen.singlets@meta.data$SampleID)
# O1   O2   O3   O4   Y1   Y2   Y3   Y4 
# 7690 6795 6377 7116 8832 7776 6630 7258


# save filtered/annotated object
save(allen.singlets, file = paste0(Sys.Date(),"_Allen_Seurat_object_SINGLETS_ONLY.RData"))
################################################################################################


#######################
sink(file = paste(Sys.Date(),"_Allen_Brain_Data_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()


