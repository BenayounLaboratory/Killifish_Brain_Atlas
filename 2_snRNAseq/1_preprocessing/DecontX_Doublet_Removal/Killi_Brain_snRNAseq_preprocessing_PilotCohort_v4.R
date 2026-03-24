setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Preprocessing/DecontX/Pilot_Cohort/')
options(stringsAsFactors = F)

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
cts.GRZ_Y_F_p        <- Read10X("/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-01-20_Killifish_tissue_scRNAseq_cohort3/CellRanger/GRZ_Brain_F_FishTEDB_NR/outs/filtered_feature_bc_matrix/")
cts.GRZ_Y_M_p        <- Read10X("/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-01-20_Killifish_tissue_scRNAseq_cohort3/CellRanger/GRZ_Brain_M_FishTEDB_NR/outs/filtered_feature_bc_matrix/")
cts_raw.GRZ_Y_F_p    <- Read10X("/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-01-20_Killifish_tissue_scRNAseq_cohort3/CellRanger/GRZ_Brain_F_FishTEDB_NR/outs/raw_feature_bc_matrix/")
cts_raw.GRZ_Y_M_p    <- Read10X("/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/2022-01-20_Killifish_tissue_scRNAseq_cohort3/CellRanger/GRZ_Brain_M_FishTEDB_NR/outs/raw_feature_bc_matrix/")

# Create SingleCellExperiment objects
sce.GRZ_Y_F_p        <- SingleCellExperiment(list(counts = cts.GRZ_Y_F_p     ))
sce.GRZ_Y_M_p        <- SingleCellExperiment(list(counts = cts.GRZ_Y_M_p     ))
sce_raw.GRZ_Y_F_p    <- SingleCellExperiment(list(counts = cts_raw.GRZ_Y_F_p ))
sce_raw.GRZ_Y_M_p    <- SingleCellExperiment(list(counts = cts_raw.GRZ_Y_M_p ))

# Run decontX
sce.GRZ_Y_F_p        <- decontX(sce.GRZ_Y_F_p, background = sce_raw.GRZ_Y_F_p )
sce.GRZ_Y_M_p        <- decontX(sce.GRZ_Y_M_p, background = sce_raw.GRZ_Y_M_p )


# get seurat objects
seurat.GRZ_Y_F_p  <- CreateSeuratObject(round(decontXcounts(sce.GRZ_Y_F_p )) )
seurat.GRZ_Y_M_p  <- CreateSeuratObject(round(decontXcounts(sce.GRZ_Y_M_p )) )

seurat.GRZ_Y_F_p  <- AddMetaData(object = seurat.GRZ_Y_F_p, colData(sce.GRZ_Y_F_p)$decontX_contamination  , col.name = "decontX_contamination"  )
seurat.GRZ_Y_M_p  <- AddMetaData(object = seurat.GRZ_Y_M_p, colData(sce.GRZ_Y_M_p)$decontX_contamination  , col.name = "decontX_contamination"  )

# Merge objects for the cohort
brain.cp.combined <- merge(seurat.GRZ_Y_F_p,
                           y =  c(seurat.GRZ_Y_M_p),
                           add.cell.ids = c("GRZ_Y_F_p"  ,
                                            "GRZ_Y_M_p"  ),
                           project = "10x_Killi_Brain_Aging")
brain.cp.combined
# An object of class Seurat 
# 27834 features across 14586 samples within 1 assay 
# Active assay: RNA (27834 features, 0 variable features)

rm(cts.GRZ_Y_F_p,cts.GRZ_Y_M_p,cts_raw.GRZ_Y_F_p,cts_raw.GRZ_Y_M_p,sce.GRZ_Y_F_p,sce.GRZ_Y_M_p,sce_raw.GRZ_Y_F_p,sce_raw.GRZ_Y_M_p,seurat.GRZ_Y_F_p,seurat.GRZ_Y_M_p)
################################################################################################


################################################################################################
#### 1b. Add key metadata to Seurat object

# create Group label
my.GRZ_Y_F_p   <- grep("GRZ_Y_F_p"  , colnames(brain.cp.combined@assays$RNA))
my.GRZ_Y_M_p   <- grep("GRZ_Y_M_p"  , colnames(brain.cp.combined@assays$RNA))

#####
Group <- rep("NA", length(colnames(brain.cp.combined@assays$RNA)))
Group[ my.GRZ_Y_F_p  ]   <- "GRZ_Y_F"
Group[ my.GRZ_Y_M_p  ]   <- "GRZ_Y_M"
Group <- data.frame(Group)
rownames(Group) <- colnames(brain.cp.combined@assays$RNA)

#####
Sex <- rep("NA", length(colnames(brain.cp.combined@assays$RNA)))
Sex[ my.GRZ_Y_F_p ]   <- "F" 
Sex[ my.GRZ_Y_M_p ]   <- "M"
Sex <- data.frame(Sex)
rownames(Sex) <- colnames(brain.cp.combined@assays$RNA)

##### Y, M, O, G
Age.gp <- rep("NA", length(colnames(brain.cp.combined@assays$RNA)))
Age.gp[ my.GRZ_Y_F_p  ]   <- "Y"
Age.gp[ my.GRZ_Y_M_p  ]   <- "Y"
Age.gp <- data.frame(Age.gp)
rownames(Age.gp) <- colnames(brain.cp.combined@assays$RNA)

##### (average of the pool) (6.2857 + 6.381)/2 = 6.33335
Age.w <- rep("NA", length(colnames(brain.cp.combined@assays$RNA)))
Age.w[ my.GRZ_Y_F_p  ]   <- 6.3
Age.w[ my.GRZ_Y_M_p  ]   <- 6.3
Age.w <- data.frame(Age.w)
rownames(Age.w) <- colnames(brain.cp.combined@assays$RNA)

#####
Strain <- rep("NA", length(colnames(brain.cp.combined@assays$RNA)))
Strain[ my.GRZ_Y_F_p ]   <- "GRZ"
Strain[ my.GRZ_Y_M_p ]   <- "GRZ"
Strain <- data.frame(Strain)
rownames(Strain) <- colnames(brain.cp.combined@assays$RNA)


#####
Batch <- rep("Pilot", length(colnames(brain.cp.combined@assays$RNA)))
Batch <- data.frame(Batch)
rownames(Batch) <- colnames(brain.cp.combined@assays$RNA)

# update Seurat with metadata
brain.cp.combined <- AddMetaData(object = brain.cp.combined, metadata = as.vector(Group)    , col.name = "Group"       )
brain.cp.combined <- AddMetaData(object = brain.cp.combined, metadata = as.vector(Sex)      , col.name = "Sex"         )
brain.cp.combined <- AddMetaData(object = brain.cp.combined, metadata = as.vector(Age.w)    , col.name = "Age_weeks"   )
brain.cp.combined <- AddMetaData(object = brain.cp.combined, metadata = as.vector(Age.gp)   , col.name = "Age_Group"   )
brain.cp.combined <- AddMetaData(object = brain.cp.combined, metadata = as.vector(Strain)   , col.name = "Strain"      )
brain.cp.combined <- AddMetaData(object = brain.cp.combined, metadata = as.vector(Batch)    , col.name = "Batch"       )
################################################################################################


################################################################################################
#### 2. Basic QC and filtering with Seurat

### No filtering on genes at this stage - only after all cohorts merged for fairness
brain.cp.combined <- SetIdent(brain.cp.combined, value = "Group")

# DecontX contamination levels for filtration
pdf(paste(Sys.Date(),"Killifish_Brain_Pilot_Cohort_violinPlots_QC_DecontX.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = brain.cp.combined, features = c("decontX_contamination"), pt.size = 0)
dev.off()

# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
brain.cp.combined[["percent.mito"]] <- PercentageFeatureSet(brain.cp.combined, pattern = "^MT-")
head(brain.cp.combined@meta.data)

pdf(paste(Sys.Date(),"Killifish_Brain_Pilot_Cohort_violinPlots_QC_gene_UMI_mito.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = brain.cp.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(brain.cp.combined, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(brain.cp.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste(Sys.Date(),"Killifish_Brain_Pilot_Cohort_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
plot1 + plot2
dev.off()

# filter dead/low Q cell, or eggregious doublets
brain.cp.combined <- subset(brain.cp.combined, subset = nFeature_RNA > 250 & nFeature_RNA < 5000 & percent.mito < 10 & nCount_RNA < 25000 & decontX_contamination < 0.25 )
brain.cp.combined
# An object of class Seurat 
# 27834 features across 13307 samples within 1 assay 
# Active assay: RNA (27834 features, 0 variable features)

### Check data after cell filtering
head(brain.cp.combined@meta.data)
#                                 orig.ident nCount_RNA nFeature_RNA decontX_contamination   Group Sex Age_weeks Age_Group Strain Batch percent.mito
# GRZ_Y_F_p_AAACCCACAGCGTAGA-1 SeuratProject       2537         1555           0.071846126 GRZ_Y_F   F       6.3         Y    GRZ Pilot   0.07883327
# GRZ_Y_F_p_AAACCCACAGGTATGG-1 SeuratProject       1616         1205           0.101865036 GRZ_Y_F   F       6.3         Y    GRZ Pilot   0.00000000
# GRZ_Y_F_p_AAACCCACATGTCTAG-1 SeuratProject        990          820           0.100146993 GRZ_Y_F   F       6.3         Y    GRZ Pilot   0.00000000
# GRZ_Y_F_p_AAACCCATCATTGTTC-1 SeuratProject       1613         1074           0.001576056 GRZ_Y_F   F       6.3         Y    GRZ Pilot   0.06199628
# GRZ_Y_F_p_AAACGAAAGACATGCG-1 SeuratProject       1359          990           0.064152904 GRZ_Y_F   F       6.3         Y    GRZ Pilot   0.14716703
# GRZ_Y_F_p_AAACGAAAGGTAAGTT-1 SeuratProject        829          643           0.052335460 GRZ_Y_F   F       6.3         Y    GRZ Pilot   0.00000000


table(brain.cp.combined@meta.data$Group)
# GRZ_Y_F GRZ_Y_M 
#      6626    6681

#### Normalize the data for doublet analysis, etc
# global-scaling normalization method 'LogNormalize' normalizes gene expression measurements for each cell 
# by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
brain.cp.combined <- NormalizeData(object = brain.cp.combined, normalization.method = "LogNormalize",  scale.factor = 10000)
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
brain.cp.combined <- CellCycleScoring(object = brain.cp.combined, s.features = s.genes.k, g2m.features = g2m.genes.k, set.ident = TRUE)

# write predictions to file
write.table(brain.cp.combined@meta.data, file = paste0(Sys.Date(),"_Killi_Brain_Pilot_Cohort_CellCycle_predictions.txt"), sep = "\t", quote = F)
################################################################################################


################################################################################################
#### 3. Find and remove doublets using doublet finder & scds workflow

# https://github.com/chris-mcginnis-ucsf/DoubletFinder
brain.cp.combined <- SCTransform(object = brain.cp.combined, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mito", "Phase"))
save(brain.cp.combined, file = paste0(Sys.Date(),"_Killi_Brain_Pilot_Cohort_Seurat_object_postSCT.RData"))

# Run first pass analysis just for doublet identification (not final clustering)
brain.cp.combined <- RunPCA(brain.cp.combined, npcs = 30)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_Killi_Brain_Pilot_Cohort_ElbowPlot.pdf"))
ElbowPlot(brain.cp.combined, ndims = 30)
dev.off()

# run dimensionality reduction
# Keep all PCs here, we'll do the clean clustering analysis on the merged object across all cohorts
brain.cp.combined <- RunUMAP(brain.cp.combined, dims = 1:30)
brain.cp.combined <- FindNeighbors(brain.cp.combined, dims = 1:30)
brain.cp.combined <- FindClusters(object = brain.cp.combined)

#### need to split by 10x sample to make sure to identify real doublets
# will run on one object at a time
cohortp.list <- SplitObject(brain.cp.combined, split.by = "Group")

## Assume doublet rate based on 10x information (add a 15% fudge factor due to nuclei being more sticky)
pred.dblt.rate <- 1.15 * predict(pred_dblt_lm, data.frame("cell_number" = unlist(lapply(cohortp.list, ncol))))/100

## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& A. Run DoubletFinder
# loop over samples
for (i in 1:length(cohortp.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list_killi <- paramSweep_v3(cohortp.list[[i]], PCs = 1:30, sct = TRUE, num.cores	 = 4)
  sweep.stats_killi    <- summarizeSweep(sweep.res.list_killi, GT = FALSE)
  bcmvn_killi          <- find.pK(sweep.stats_killi)
  
  # need some R gymnastics since the Pk is stored as a factor for some reason
  # to get the pK number, need to first convert to character and THEN to numeric
  # numeric first yield row number
  pk.killi <- as.numeric(as.character(bcmvn_killi[as.numeric(bcmvn_killi$pK[bcmvn_killi$BCmetric == max(bcmvn_killi$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(cohortp.list[[i]]@meta.data$seurat_clusters)             ## ex: annotations
  nExp_poi       <- round((pred.dblt.rate[i]) *length(cohortp.list[[i]]@meta.data$Group))   ## Assume doublets based on nuclei isolation protocol performance
  
  ## Run DoubletFinder with varying classification stringencies
  cohortp.list[[i]] <- doubletFinder_v3(cohortp.list[[i]], PCs = 1:30, pN = 0.25, pK = pk.killi, nExp = nExp_poi,     reuse.pANN = FALSE, sct = T)
  
  # get classification name
  my.DF.res.col <- colnames(cohortp.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(cohortp.list[[i]]@meta.data))]
  
  # rename column to enable subsetting
  colnames(cohortp.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(cohortp.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# run UMAP plots
for (i in 1:length(cohortp.list)) {
  
  pdf(paste(Sys.Date(),"Killifish_Tissue",names(cohortp.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(cohortp.list[[i]], reduction = "umap", group.by = "DoubletFinder"))
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
killi.singlets.annot.cp <- merge(cohortp.list[[1]],
                                 y = c(cohortp.list[[2]]),
                                 project = "Killi_Brain_Pilot_Cohort")
killi.singlets.annot.cp
# An object of class Seurat 
# 51150 features across 14572 samples within 2 assays 
# Active assay: SCT (23316 features, 0 variable features)
# 1 other assay present: RNA

# remove pANN columns that are 10xGenomics library lane specific
killi.singlets.annot.cp@meta.data <- killi.singlets.annot.cp@meta.data[,-grep("pANN",colnames(killi.singlets.annot.cp@meta.data))]


## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& B. Run scds:single cell doublet scoring (hybrid method)
# cxds is based on co-expression of gene pairs and works with absence/presence calls only, 
# bcds uses the full count information and a binary classification approach using artificially generated doublets. 
# cxds_bcds_hybrid combines both approaches

# create scds working object - convert list to SingleCellExperiment
cohortp.list.scds        <- lapply(cohortp.list, as.SingleCellExperiment)

# loop over sample
for (i in 1:length(cohortp.list.scds)) {
  
  # Annotate doublets using co-expression based doublet scoring:
  cohortp.list.scds[[i]] <- cxds_bcds_hybrid(cohortp.list.scds[[i]])
  
  # predicted doublet rate
  n.db <- round((pred.dblt.rate[i])*ncol(cohortp.list.scds[[i]]))                         ## Assume doublets based on nuclei isolation protocol performance
  
  # sort prediction, get top n.db cells
  srt.db.score <- sort(cohortp.list.scds[[i]]$hybrid_score, index.return = T, decreasing = T)
  cohortp.list.scds[[i]]$scds <- "Singlet"
  cohortp.list.scds[[i]]$scds[srt.db.score$ix[1:n.db]] <- "Doublet"
  
}

# run UMAP plots
for (i in 1:length(cohortp.list.scds)) {
  
  p <- plotReducedDim(cohortp.list.scds[[i]], dimred = "UMAP", colour_by = "scds")
  
  pdf(paste(Sys.Date(),"Killifish_Tissue",names(cohortp.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to doubletFinder annotated Seurat object
killi.singlets.annot.cp@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(cohortp.list.scds)) {
  
  # for each object compare and move doublet annotations over
  killi.singlets.annot.cp@meta.data[colnames(cohortp.list.scds[[i]]), ]$scds_hybrid <- cohortp.list.scds[[i]]$scds
  
}

## &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& C. Merge and summarize doublet findings

table(killi.singlets.annot.cp@meta.data$DoubletFinder, killi.singlets.annot.cp@meta.data$scds_hybrid)
#            Doublet Singlet
#    Doublet     155     660
#    Singlet     660   11832

# Union (more conservative)
killi.singlets.annot.cp@meta.data$DoubletCall <- ifelse( bitOr(killi.singlets.annot.cp@meta.data$DoubletFinder == "Doublet", killi.singlets.annot.cp@meta.data$scds_hybrid == "Doublet") > 0, 
                                                           "Doublet", "Singlet")
table(killi.singlets.annot.cp@meta.data$DoubletCall)
# Doublet Singlet 
#    1475   11832

# re-run dimensionality reduction for plotting purposes
killi.singlets.annot.cp <- SCTransform(object = killi.singlets.annot.cp, vars.to.regress =  c("nFeature_RNA", "nCount_RNA", "percent.mito"))
killi.singlets.annot.cp <- RunPCA(killi.singlets.annot.cp, npcs = 30)
killi.singlets.annot.cp <- RunUMAP(killi.singlets.annot.cp, dims = 1:30)

pdf(paste0(Sys.Date(),"_Killi_Brain_Pilot_Cohort_UMAP_Singlets_labelled_UNION.pdf"), width = 6, height = 5)
DimPlot(killi.singlets.annot.cp, reduction = "umap", group.by = "DoubletCall", raster = T)
dev.off()

pdf(paste0(Sys.Date(),"_Killi_Brain_Pilot_Cohort_Key_Marker_Gene_expression_plots.pdf"), width = 3.5, height = 3)
FeaturePlot(killi.singlets.annot.cp, features = "olig2", raster = T)
FeaturePlot(killi.singlets.annot.cp, features = "olig1", raster = T)
FeaturePlot(killi.singlets.annot.cp, features = "mpz"  , raster = T)
FeaturePlot(killi.singlets.annot.cp, features = "csf1r", raster = T)
dev.off()

# DimPlot(killi.singlets.annot.cp, reduction = "umap", group.by = "DoubletFinder")
# DimPlot(killi.singlets.annot.cp, reduction = "umap", group.by = "scds_hybrid")

# save annotated object
save(killi.singlets.annot.cp, file = paste0(Sys.Date(),"_Killifish_Brain_Pilot_Cohort_Seurat_object_with_AnnotatedDoublets.RData"))


### extract/subset only singlets
# save data for singlets df
killi.singlets.cp   <- subset(killi.singlets.annot.cp, subset = DoubletCall %in% "Singlet")  # only keep singlets
killi.singlets.cp
# An object of class Seurat
# 50859 features across 11864 samples within 2 assays
# Active assay: SCT (23025 features, 3000 variable features)
#  1 other assay present: RNA
#  2 dimensional reductions calculated: pca, umap

pdf(paste0(Sys.Date(),"_Killi_Brain_Pilot_Cohort_UMAP_Singlets_ONLY_UNION.pdf"), width = 6, height = 5)
DimPlot(killi.singlets.cp, reduction = "umap", group.by = "DoubletCall", raster = T)
dev.off()

table(killi.singlets.cp@meta.data$Group)
# GRZ_Y_F GRZ_Y_M 
#    5905    5959

# save filtered/annotated object
save(killi.singlets.cp, file = paste0(Sys.Date(),"_Killifish_Brain_Pilot_Cohort_Seurat_object_SINGLETS_ONLY.RData"))
################################################################################################


#######################
sink(file = paste(Sys.Date(),"_Killifish_Brain_Pilot_Cohort_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()
