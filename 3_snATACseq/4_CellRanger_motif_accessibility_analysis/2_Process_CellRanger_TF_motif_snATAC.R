setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/snATAC_Brain_Aging_Meta/Downstream_Analyses/Process_CR_TF')
options(stringsAsFactors = F)

#### Packages
library(Seurat)        # single cell general package
library(Signac)        # scATAC processing
library(bitops)        # 
library(Matrix)
library(Polychrome)

# 2024-04-24
# Process scATAC to get lists of barcodes per cell type/per replicate
# to generate sample specific cell type specific bam files
# https://github.com/10XGenomics/subset-bam


#########################################################################################
# 0. Load annotated object
load("../../Signac/2023-08-23_Killi_brain_ATAC_2Cohorts_FiltNorm_Singlets_ANNOTATEDfromRNA.RData")

brain.atac.singlets
# An object of class Seurat 
# 148827 features across 27171 samples within 2 assays 
# Active assay: RNA (22135 features, 0 variable features)
# 1 other assay present: ATAC
# 2 dimensional reductions calculated: lsi, umap

#########################################################################################
# 1. Load TF data from Cell Ranger

#https://bioinformatics.stackexchange.com/questions/9213/can-i-change-the-name-of-file-features-tsv-to-genes-tsv
import_tf_data <- function (matrix_dir, project_info) {
  mat           <- readMM(file = paste0(matrix_dir, "matrix.mtx.gz"))
  feature.names <- read.delim(paste0(matrix_dir, "motifs.tsv")     , header = FALSE,stringsAsFactors = FALSE)
  barcode.names <- read.delim(paste0(matrix_dir, "barcodes.tsv.gz"), header = FALSE,stringsAsFactors = FALSE)
  colnames(mat) <- barcode.names$V1
  rownames(mat) <- feature.names$V1
  object        <- CreateSeuratObject(mat, project=project_info,  assay = "TF")
  return(object)
}

YF1.tf <- import_tf_data(matrix_dir = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1/CellRanger/Killi_Brain_YF1_FishTEDB_NR/outs/filtered_tf_bc_matrix/", project_info = "CellRanger_TFData")
OF1.tf <- import_tf_data(matrix_dir = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1/CellRanger/Killi_Brain_OF1_FishTEDB_NR/outs/filtered_tf_bc_matrix/", project_info = "CellRanger_TFData")
YM1.tf <- import_tf_data(matrix_dir = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1/CellRanger/Killi_Brain_YM1_FishTEDB_NR/outs/filtered_tf_bc_matrix/", project_info = "CellRanger_TFData")
OM1.tf <- import_tf_data(matrix_dir = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1/CellRanger/Killi_Brain_OM1_FishTEDB_NR/outs/filtered_tf_bc_matrix/", project_info = "CellRanger_TFData")
YF2.tf <- import_tf_data(matrix_dir = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2/CellRanger/Killi_Brain_YF2_FishTEDB_NR/outs/filtered_tf_bc_matrix/", project_info = "CellRanger_TFData")
OF2.tf <- import_tf_data(matrix_dir = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2/CellRanger/Killi_Brain_OF2_FishTEDB_NR/outs/filtered_tf_bc_matrix/", project_info = "CellRanger_TFData")
YM2.tf <- import_tf_data(matrix_dir = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2/CellRanger/Killi_Brain_YM2_FishTEDB_NR/outs/filtered_tf_bc_matrix/", project_info = "CellRanger_TFData")
OM2.tf <- import_tf_data(matrix_dir = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2/CellRanger/Killi_Brain_OM2_FishTEDB_NR/outs/filtered_tf_bc_matrix/", project_info = "CellRanger_TFData")


# merge all datasets, adding a cell ID to make sure cell names are unique
brain.atac.tf <- merge(x = YF1.tf,
                       y = list(OF1.tf, YM1.tf, OM1.tf,
                                YF2.tf,OF2.tf, YM2.tf, OM2.tf),
                       add.cell.ids = c("YF1", "OF1", "YM1", "OM1",
                                        "YF2", "OF2", "YM2", "OM2"))
brain.atac.tf
# An object of class Seurat 
# 579 features across 40921 samples within 1 assay 
# Active assay: TF (579 features, 0 variable features)
#########################################################################################

#########################################################################################
# 2. Clean up for merging
# Select cell barcodes detected by both ATAC and TF
joint.bcs <- intersect(colnames(brain.atac.tf), colnames(brain.atac.singlets)) # 26942

brain.atac.tf.sub <- brain.atac.tf[,joint.bcs]
brain.atac.tf.sub
# An object of class Seurat 
# 579 features across 26942 samples within 1 assay 
# Active assay: TF (579 features, 0 variable features)

# create a Seurat object copy
brain.atac.tf.cp <- brain.atac.tf.sub

# Normalize RNA data with log normalization
brain.atac.tf.cp <- NormalizeData(brain.atac.tf.cp, normalization.method = "LogNormalize",  scale.factor = 10000)

# Normalization and linear dimensional reduction
# The combined steps of TF-IDF followed by SVD are known as latent semantic indexing (LSI),
# and were first introduced for the analysis of scATAC-seq data by Cusanovich et al. 2015.
brain.atac.tf.cp <- RunTFIDF(brain.atac.tf.cp)
brain.atac.tf.cp <- FindTopFeatures(brain.atac.tf.cp, min.cutoff = 'q50') # top 30% features for UMAP
brain.atac.tf.cp <- RunSVD(brain.atac.tf.cp)

# run non linear dim reduction / use number of dimensions from RNAseq (19)
brain.atac.tf.cp <- RunUMAP(brain.atac.tf.cp, dims = 1:19, reduction = 'lsi')

# plot by Group and by doublet call
pdf(paste0(Sys.Date(),"_UMAP_by_sample_Killi_brain_ATAC_2Cohorts_TF_data.pdf"), height = 5, width = 6)
DimPlot(brain.atac.tf.cp, raster = T, raster.dpi = c(600,600))
dev.off()

# Subset ATAC object by joint cell barcodes
brain.atac.singlets.filt <- brain.atac.singlets[,joint.bcs]
brain.atac.singlets.filt
# An object of class Seurat 
# 148827 features across 26942 samples within 2 assays 
# Active assay: RNA (22135 features, 0 variable features)
# 1 other assay present: ATAC
# 2 dimensional reductions calculated: lsi, umap

# Add TF data as a new assay independent from ATAC/RNA
brain.atac.singlets.filt[["TF"]] <- CreateAssayObject(counts = brain.atac.tf.cp@assays$TF@counts)
brain.atac.singlets.filt         <- NormalizeData(brain.atac.singlets.filt, assay = "TF", normalization.method = "CLR")

# Calculate a UMAP embedding of the HTO data
DefaultAssay(brain.atac.singlets.filt) <- "TF"
brain.atac.singlets.filt <- RunTFIDF(brain.atac.singlets.filt)
brain.atac.singlets.filt <- FindTopFeatures(brain.atac.singlets.filt, min.cutoff = 'q50') # top 30% features for UMAP
brain.atac.singlets.filt <- RunSVD(brain.atac.singlets.filt)

# run non linear dim reduction / use number of dimensions from RNAseq (19)
brain.atac.singlets.filt <- RunUMAP(brain.atac.singlets.filt, dims = 1:19, reduction = 'lsi')

# Check representation
set.seed(123123) # stabilize
P16 = createPalette(16+3,  c("#ff0000", "#00ff00", "#0000ff"))
swatch(P16)
brain.atac.singlets.filt@meta.data$predicted.id <- factor(x = brain.atac.singlets.filt@meta.data$predicted.id, levels = sort(unique(brain.atac.singlets.filt@meta.data$predicted.id)))

pdf(paste0(Sys.Date(),"_UMAP_by_sample_Killi_brain_ATAC_2Cohorts_TF_data_byCellType.pdf"), height = 6, width = 10)
DimPlot(brain.atac.singlets.filt, group.by = 'predicted.id',  cols = as.vector(P16[-c(1:3)]), raster = T, raster.dpi = c(600,600))
dev.off()

pdf(paste0(Sys.Date(),"_UMAP_by_sample_Killi_brain_ATAC_2Cohorts_TF_data_byGroup.pdf"), height = 6, width = 10)
DimPlot(brain.atac.singlets.filt, group.by = 'Group', shuffle = T, cols = c("deeppink4","deepskyblue4","deeppink","deepskyblue"), raster = T, raster.dpi = c(600,600))
dev.off()


##########################################################################################
######### go back to the representation based on accessibility

# Subset ATAC object by joint cell barcodes
brain.atac.singlets.filt <- brain.atac.singlets[,joint.bcs]
brain.atac.singlets.filt
# An object of class Seurat 
# 148827 features across 26942 samples within 2 assays 
# Active assay: RNA (22135 features, 0 variable features)
# 1 other assay present: ATAC
# 2 dimensional reductions calculated: lsi, umap

# Add TF data as a new assay independent from ATAC/RNA
brain.atac.singlets.filt[["TF"]] <- CreateAssayObject(counts = brain.atac.tf.cp@assays$TF@counts)
brain.atac.singlets.filt         <- NormalizeData(brain.atac.singlets.filt, assay = "TF", normalization.method = "CLR")

# p2 <- FeaturePlot(object = brain.atac.singlets.filt,
#                   features = "MA0137.3-STAT1",
#                   min.cutoff = 'q10',
#                   max.cutoff = 'q90',
#                   pt.size = 0.1)
####
DefaultAssay(brain.atac.singlets.filt) <- "TF"

save(brain.atac.singlets.filt, file = paste0(Sys.Date(),"_Killi_brain_ATAC_2Cohorts_FiltNorm_Singlets_with_TF_activity.RData") )

#########################################################################################
# 3. Perform differential analysis
brain.atac.singlets.filt <- SetIdent(brain.atac.singlets.filt, value = 'Group')


brain.atac.singlets.filt$Age_Group <- ifelse(brain.atac.singlets.filt$age_weeks < 7, "Y","O")
brain.atac.singlets.filt           <- SetIdent(brain.atac.singlets.filt, value = 'Age_Group')


# split by cell type
brain.atac.list <- SplitObject(brain.atac.singlets.filt, split.by = "predicted.id")

# nb. of cells per cell type/sample
cell.per.samp.tab <- table(brain.atac.singlets.filt$sample, brain.atac.singlets.filt$predicted.id)

# cell types with at least 10 cells from every each sex/cohort sample
celltype.qc <- colnames(cell.per.samp.tab)[apply(cell.per.samp.tab >= 10, 2, sum) == 8]
celltype.qc <- celltype.qc[-grep("misc",celltype.qc)] # REMOVE "MISC" CELL TYPES
celltype.qc
# [1] "Astrocytes_Radial_Glia"     "Ependymal_cells"            "GABAergic_neurons"          "Granule_Excitatory_Neurons" "Microglia"                 
# [6] "NSPCs"                      "Oligodendrocytes"           "OPCs"                       "PV_interneurons"     

### extract data for cell types that pass the cell number cutoff
brain.atac.list <- brain.atac.list[celltype.qc]

###############################################
#######   DAP analysis   ++++   GRZ     #######
###############################################

# Create list object to receive DESeq2 results
tf.res.list        <- vector(mode = "list", length = length(brain.atac.list))
names(tf.res.list) <- names(brain.atac.list)

# loop over pseudobulk data
for  (i in 1:length(brain.atac.list)) {
  # For sparse data (such as scATAC-seq), we find it is often necessary 
  # to lower the min.pct threshold in FindMarkers() from the default 
  # (0.1, which was designed for scRNA-seq data).
  tf.res.list[[i]]  <- FindMarkers(object = brain.atac.list[[i]],
                                ident.1 = 'O', ident.2 = 'Y',
                                only.pos = F,
                                min.pct = 0.05,
                                logfc.threshold = 0)
}

save(tf.res.list, file = paste0(Sys.Date(),"_findMarkers_DA_TFs_aging.RData"))

######### Make jitter plot of DA peaks #########
## Order by pvalue:
age.results <- lapply(tf.res.list,function(x) {x[order(x$p_val_adj),]})
n        <- sapply(age.results, nrow)
names(n) <- names(age.results)

cols <- list()
xlab <- character(length = length(age.results))
for(i in seq(along = age.results)){
  cols[[i]]      <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), n[i]) # grey60
  ind.sig.i      <- age.results[[i]]$p_val_adj < 0.05
  ind.sig.i.up   <- bitAnd(age.results[[i]]$p_val_adj < 0.05, age.results[[i]]$avg_log2FC >0)>0
  ind.sig.i.down <- bitAnd(age.results[[i]]$p_val_adj < 0.05, age.results[[i]]$avg_log2FC <0)>0
  
  cols[[i]][ind.sig.i.up]   <- "#CC3333"
  cols[[i]][ind.sig.i.down] <- "#333399"
  xlab[i] <- paste(names(age.results)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}
names(cols) <- names(age.results)

pdf(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_snATAC_FindMarkers_TF_with_reg_colors_FDR5.pdf"), width = 6, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 9.5),
     ylim = c(-0.5, 0.5),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change with aging"
)
abline(h = 0)
abline(h = seq(-1, 1, by = 0.25)[-5],
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
         y = rev(age.results[[i]]$avg_log2FC),
         pch = 16,
         col = rev(cols[[i]]),
         bg = rev(cols[[i]]),
         cex = 0.75)
}
axis(1,
     at = 1:9,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-1, 1, by = 0.25))
box()
dev.off()

png(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_snATAC_FindMarkers_TF_with_reg_colors_FDR5.png"), width = 2400, height = 2200, res = 300, units = "px")
par(mar = c(3.1, 4.1, 2, 2))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 9.5),
     ylim = c(-0.5, 0.5),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change  with aging"
)
abline(h = 0)
abline(h = seq(-1, 1, by = 0.25)[-5],
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
         y = rev(age.results[[i]]$avg_log2FC),
         pch = 16,
         col = rev(cols[[i]]),
         bg = rev(cols[[i]]),
         cex = 0.75)
}
axis(1,
     at = 1:9,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-1, 1, by = 0.25))
box()
dev.off()
###########################################################################################

###########################################################################################
# 4. parse across cell types

summary.mtx <- data.frame(matrix(0,0,4))
colnames(summary.mtx) <- c("TF","avg_log2FC","p_val_adj","Cell_Type")

for (i in 1:length(tf.res.list)) {
  
  # extract FDR < 5%
  tmp.mat <- tf.res.list[[i]]
  tmp.mat <- tmp.mat[tmp.mat$p_val_adj < 0.05,]
  tmp.mat$TF <- rownames(tmp.mat)
  
  # extract columns
  tmp.mat <- tmp.mat[,c("TF","avg_log2FC","p_val_adj")]
  
  # populate Cell type  column
  tmp.mat$Cell_Type <- names(tf.res.list)[i]
  
  summary.mtx <- rbind(summary.mtx,tmp.mat)
}

write.table(summary.mtx, file = paste0(Sys.Date(),"_CellRangerTF_findMarker_DA_footprints_by_CellType_GRZ_FDR5.txt"), sep = "\t", quote = F, row.names = F)

# Reccurrent TF in more than half the cell types
concat <- function (x) {
  paste0(x, collapse = ",")
}

summary.mtx$sign <- sign(summary.mtx$avg_log2FC)
my.num.grz   <- aggregate(summary.mtx$Cell_Type, by = list("TF_bind" = summary.mtx$TF), FUN = length)
colnames(my.num.grz)[2] <- "N_CellType"
my.signs.grz <- aggregate(summary.mtx$sign, by = list("TF_bind" = summary.mtx$TF), FUN = concat)
colnames(my.signs.grz)[2] <- "Signs"
my.sum.grz <- merge(my.num.grz,my.signs.grz)

# recurrent perturbed regulons
my.sum.grz[my.sum.grz$N_CellType == 9,]
#                   TF_bind N_CellType             Signs
# 5            MA0007.3-Ar          9 1,1,1,1,1,1,1,1,1
# 83        MA0113.3-NR3C1          9 1,1,1,1,1,1,1,1,1
# 86       MA0116.1-Znf423          9 1,1,1,1,1,1,1,1,1
# 93        MA0131.2-HINFP          9 1,1,1,1,1,1,1,1,1
# 163       MA0497.1-MEF2C          9 1,1,1,1,1,1,1,1,1
# 191        MA0526.2-USF2          9 1,1,1,1,1,1,1,1,1
# 234        MA0631.1-Six3          9 1,1,1,1,1,1,1,1,1
# 252        MA0649.1-HEY2          9 1,1,1,1,1,1,1,1,1
# 255        MA0652.1-IRF8          9 1,1,1,1,1,1,1,1,1
# 256        MA0653.1-IRF9          9 1,1,1,1,1,1,1,1,1
# 260       MA0657.1-KLF13          9 1,1,1,1,1,1,1,1,1
# 290        MA0687.1-SPIC          9 1,1,1,1,1,1,1,1,1
# 330       MA0727.1-NR3C2          9 1,1,1,1,1,1,1,1,1
# 332        MA0729.1-RARA          9 1,1,1,1,1,1,1,1,1
# 338       MA0735.1-GLIS1          9 1,1,1,1,1,1,1,1,1
# 418       MA0816.1-Ascl2          9 1,1,1,1,1,1,1,1,1
# 421       MA0819.1-CLOCK          9 1,1,1,1,1,1,1,1,1
# 453       MA0851.1-Foxj3          9 1,1,1,1,1,1,1,1,1
# 459        MA0857.1-Rarb          9 1,1,1,1,1,1,1,1,1
# 525        MA1107.1-KLF9          9 1,1,1,1,1,1,1,1,1
# 566 MA1148.1-PPARA::RXRA          9 1,1,1,1,1,1,1,1,1
# 575        MA1419.1-IRF4          9 1,1,1,1,1,1,1,1,1
# 576        MA1420.1-IRF5          9 1,1,1,1,1,1,1,1,1

pdf(paste0(Sys.Date(),"_CellRangerTF_findMarker_DA_footprints_CellTypeSharing_Boxplot_GRZ.pdf"), height = 5, width = 3.5)
boxplot(list("GRZ"= my.sum.grz$N_CellType), col = "goldenrod1", main = "DecoupleR Regulon Sharing", ylab = "Number of cell types", outline = F, ylim = c(0,10))
beeswarm::beeswarm(my.sum.grz$N_CellType[my.sum.grz$N_CellType == 9], method = "compactswarm",corral = "gutter", add = T, pch = 16)
text(1.1, 9.5, "nr3c1"          , pos = 4, cex = 0.5)
dev.off()

################################################
# get summary plot
library(ggplot2) 
library(scales) 
theme_set(theme_bw())

my.all.cell.types <- my.sum.grz$TF_bind[my.sum.grz$N_CellType == 9]

summary.mtx.sub <- summary.mtx[summary.mtx$TF %in% my.all.cell.types,]
summary.mtx.sub$minusLog10padj <- -log10(summary.mtx.sub$p_val_adj)

av.change <- data.frame(aggregate(summary.mtx.sub$avg_log2FC,by = list("TF"=summary.mtx.sub$TF), FUN = "mean"))
av.change.srt <- sort(av.change$x, decreasing = T, index.return = T)

#### 
my.max <- max(summary.mtx.sub$avg_log2FC)
my.min <- -my.max
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
summary.mtx.sub$TF  <- factor(summary.mtx.sub$TF, levels = rev(av.change$TF[av.change.srt$ix]))

pdf(paste0(Sys.Date(),"JASPAR_Footprints_snATAC_DA_KilliBrainAging_all_FDR5.pdf"), height = 8, width=5.5)
my.plot <- ggplot(summary.mtx.sub,aes(x=Cell_Type,y=TF,colour=avg_log2FC,size=minusLog10padj))+ theme_bw() + geom_point(shape = 16)
my.plot <- my.plot + ggtitle("JASPAR footprints") + labs(x = "Killifish Brain Aging", y = "Aging DA snATAC footprints (FDR < 5%)")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(my.min, my.max))
my.plot <- my.plot + scale_x_discrete(guide = guide_axis(angle = 45))
print(my.plot)
dev.off()

# just top 10
summary.mtx.sub10 <- summary.mtx.sub[summary.mtx.sub$TF %in% av.change$TF[av.change.srt$ix[1:10]],]

pdf(paste0(Sys.Date(),"JASPAR_Footprints_snATAC_DA_KilliBrainAging_Top10_FDR5.pdf"), height = 5, width=5.5)
my.plot <- ggplot(summary.mtx.sub10,aes(x=Cell_Type,y=TF,colour=avg_log2FC,size=minusLog10padj))+ theme_bw() + geom_point(shape = 16)
my.plot <- my.plot + ggtitle("Top 10 DA JASPAR footprints") + labs(x = "Killifish Brain Aging", y = "Aging DA snATAC footprints (FDR < 5%)")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(my.min, my.max))
my.plot <- my.plot + scale_x_discrete(guide = guide_axis(angle = 45))
print(my.plot)
dev.off()
# differential.activity["MA0113.3",]
###########################################################################################
