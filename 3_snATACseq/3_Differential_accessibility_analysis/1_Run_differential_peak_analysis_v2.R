setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/snATAC_Brain_Aging_Meta/Downstream_Analyses/Differential_Accessiblity')
options(stringsAsFactors = F)

#### Packages
library(Seurat)        # single cell general package
library(Signac)        # scATAC processing
library(sctransform)      # 
library("singleCellTK")   # 

library(GenomeInfoDb)  # for genome info
library(GenomicRanges) # for genome info
library(seqinr)        # for genome info

library(Polychrome)
library(ggplot2)

library('muscat')         # 
library('DESeq2')         # 
library('sva')            # 
library('limma')          # 

library('bitops')          # 

theme_set(theme_bw())   

# 2023-08-24
# Process scATAC cohorts for differential accessiblity analysis

# 2024-04012
# clean up analysis to match RNA analysis

#########################################################################################
# 0. Load annotated object
load("../../Signac/2023-08-23_Killi_brain_ATAC_2Cohorts_FiltNorm_Singlets_ANNOTATEDfromRNA.RData")

brain.atac.singlets
# An object of class Seurat 
# 148827 features across 27171 samples within 2 assays 
# Active assay: RNA (22135 features, 0 variable features)
# 1 other assay present: ATAC
# 2 dimensional reductions calculated: lsi, umap

# change back to working with peaks instead of gene activities
DefaultAssay(brain.atac.singlets) <- 'ATAC'
brain.atac.singlets <- SetIdent(brain.atac.singlets, value = 'Group')

# convert to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
killi.atac.sce <- as.SingleCellExperiment(brain.atac.singlets)
###########################################################################################

###########################################################################################
# 1. Find differentially accessible peaks

killi.atac.sce.cl <- prepSCE(killi.atac.sce, 
                             kid   = "predicted.id"   ,  # cell population assignments
                             gid   = "Group"          ,  # group IDs (ctrl/stim)
                             sid   = "sample"         ,  # sample IDs (ctrl/stim.1234)
                             drop  = TRUE             )  # drop all other colData columns


# store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk          <- length(kids <- levels(killi.atac.sce.cl$cluster_id))
ns          <- length(sids <- levels(killi.atac.sce.cl$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
t(table(killi.atac.sce.cl$cluster_id, killi.atac.sce.cl$sample_id))
#     Astrocytes_Radial_Glia Ependymal_cells Erythrocytes GABAergic_neurons Granule_Excitatory_Neurons Microglia Neurons_misc_1
# OF1                    109              34            3                76                       1615        69             16
# OF2                    147              29           19                61                       1832       131              9
# OM1                     27              13            3                31                        661        30              3
# OM2                    140              33           15                81                       2024        85             10
# YF1                    154              28           11               128                        875        33             22
# YF2                    205              33           13               160                       1387        43             23
# YM1                    115              27            2               143                       1077        37             18
# YM2                    187              32           19               135                       1196       138             23
# 
#     Neurons_misc_2 Neurons_misc_3 Neurons_misc_4 NSPCs Oligodendrocytes OPCs Purkinje_cells PV_interneurons Vascular_smooth_muscle_cells
# OF1            136             91           1192    53              219   39              6              63                           33
# OF2            102             48            897    36              251   26              5              45                           28
# OM1             57             27            431    24               99   13              0              29                            6
# OM2            126            100           1155    64              311   37              3              53                           31
# YF1             40             64           1482    90              121   37              9              66                           12
# YF2            141             66           1641    80              167   52             12              61                           22
# YM1             31             46           1282    82               94   40              9              69                           29
# YM2            125             68           1278    67              229   39             13              63                           38

# Aggregation of single-cell to pseudobulk data
pb <- aggregateData(killi.atac.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# one sheet per subpopulation
assayNames(pb)
# [1] "Astrocytes_Radial_Glia"       "Ependymal_cells"              "Erythrocytes"                 "GABAergic_neurons"           
# [5] "Granule_Excitatory_Neurons"   "Microglia"                    "Neurons_misc_1"               "Neurons_misc_2"              
# [9] "Neurons_misc_3"               "Neurons_misc_4"               "NSPCs"                        "Oligodendrocytes"            
# [13] "OPCs"                         "Purkinje_cells"               "PV_interneurons"              "Vascular_smooth_muscle_cells"

# Pseudobulk-level MDS plot
pb_mds <- pbMDS(pb)

# output global MDS
pdf(paste0(Sys.Date(),"_10x_Killi_Brain_Aging_scATAC_Muscat_PB_MDS.pdf"))
pb_mds
dev.off()

# nb. of cells per cluster-sample
cell.per.samp.tab <- t(table(killi.atac.sce.cl$cluster_id, killi.atac.sce.cl$sample_id))

# cell types with at least 10 cells from every each sex/cohort sample
celltype.qc <- colnames(cell.per.samp.tab)[apply(cell.per.samp.tab >= 10, 2, sum) == 8]

# REMOVE "MISC" CELL TYPES
celltype.qc <- celltype.qc[-grep("misc",celltype.qc)]
celltype.qc
# keep:
# [1] "Astrocytes_Radial_Glia"     "Ependymal_cells"            "GABAergic_neurons"          "Granule_Excitatory_Neurons" "Microglia"                 
# [6] "NSPCs"                      "Oligodendrocytes"           "OPCs"                       "PV_interneurons"    

### extract pseudobulk information for samples that pass the cell number cutoff
counts.pb <- pb@assays@data[celltype.qc]

#### save counts
save(counts.pb, celltype.qc,
     file = paste0(Sys.Date(),"_muscat_PB_GRZ_ATAC_objects_QC_Clean.RData"))
###############################################################################################


##############################################################################################
# 2. Use SVA to clean up batch effects on expression and DEseq2 for DA analysis

# clean up memory and reload only muscat PBs
load('2024-04-15_muscat_PB_GRZ_ATAC_objects_QC_Clean.RData')

# will run SVA to clean up noise
# run for the cell types with at least 10 cells from every each sample

###############################################
#######   DEG analysis   ++++   GRZ     #######
###############################################

# import metadata and order it
my.grz.meta           <- data.frame(matrix(0,8,6))
colnames(my.grz.meta) <- c("SampleID","Age_Group","Age_Weeks","Batch","Sex","Group")
my.grz.meta$SampleID  <- c("YF1", "YF2", "YM1", "YM2","OF1", "OF2", "OM1", "OM2" )
my.grz.meta$Age_Group <- factor(c(rep("Y",4),rep("O",4)) , levels = c("Y","O"))
my.grz.meta$Age_Weeks <- c(6.29, 6.14, 6.29, 6.14, 15.14, 16.23, 15.14,16.57)
my.grz.meta$Batch     <- rep(c("c1","c2"),4)
my.grz.meta$Sex       <- factor(c(rep("F",2),rep("M",2),rep("F",2),rep("M",2)) , levels = c("F","M"))
my.grz.meta$Group     <- paste0(my.grz.meta$Age_Group,my.grz.meta$Sex)
rownames(my.grz.meta) <- my.grz.meta$SampleID

# reorder count tables in sensical order
for  (i in 1:length(counts.pb)) {
  counts.pb[[i]] <- counts.pb[[i]][,my.grz.meta$SampleID]
}

# Create list object to receive clean SVA counts
sva.cts.grz        <- vector(mode = "list", length = length(counts.pb))
names(sva.cts.grz) <- names(counts.pb)

# Create list object to receive VST normalized counts
vst.cts.grz        <- vector(mode = "list", length = length(counts.pb))
names(vst.cts.grz) <- names(counts.pb)

# Create list object to receive DESeq2 results
deseq.res.list.grz        <- vector(mode = "list", length = length(counts.pb))
names(deseq.res.list.grz) <- names(counts.pb)

# loop over pseudobulk data
for  (i in 1:length(counts.pb)) {
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_snATAC_Pseudobulk_GRZ_",names(counts.pb)[[i]])
  
  ###################################
  #######       Run SVA      #######
  
  # build design matrix
  sva.dataDesign = data.frame( row.names = my.grz.meta$SampleID , 
                               sex       = my.grz.meta$Sex      ,
                               age       = my.grz.meta$Age_Weeks,
                               batch     = my.grz.meta$Batch)
  
  # Set null and alternative models (ignore batch)
  mod1    = model.matrix(~ sex + age + batch, data = sva.dataDesign)
  n.sv.be = num.sv(counts.pb[[i]], mod1, method="be") # 
  
  # apply SVAseq algortihm
  my.svseq = svaseq(as.matrix(counts.pb[[i]]), mod1, n.sv=n.sv.be, constant = 0.1)
  
  # remove RIN and SV, preserve age and sex
  my.clean <- removeBatchEffect(log2(counts.pb[[i]] + 0.1), 
                                batch      = my.grz.meta$Batch, 
                                covariates = cbind(my.svseq$sv),
                                design     = mod1[,1:3])
  
  # delog and round data for DEseq2 processing
  my.filtered.sva <- round(2^my.clean-0.1)
  
  # keep only robustly expressed genes
  sva.cts.grz[[i]] <- my.filtered.sva
  
  # legend
  my.cols  <- rep("",nrow(my.grz.meta))
  my.cols[my.grz.meta$Group %in% "YF"] <- "deeppink"
  my.cols[my.grz.meta$Group %in% "OF"] <- "deeppink4"
  my.cols[my.grz.meta$Group %in% "YM"] <- "deepskyblue"
  my.cols[my.grz.meta$Group %in% "OM"] <- "deepskyblue4"
  
  my.pch  <- rep(0,nrow(my.grz.meta))
  my.pch[my.grz.meta$Group %in% "YF"] <- 16
  my.pch[my.grz.meta$Group %in% "OF"] <- 15
  my.pch[my.grz.meta$Group %in% "YM"] <- 16
  my.pch[my.grz.meta$Group %in% "OM"] <- 15
  
  # get matrix using age as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = sva.cts.grz[[i]],
                                colData   = my.grz.meta,
                                design    = ~ Age_Weeks + Sex)
  
  # run DESeq normalizations and export results
  dds.deseq <- DESeq(dds)
  
  # plot dispersion
  my.disp.out <- paste(my.outprefix,"_dispersion_plot.pdf")
  
  pdf(my.disp.out)
  plotDispEsts(dds.deseq)
  dev.off()
  
  # get DESeq2 normalized expression value
  vst.cts.grz[[i]] <- getVarianceStabilizedData(dds.deseq)
  
  # MDS analysis
  mds.result <- cmdscale(1-cor(vst.cts.grz[[i]],method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  x <- mds.result[, 1]
  y <- mds.result[, 2]
  
  pdf(paste0(my.outprefix,"_MDS_plot.pdf"))
  plot(x, y,
       xlab = "MDS dimension 1", ylab = "MDS dimension 2",
       main= paste0(names(counts.pb)[[i]]," (MDS)"),
       cex=3, col = my.cols, pch = my.pch,
       cex.lab = 1.25,
       cex.axis = 1.25, las = 1)
  dev.off()
  
  # extract peaks significance by DEseq2
  res.age <- results(dds.deseq, name = "Age_Weeks") # FC per week
  
  # exclude peaks with NA FDR value
  res.age <- res.age[!is.na(res.age$padj),]
  
  # store results
  deseq.res.list.grz[[i]]       <- data.frame(res.age)

  ### get sex dimorphic changes at FDR5
  peaks.age <- rownames(res.age)[res.age$padj < 0.05]
  my.num.age <- length(peaks.age)
  
  if (my.num.age > 2) {
    # heatmap drawing - only if there is at least 2 gene
    my.heatmap.out <- paste0(my.outprefix,"_AGING_Heatmap_DA_ATC_peaks_FDR5.pdf")
    
    pdf(my.heatmap.out, onefile = F, height = 10, width = 10)
    my.heatmap.title <- paste0(names(counts.pb)[[i]], " aging significant (FDR<5%), ", my.num.age, " genes")
    pheatmap::pheatmap(vst.cts.grz[[i]][peaks.age,],
                       cluster_cols = F,
                       cluster_rows = T,
                       colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                       show_rownames = F, scale="row",
                       main = my.heatmap.title,
                       cellwidth = 15,
                       border    = NA,
                       cellheight = 0.15 )
    dev.off()
  }
  
  # output result tables of combined analysis to text files
  my.out.ct.mat <- paste0(my.outprefix,"_AGING_VST_log2_counts_matrix.txt")
  write.table(vst.cts.grz[[i]], file = my.out.ct.mat , sep = "\t" , row.names = T, quote = F)
  
  my.out.stats.age <- paste0(my.outprefix,"_AGING_all_peaks_statistics.txt")
  write.table(deseq.res.list.grz[[i]], file = my.out.stats.age , sep = "\t" , row.names = T, quote = F)
  
  my.out.fdr5.age <- paste0(my.outprefix,"_AGING_FDR5_peaks_statistics.txt")
  write.table(deseq.res.list.grz[[i]][peaks.age,], file = my.out.fdr5.age, sep = "\t" , row.names = T, quote = F)
  
}

# save R object with all DEseq2 results
my.rdata.age <- paste0(Sys.Date(),"_pseudobulk_killi_cell_types_snATAC_AGING_GRZ_DEseq2_objects.RData")
save(deseq.res.list.grz, file = my.rdata.age)

my.vst.age <- paste0(Sys.Date(),"_pseudobulk_killi_cell_types_AGING_GRZ_VST_data_objects.RData")
save(vst.cts.grz, file = my.vst.age)

######### Make jitter plot of DA peaks #########
## Order by pvalue:
age.results <- lapply(deseq.res.list.grz,function(x) {x[order(x$padj),]})
n        <- sapply(age.results, nrow)
names(n) <- names(age.results)

cols <- list()
xlab <- character(length = length(age.results))
for(i in seq(along = age.results)){
  cols[[i]]      <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), n[i]) # grey60
  ind.sig.i      <- age.results[[i]]$padj < 0.05
  ind.sig.i.up   <- bitAnd(age.results[[i]]$padj < 0.05, age.results[[i]]$log2FoldChange >0)>0
  ind.sig.i.down <- bitAnd(age.results[[i]]$padj < 0.05, age.results[[i]]$log2FoldChange <0)>0
  
  cols[[i]][ind.sig.i.up]   <- "#CC3333"
  cols[[i]][ind.sig.i.down] <- "#333399"
  xlab[i] <- paste(names(age.results)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}
names(cols) <- names(age.results)

pdf(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_snATAC_DESeq2_with_reg_colors_FDR5.pdf"), width = 6, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 9.5),
     ylim = c(-0.75, 0.75),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change per week of life"
)
abline(h = 0)
abline(h = seq(-1, 1, by = 0.25)[-5],
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
         y = rev(age.results[[i]]$log2FoldChange),
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
     at = seq(-1, 1, by = 0.5))
box()
dev.off()

png(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_snATAC_DESeq2_with_reg_colors_FDR5.png"), width = 2400, height = 2200, res = 300, units = "px")
par(mar = c(3.1, 4.1, 2, 2))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 9.5),
     ylim = c(-0.75, 0.75),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change per week of life"
)
abline(h = 0)
abline(h = seq(-1, 1, by = 0.25)[-5],
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
         y = rev(age.results[[i]]$log2FoldChange),
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
     at = seq(-1, 1, by = 0.5))
box()
dev.off()


######### Make jitter plot of DA peaks - FDR10% #########
## Order by pvalue:
age.results <- lapply(deseq.res.list.grz,function(x) {x[order(x$padj),]})
n        <- sapply(age.results, nrow)
names(n) <- names(age.results)

cols <- list()
xlab <- character(length = length(age.results))
for(i in seq(along = age.results)){
  cols[[i]]      <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), n[i]) # grey60
  ind.sig.i      <- age.results[[i]]$padj < 0.1
  ind.sig.i.up   <- bitAnd(age.results[[i]]$padj < 0.1, age.results[[i]]$log2FoldChange >0)>0
  ind.sig.i.down <- bitAnd(age.results[[i]]$padj < 0.1, age.results[[i]]$log2FoldChange <0)>0
  
  cols[[i]][ind.sig.i.up]   <- "#CC3333"
  cols[[i]][ind.sig.i.down] <- "#333399"
  xlab[i] <- paste(names(age.results)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}
names(cols) <- names(age.results)

pdf(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_snATAC_DESeq2_with_reg_colors_FDR10.pdf"), width = 6, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 9.5),
     ylim = c(-0.75, 0.75),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change per week of life"
)
abline(h = 0)
abline(h = seq(-1, 1, by = 0.25)[-5],
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
         y = rev(age.results[[i]]$log2FoldChange),
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
     at = seq(-1, 1, by = 0.5))
box()
dev.off()

png(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_snATAC_DESeq2_with_reg_colors_FDR10.png"), width = 2400, height = 2200, res = 300, units = "px")
par(mar = c(3.1, 4.1, 2, 2))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 9.5),
     ylim = c(-0.75, 0.75),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change per week of life"
)
abline(h = 0)
abline(h = seq(-1, 1, by = 0.25)[-5],
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
         y = rev(age.results[[i]]$log2FoldChange),
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
     at = seq(-1, 1, by = 0.5))
box()
dev.off()

######### Make jitter plot of DA peaks - FDR15% #########
## Order by pvalue:
age.results <- lapply(deseq.res.list.grz,function(x) {x[order(x$padj),]})
n        <- sapply(age.results, nrow)
names(n) <- names(age.results)

cols <- list()
xlab <- character(length = length(age.results))
for(i in seq(along = age.results)){
  cols[[i]]      <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), n[i]) # grey60
  ind.sig.i      <- age.results[[i]]$padj < 0.15
  ind.sig.i.up   <- bitAnd(age.results[[i]]$padj < 0.15, age.results[[i]]$log2FoldChange >0)>0
  ind.sig.i.down <- bitAnd(age.results[[i]]$padj < 0.15, age.results[[i]]$log2FoldChange <0)>0
  
  cols[[i]][ind.sig.i.up]   <- "#CC3333"
  cols[[i]][ind.sig.i.down] <- "#333399"
  xlab[i] <- paste(names(age.results)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}
names(cols) <- names(age.results)

pdf(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_snATAC_DESeq2_with_reg_colors_FDR15.pdf"), width = 6, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 9.5),
     ylim = c(-0.75, 0.75),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change per week of life"
)
abline(h = 0)
abline(h = seq(-1, 1, by = 0.25)[-5],
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
         y = rev(age.results[[i]]$log2FoldChange),
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
     at = seq(-1, 1, by = 0.5))
box()
dev.off()

png(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_snATAC_DESeq2_with_reg_colors_FDR15.png"), width = 2400, height = 2200, res = 300, units = "px")
par(mar = c(3.1, 4.1, 2, 2))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 9.5),
     ylim = c(-0.75, 0.75),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change per week of life"
)
abline(h = 0)
abline(h = seq(-1, 1, by = 0.25)[-5],
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
         y = rev(age.results[[i]]$log2FoldChange),
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
     at = seq(-1, 1, by = 0.5))
box()
dev.off()
###############################################################

###############################################################
options(java.parameters = "-Xmx16g" )
require(openxlsx)

deseq.out <-  paste0(Sys.Date(),"_GRZ_snATAC_DESeq2_Aging_Results.xlsx")

write.xlsx(deseq.res.list.grz, rowNames = TRUE, file = deseq.out)
###############################################################

#######################
sink(file = paste(Sys.Date(),"_Seurat_Signac_MuscatDEseq2_PB_scATAC_session_Info.txt", sep =""))
sessionInfo()
sink()