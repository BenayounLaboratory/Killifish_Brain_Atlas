setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Differential_Expression_Age/DGE_Analysis')
options(stringsAsFactors = F)

#### Packages
library('Seurat')         # 
library(sctransform)      # 
library("singleCellTK")   # 

library('muscat')         # 
library('DESeq2')         # 
library('sva')            # 
library('limma')          # 

library(ggplot2)          # 
library(scales)           # 
library("bitops")         # 
library(Vennerable)       # 
library(data.table)       #

library(ComplexHeatmap)   #
library(circlize)         #



theme_set(theme_bw())   

# 2024-02-16
# Process scRNAseq brain aging cohorts for differential gene analysis

# 2024-02-21
# rerun with same set of cell types for both strains (only run cell types passing QC for both strains)

# 2024-02-26
# clean up output to wrap text

# 2023-02-28
# Only run with well defined cell types, since neuron misc is difficult to interpret

# 2024-03-20
# Separate DGE analysis and GSEA for ease of reading

###############################################################################################
# 0. preprocess Seurat object for use with muscat

# Import final annotation
load('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Atlas_Analysis/2023-08-23_Seurat_object_with_Manual_annotation_FINAL.RData')
killi.brain.clean
# An object of class Seurat 
# 21160 features across 209939 samples within 1 assay 
# Active assay: RNA (21160 features, 5000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony

# Add in a sample ID
killi.brain.clean@meta.data$SampleID <- paste0(killi.brain.clean@meta.data$Group,"_",killi.brain.clean@meta.data$Batch)

# subset per strains
killi.brain.grz <- subset(killi.brain.clean, subset = Strain %in% "GRZ")    # 104,665 cells
killi.brain.zmz <- subset(killi.brain.clean, subset = Strain %in% "ZMZ")    # 105,274 cells
save(killi.brain.grz, killi.brain.zmz, 
     file = paste(Sys.Date(),"Seurat_objects_SPLIT_PER_STRAIN_withSampleID_with_Manual_annotation_FINAL.RData",sep = "_"))


#### Make covariate table
grz.meta <- unique(killi.brain.grz@meta.data[,c("SampleID", "Age_Group" , "Age_weeks", "Group", "Sex" )])
zmz.meta <- unique(killi.brain.zmz@meta.data[,c("SampleID", "Age_Group" , "Age_weeks", "Group", "Sex" )])

write.table(grz.meta, file = paste0(Sys.Date(),"_GRZ_sample_metadata_table.txt"), sep = "\t", row.names = F, quote = F)
write.table(zmz.meta, file = paste0(Sys.Date(),"_ZMZ_sample_metadata_table.txt"), sep = "\t", row.names = F, quote = F)


# bring RNA as main assay fro processing
DefaultAssay(killi.brain.grz) <- "RNA"
DefaultAssay(killi.brain.zmz) <- "RNA"

# convert to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
killi.brain.clean.GRZ.sce <- as.SingleCellExperiment(killi.brain.grz)
killi.brain.clean.ZMZ.sce <- as.SingleCellExperiment(killi.brain.zmz)
save(killi.brain.clean.GRZ.sce, killi.brain.clean.ZMZ.sce, 
     file = paste(Sys.Date(),"SingleCellExperimnents_objects_PER_STRAIN.RData",sep = "_"))

rm(killi.brain.clean) # free up some memory
rm(killi.brain.grz, killi.brain.zmz) # free up some memory
###############################################################################################


###############################################################################################
# 1. Run muscat for pseudobulking and extraction of samples

# Clean workspace and reload necessary objects
load('2024-02-16_SingleCellExperimnents_objects_PER_STRAIN.RData')

###############################################
####### Data preparation   ++++   GRZ   #######
###############################################

brain.grz.sce.cl <- prepSCE(killi.brain.clean.GRZ.sce, 
                            kid    = "Cell_Identity",  # population assignments
                            gid    = "Sex"          ,  # group IDs (ctrl/stim)
                            sid    = "SampleID"     ,  # sample IDs (ctrl/stim.1234)
                            drop   = TRUE           )  # drop all other colData columns

# store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk  <- length(kids <- levels(brain.grz.sce.cl$cluster_id))
ns  <- length(sids <- levels(brain.grz.sce.cl$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
t(table(brain.grz.sce.cl$cluster_id, brain.grz.sce.cl$sample_id))

# Aggregation of single-cell to pseudobulk data
pb.grz <- aggregateData(brain.grz.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# one list item per cell type
assayNames(pb.grz)
# [1] "Astrocytes_Radial_Glia"       "Ependymal_cells"              "Erythrocytes"                 "Microglia"                   
# [5] "NSPCs"                        "Oligodendrocytes"             "OPCs"                         "Vascular_smooth_muscle_cells"
# [9] "GABAergic_neurons"            "Granule_Excitatory_Neurons"   "Neurons_misc_1"               "Neurons_misc_2"              
# [13] "Neurons_misc_3"               "Neurons_misc_4"               "Purkinje_cells"               "PV_interneurons"    


# Number of cells in each sample and cell type
cell.per.samp.tab.grz <- t(table(brain.grz.sce.cl$cluster_id, brain.grz.sce.cl$sample_id))

# extract pseudobulk information
counts.pb.tmp.grz <- pb.grz@assays@data

# some genes have "gene-" in their name: remove for processing
for (i in 1:length(counts.pb.tmp.grz)) {
  rownames(counts.pb.tmp.grz[[i]]) <- gsub("gene-", "",rownames(counts.pb.tmp.grz[[i]]))
}

# get the genes with no reads in at least half the samples out, they mess up the algorithm
for (i in 1:length(counts.pb.tmp.grz)) {
  my.good <- which(apply(counts.pb.tmp.grz[[i]]>0, 1, sum) >= nrow(cell.per.samp.tab.grz)/2) # see deseq2 vignette, need to remove too low genes
  counts.pb.tmp.grz[[i]] <- counts.pb.tmp.grz[[i]][my.good,]
}


###############################################
####### Data preparation   ++++   ZMZ   #######
###############################################

brain.zmz.sce.cl <- prepSCE(killi.brain.clean.ZMZ.sce, 
                            kid    = "Cell_Identity",  # population assignments
                            gid    = "Sex"          ,  # group IDs (ctrl/stim)
                            sid    = "SampleID"     ,  # sample IDs (ctrl/stim.1234)
                            drop   = TRUE           )  # drop all other colData columns

# store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk  <- length(kids <- levels(brain.zmz.sce.cl$cluster_id))
ns  <- length(sids <- levels(brain.zmz.sce.cl$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
t(table(brain.zmz.sce.cl$cluster_id, brain.zmz.sce.cl$sample_id))

# Aggregation of single-cell to pseudobulk data
pb.zmz <- aggregateData(brain.zmz.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# one list item per cell type
assayNames(pb.zmz)
# [1] "Astrocytes_Radial_Glia"       "Ependymal_cells"              "Erythrocytes"                 "Microglia"                   
# [5] "NSPCs"                        "Oligodendrocytes"             "OPCs"                         "Vascular_smooth_muscle_cells"
# [9] "GABAergic_neurons"            "Granule_Excitatory_Neurons"   "Neurons_misc_1"               "Neurons_misc_2"              
# [13] "Neurons_misc_3"               "Neurons_misc_4"               "Purkinje_cells"               "PV_interneurons"    

# Number of cells in each sample and cell type
cell.per.samp.tab.zmz <- t(table(brain.zmz.sce.cl$cluster_id, brain.zmz.sce.cl$sample_id))

# extract pseudobulk information for samples that pass the cell number cutoff
counts.pb.tmp.zmz <- pb.zmz@assays@data


# some genes have "gene-" in their name: remove for processing
for (i in 1:length(counts.pb.tmp.zmz)) {
  rownames(counts.pb.tmp.zmz[[i]]) <- gsub("gene-", "",rownames(counts.pb.tmp.zmz[[i]]))
}

# get the genes with no reads in at least half the samples out, they mess up the algorithm
for (i in 1:length(counts.pb.tmp.zmz)) {
  my.good <- which(apply(counts.pb.tmp.zmz[[i]]>0, 1, sum) >= nrow(cell.per.samp.tab.zmz)/2) # see deseq2 vignette, need to remove too low genes
  counts.pb.tmp.zmz[[i]] <- counts.pb.tmp.zmz[[i]][my.good,]
}
########################################################


####################################################################
#######   Data preparation   ++++   GRZ/ZMZ QC Cell types    #######
####################################################################

# cell types with at least 10 cells in all samples (20 samples for GRZ)
grz.celltype.qc <- colnames(cell.per.samp.tab.grz)[colSums(cell.per.samp.tab.grz  >= 10) == nrow(cell.per.samp.tab.grz)]
grz.celltype.qc
# [1] "Astrocytes_Radial_Glia"       "Microglia"                    "NSPCs"                        "Oligodendrocytes"            
# [5] "OPCs"                         "Vascular_smooth_muscle_cells" "GABAergic_neurons"            "Granule_Excitatory_Neurons"  
# [9] "Neurons_misc_3"               "Neurons_misc_4"               "PV_interneurons"             

# cell types with at least 10 cells in all samples (24 samples for ZMZ)
zmz.celltype.qc <- colnames(cell.per.samp.tab.zmz)[colSums(cell.per.samp.tab.zmz  >= 10) == nrow(cell.per.samp.tab.zmz)]
zmz.celltype.qc
# [1] "Astrocytes_Radial_Glia"     "Ependymal_cells"            "Microglia"                  "NSPCs"                      "Oligodendrocytes"          
# [6] "OPCs"                       "GABAergic_neurons"          "Granule_Excitatory_Neurons" "Neurons_misc_3"             "Neurons_misc_4"            
# [11] "PV_interneurons"           



##### QC cell types
my.criteria <- list("GRZ QC Cells"    = grz.celltype.qc,  # "Vascular_smooth_muscle_cells"
                    "ZMZ QC Cells"    = zmz.celltype.qc)  # "Ependymal_cells"
my.Venn <- Venn(my.criteria)

pdf(paste0(Sys.Date(),"_QC_Cell_Types_Overlap_ZMZ_GRZ.pdf"))
plot(my.Venn, doWeights=T)
dev.off()

#### study cell types passing QC in both strains
both.celltype.qc   <- sort(intersect(grz.celltype.qc, zmz.celltype.qc))

# Exclude poorly defined Neurons_misc cell types
both.celltype.qc   <- both.celltype.qc[-grep("Neurons_misc", both.celltype.qc)] 
both.celltype.qc
# [1] "Astrocytes_Radial_Glia"     "GABAergic_neurons"          "Granule_Excitatory_Neurons" "Microglia"                  "NSPCs"                           
# [6] "Oligodendrocytes"           "OPCs"                       "PV_interneurons"       

# extract pseudobulk information for samples that pass the cell number cutoff
counts.pb.grz <- counts.pb.tmp.grz[both.celltype.qc]

# extract pseudobulk information for samples that pass the cell number cutoff
counts.pb.zmz <- counts.pb.tmp.zmz[both.celltype.qc]

#### save counts
save(counts.pb.grz, counts.pb.zmz, 
     grz.celltype.qc, zmz.celltype.qc,
     file = paste0(Sys.Date(),"_muscat_PB_GRZ_ZMZ_objects_QC_Clean.RData"))
###############################################################################################


##############################################################################################
# 2. Use SVA to clean up batch effects on expression and DEseq2 for DE analysis

# clean up memory and reload only muscat PBs
load('2024-02-28_muscat_PB_GRZ_ZMZ_objects_QC_Clean.RData')

# will run SVA to clean up noise
# run for the cell types with at least 10 cells from every each sample

###############################################
#######   DEG analysis   ++++   GRZ     #######
###############################################

# import metadata and order it
my.grz.meta           <- read.csv("2024-02-19_GRZ_sample_metadata_table.txt", sep = "\t")
my.grz.meta$Age_Group <- factor(my.grz.meta$Age_Group , levels = c("Y", "M","O"))
my.grz.meta           <- setorder(my.grz.meta,Age_Group, SampleID)
my.grz.meta$Batch     <- substring(my.grz.meta$SampleID,9)
rownames(my.grz.meta) <- my.grz.meta$SampleID

# reorder count tables in sensical order
for  (i in 1:length(counts.pb.grz)) {
  counts.pb.grz[[i]] <- counts.pb.grz[[i]][,my.grz.meta$SampleID]
}

# Create list object to receive clean SVA counts
sva.cts.grz        <- vector(mode = "list", length = length(counts.pb.grz))
names(sva.cts.grz) <- names(counts.pb.grz)

# Create list object to receive VST normalized counts
vst.cts.grz        <- vector(mode = "list", length = length(counts.pb.grz))
names(vst.cts.grz) <- names(counts.pb.grz)

# Create list object to receive DESeq2 results
deseq.res.list.grz        <- vector(mode = "list", length = length(counts.pb.grz))
names(deseq.res.list.grz) <- names(counts.pb.grz)

# Create list object to receive DESeq2 results (just genes)
deseq.res.list.genes.grz        <- vector(mode = "list", length = length(counts.pb.grz))
names(deseq.res.list.genes.grz) <- names(counts.pb.grz)


# loop over pseudobulk data
for  (i in 1:length(counts.pb.grz)) {
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_Pseudobulk_GRZ_",names(counts.pb.grz)[[i]])
  
  ###################################
  #######       Run SVA      #######
  
  # build design matrix
  sva.dataDesign = data.frame( row.names = my.grz.meta$SampleID , 
                               sex       = my.grz.meta$Sex      ,
                               age       = my.grz.meta$Age_weeks,
                               batch     = my.grz.meta$Batch)
  
  # Set null and alternative models (ignore batch)
  mod1    = model.matrix(~ sex + age + batch, data = sva.dataDesign)
  n.sv.be = num.sv(counts.pb.grz[[i]], mod1, method="be") # microglia is 2
  
  # apply SVAseq algortihm
  my.svseq = svaseq(as.matrix(counts.pb.grz[[i]]), mod1, n.sv=n.sv.be, constant = 0.1)
  
  # remove RIN and SV, preserve age and sex
  my.clean <- removeBatchEffect(log2(counts.pb.grz[[i]] + 0.1), 
                                batch      = my.grz.meta$Batch, 
                                covariates = cbind(my.svseq$sv),
                                design     = mod1[,1:3])
  
  # delog and round data for DEseq2 processing
  my.filtered.sva <- round(2^my.clean-0.1)
  
  # keep only robustly expressed genes
  sva.cts.grz[[i]] <- my.filtered.sva
  
  # legend
  my.cols  <- rep("",nrow(my.grz.meta))
  my.cols[my.grz.meta$Group %in% "GRZ_Y_F"] <- "deeppink"
  my.cols[my.grz.meta$Group %in% "GRZ_M_F"] <- "deeppink3"
  my.cols[my.grz.meta$Group %in% "GRZ_O_F"] <- "deeppink4"
  my.cols[my.grz.meta$Group %in% "GRZ_Y_M"] <- "deepskyblue"
  my.cols[my.grz.meta$Group %in% "GRZ_M_M"] <- "deepskyblue3"
  my.cols[my.grz.meta$Group %in% "GRZ_O_M"] <- "deepskyblue4"
  
  my.pch  <- rep(0,nrow(my.grz.meta))
  my.pch[my.grz.meta$Group %in% "GRZ_Y_F"] <- 16
  my.pch[my.grz.meta$Group %in% "GRZ_M_F"] <- 17
  my.pch[my.grz.meta$Group %in% "GRZ_O_F"] <- 15
  my.pch[my.grz.meta$Group %in% "GRZ_Y_M"] <- 16
  my.pch[my.grz.meta$Group %in% "GRZ_M_M"] <- 17
  my.pch[my.grz.meta$Group %in% "GRZ_O_M"] <- 15
  
  # get matrix using age as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = sva.cts.grz[[i]],
                                colData   = my.grz.meta,
                                design    = ~ Age_weeks + Sex)
  
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
       main= paste0(names(counts.pb.grz)[[i]]," (MDS)"),
       cex=3, col = my.cols, pch = my.pch,
       cex.lab = 1.25,
       cex.axis = 1.25, las = 1)
  dev.off()
  
  # extract gene significance by DEseq2
  res.age <- results(dds.deseq, name = "Age_weeks") # FC per week
  
  # exclude genes with NA FDR value
  res.age <- res.age[!is.na(res.age$padj),]
  
  # restrict analysis to only genes, not TEs
  res.age.genes <- res.age[!grepl("NotFur", rownames(res.age)), ]
  
  # store results
  deseq.res.list.grz[[i]]       <- data.frame(res.age)
  deseq.res.list.genes.grz[[i]] <- data.frame(res.age.genes)
  
  ### get sex dimorphic changes at FDR5
  genes.age <- rownames(res.age.genes)[res.age.genes$padj < 0.05]
  my.num.age <- length(genes.age)
  
  if (my.num.age > 2) {
    # heatmap drawing - only if there is at least 2 gene
    my.heatmap.out <- paste0(my.outprefix,"_AGING_Heatmap_FDR5_GENES.pdf")
    
    pdf(my.heatmap.out, onefile = F, height = 10, width = 10)
    my.heatmap.title <- paste0(names(counts.pb.grz)[[i]], " aging significant (FDR<5%), ", my.num.age, " genes")
    pheatmap::pheatmap(vst.cts.grz[[i]][genes.age,],
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
  
  my.out.stats.age <- paste0(my.outprefix,"_AGING_all_genes_TEs_statistics.txt")
  write.table(deseq.res.list.grz[[i]], file = my.out.stats.age , sep = "\t" , row.names = T, quote = F)
  
  my.out.fdr5.age <- paste0(my.outprefix,"_AGING_FDR5_genes_TEs_statistics.txt")
  write.table(deseq.res.list.grz[[i]][genes.age,], file = my.out.fdr5.age, sep = "\t" , row.names = T, quote = F)
  
  
  #### now do TEs
  te.age <- rownames(res.age)[res.age$padj < 0.05][grep("NotFur1-",rownames(res.age)[res.age$padj < 0.05])]
  
  if (length(te.age) > 2) {
    # heatmap drawing - only if there is at least 2 gene
    my.heatmap.out <- paste0(my.outprefix,"_AGING_Heatmap_FDR5_TEs.pdf")
    
    pdf(my.heatmap.out, onefile = F, height = 10, width = 15)
    my.heatmap.title <- paste0(names(counts.pb.grz)[[i]], " aging significant (FDR<5%), ", length(te.age), " TEs")
    pheatmap::pheatmap(vst.cts.grz[[i]][te.age,],
                       cluster_cols = F,
                       cluster_rows = T,
                       colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                       show_rownames = F, scale="row",
                       main = my.heatmap.title,
                       cellwidth = 15,
                       border = NA,
                       cellheight = 1 )
    dev.off()
  }
  
}

# save R object with all DEseq2 results
my.rdata.age <- paste0(Sys.Date(),"_pseudobulk_killi_cell_types_AGING_GRZ_DEseq2_objects.RData")
save(deseq.res.list.grz, deseq.res.list.genes.grz, file = my.rdata.age)

my.vst.age <- paste0(Sys.Date(),"_pseudobulk_killi_cell_types_AGING_GRZ_VST_data_objects.RData")
save(vst.cts.grz, file = my.vst.age)


######### Make jitter plot of DE genes and TEs #########
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

pdf(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_DESeq2_with_reg_colors_FDR5.pdf"), width = 6, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 8.5),
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
     at = 1:8,
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

png(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_DESeq2_with_reg_colors_FDR5.png"), width = 2400, height = 2200, res = 300, units = "px")
par(mar = c(3.1, 4.1, 2, 2))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 8.5),
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
     at = 1:8,
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


######### Make jitter plot of DE GENES (not TEs) #########
## Order by pvalue:
age.results <- lapply(deseq.res.list.genes.grz,function(x) {x[order(x$padj),]})
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

pdf(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_DESeq2_with_reg_colors_FDR5_GENES_ONLY.pdf"), width = 6, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 8.5),
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
     at = 1:8,
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

png(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_DESeq2_with_reg_colors_FDR5_GENES_ONLY.png"), width = 2400, height = 2200, res = 300, units = "px")
par(mar = c(3.1, 4.1, 2, 2))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 8.5),
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
     at = 1:8,
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
###############################################



###############################################
#######   DEG analysis   ++++   ZMZ     #######
###############################################

# import metadata and order it
my.zmz.meta           <- read.csv("2024-02-19_ZMZ_sample_metadata_table.txt", sep = "\t")
my.zmz.meta$Age_Group <- factor(my.zmz.meta$Age_Group , levels = c("Y", "M","O", "G"))
my.zmz.meta           <- setorder(my.zmz.meta,Age_Group, SampleID)
my.zmz.meta$Batch     <- substring(my.zmz.meta$SampleID,9)
rownames(my.zmz.meta) <- my.zmz.meta$SampleID

# reorder count tables in sensical order
for  (i in 1:length(counts.pb.zmz)) {
  counts.pb.zmz[[i]] <- counts.pb.zmz[[i]][,my.zmz.meta$SampleID]
}

# Create list object to receive clean SVA counts
sva.cts.zmz        <- vector(mode = "list", length = length(counts.pb.zmz))
names(sva.cts.zmz) <- names(counts.pb.zmz)

# Create list object to receive VST normalized counts
vst.cts.zmz        <- vector(mode = "list", length = length(counts.pb.zmz))
names(vst.cts.zmz) <- names(counts.pb.zmz)

# Create list object to receive DESeq2 results
deseq.res.list.zmz        <- vector(mode = "list", length = length(counts.pb.zmz))
names(deseq.res.list.zmz) <- names(counts.pb.zmz)

# Create list object to receive DESeq2 results (just genes)
deseq.res.list.genes.zmz        <- vector(mode = "list", length = length(counts.pb.zmz))
names(deseq.res.list.genes.zmz) <- names(counts.pb.zmz)


# loop over pseudobulk data
for  (i in 1:length(counts.pb.zmz)) {
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_Pseudobulk_ZMZ_",names(counts.pb.zmz)[[i]])
  
  ###################################
  #######       Run SVA      #######
  
  # build design matrix
  sva.dataDesign = data.frame( row.names = my.zmz.meta$SampleID , 
                               sex       = my.zmz.meta$Sex      ,
                               age       = my.zmz.meta$Age_weeks,
                               batch     = my.zmz.meta$Batch)
  
  # Set null and alternative models (ignore batch)
  mod1    = model.matrix(~ sex + age + batch, data = sva.dataDesign)
  n.sv.be = num.sv(counts.pb.zmz[[i]], mod1, method="be") # microglia is 1
  
  # apply SVAseq algortihm
  my.svseq = svaseq(as.matrix(counts.pb.zmz[[i]]), mod1, n.sv=n.sv.be, constant = 0.1)
  
  # remove RIN and SV, preserve age and sex
  my.clean <- removeBatchEffect(log2(counts.pb.zmz[[i]] + 0.1), 
                                batch      = my.zmz.meta$Batch, 
                                covariates = cbind(my.svseq$sv),
                                design     = mod1[,1:3])
  
  # delog and round data for DEseq2 processing
  my.filtered.sva <- round(2^my.clean-0.1)
  
  # keep only robustly expressed genes
  sva.cts.zmz[[i]] <- my.filtered.sva
  
  # legend
  my.cols  <- rep("",nrow(my.zmz.meta))
  my.cols[my.zmz.meta$Group %in% "ZMZ_Y_F"] <- "deeppink"
  my.cols[my.zmz.meta$Group %in% "ZMZ_M_F"] <- "deeppink3"
  my.cols[my.zmz.meta$Group %in% "ZMZ_O_F"] <- "deeppink4"
  my.cols[my.zmz.meta$Group %in% "ZMZ_G_F"] <- "magenta4"
  my.cols[my.zmz.meta$Group %in% "ZMZ_Y_M"] <- "deepskyblue"
  my.cols[my.zmz.meta$Group %in% "ZMZ_M_M"] <- "deepskyblue3"
  my.cols[my.zmz.meta$Group %in% "ZMZ_O_M"] <- "deepskyblue4"
  my.cols[my.zmz.meta$Group %in% "ZMZ_G_M"] <- "royalblue4"
  
  my.pch  <- rep(0,nrow(my.zmz.meta))
  my.pch[my.zmz.meta$Group %in% "ZMZ_Y_F"] <- 16
  my.pch[my.zmz.meta$Group %in% "ZMZ_M_F"] <- 17
  my.pch[my.zmz.meta$Group %in% "ZMZ_O_F"] <- 15
  my.pch[my.zmz.meta$Group %in% "ZMZ_G_F"] <- 18
  my.pch[my.zmz.meta$Group %in% "ZMZ_Y_M"] <- 16
  my.pch[my.zmz.meta$Group %in% "ZMZ_M_M"] <- 17
  my.pch[my.zmz.meta$Group %in% "ZMZ_O_M"] <- 15
  my.pch[my.zmz.meta$Group %in% "ZMZ_G_M"] <- 18
  
  # get matrix using age as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = sva.cts.zmz[[i]],
                                colData   = my.zmz.meta,
                                design    = ~ Age_weeks + Sex)
  
  # run DESeq normalizations and export results
  dds.deseq <- DESeq(dds)
  
  # plot dispersion
  my.disp.out <- paste(my.outprefix,"_dispersion_plot.pdf")
  
  pdf(my.disp.out)
  plotDispEsts(dds.deseq)
  dev.off()
  
  # get DESeq2 normalized expression value
  vst.cts.zmz[[i]] <- getVarianceStabilizedData(dds.deseq)
  
  # MDS analysis
  mds.result <- cmdscale(1-cor(vst.cts.zmz[[i]],method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  x <- mds.result[, 1]
  y <- mds.result[, 2]
  
  pdf(paste0(my.outprefix,"_MDS_plot.pdf"))
  plot(x, y,
       xlab = "MDS dimension 1", ylab = "MDS dimension 2",
       main= paste0(names(counts.pb.zmz)[[i]]," (MDS)"),
       cex=3, col = my.cols, pch = my.pch,
       cex.lab = 1.25,
       cex.axis = 1.25, las = 1)
  dev.off()
  
  # extract gene significance by DEseq2
  res.age <- results(dds.deseq, name = "Age_weeks") # FC per week
  
  # exclude genes with NA FDR value
  res.age <- res.age[!is.na(res.age$padj),]
  
  # restrict analysis to only genes, not TEs
  res.age.genes <- res.age[!grepl("NotFur", rownames(res.age)), ]
  
  # store results
  deseq.res.list.zmz[[i]]       <- data.frame(res.age)
  deseq.res.list.genes.zmz[[i]] <- data.frame(res.age.genes)
  
  ### get sex dimorphic changes at FDR5
  genes.age <- rownames(res.age.genes)[res.age.genes$padj < 0.05]
  my.num.age <- length(genes.age)
  
  if (my.num.age > 2) {
    # heatmap drawing - only if there is at least 2 gene
    my.heatmap.out <- paste0(my.outprefix,"_AGING_Heatmap_FDR5_GENES.pdf")
    
    pdf(my.heatmap.out, onefile = F, height = 10, width = 10)
    my.heatmap.title <- paste0(names(counts.pb.zmz)[[i]], " aging significant (FDR<5%), ", my.num.age, " genes")
    pheatmap::pheatmap(vst.cts.zmz[[i]][genes.age,],
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
  write.table(vst.cts.zmz[[i]], file = my.out.ct.mat , sep = "\t" , row.names = T, quote = F)
  
  my.out.stats.age <- paste0(my.outprefix,"_AGING_all_genes_TEs_statistics.txt")
  write.table(deseq.res.list.zmz[[i]], file = my.out.stats.age , sep = "\t" , row.names = T, quote = F)
  
  my.out.fdr5.age <- paste0(my.outprefix,"_AGING_FDR5_genes_TEs_statistics.txt")
  write.table(deseq.res.list.zmz[[i]][genes.age,], file = my.out.fdr5.age, sep = "\t" , row.names = T, quote = F)
  
  
  #### now do TEs
  te.age <- rownames(res.age)[res.age$padj < 0.05][grep("NotFur1-",rownames(res.age)[res.age$padj < 0.05])]
  
  if (length(te.age) > 2) {
    # heatmap drawing - only if there is at least 2 gene
    my.heatmap.out <- paste0(my.outprefix,"_AGING_Heatmap_FDR5_TEs.pdf")
    
    pdf(my.heatmap.out, onefile = F, height = 10, width = 15)
    my.heatmap.title <- paste0(names(counts.pb.zmz)[[i]], " aging significant (FDR<5%), ", length(te.age), " TEs")
    pheatmap::pheatmap(vst.cts.zmz[[i]][te.age,],
                       cluster_cols = F,
                       cluster_rows = T,
                       colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                       show_rownames = F, scale="row",
                       main = my.heatmap.title,
                       cellwidth = 15,
                       border = NA,
                       cellheight = 1 )
    dev.off()
  }
  
}

# save R object with all DEseq2 results
my.rdata.age <- paste0(Sys.Date(),"_pseudobulk_killi_cell_types_AGING_ZMZ_DEseq2_objects.RData")
save(deseq.res.list.zmz, deseq.res.list.genes.zmz, file = my.rdata.age)

my.vst.age <- paste0(Sys.Date(),"_pseudobulk_killi_cell_types_AGING_ZMZ_VST_data_objects.RData")
save(vst.cts.zmz, file = my.vst.age)


######### Make jitter plot of DE genes and TEs #########
## Order by pvalue:
age.results <- lapply(deseq.res.list.zmz,function(x) {x[order(x$padj),]})
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

pdf(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_ZMZ_stripplot_DESeq2_with_reg_colors_FDR5.pdf"), width = 6, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 8.5),
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
     at = 1:8,
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

png(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_ZMZ_stripplot_DESeq2_with_reg_colors_FDR5.png"), width = 2400, height = 2200, res = 300, units = "px")
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 8.5),
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
     at = 1:8,
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


######### Make jitter plot of DE GENES (not TEs) #########
## Order by pvalue:
age.results <- lapply(deseq.res.list.genes.zmz,function(x) {x[order(x$padj),]})
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

pdf(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_ZMZ_stripplot_DESeq2_with_reg_colors_FDR5_GENES_ONLY.pdf"), width = 6, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 8.5),
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
     at = 1:8,
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

png(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_ZMZ_stripplot_DESeq2_with_reg_colors_FDR5_GENES_ONLY.png"), width = 2400, height = 2200, res = 300, units = "px")
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 8.5),
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
     at = 1:8,
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
###############################################################################################


###############################################################################################
# 3. cross cell type comparison of genes with aging gene expression

# Clean workspace and reload necessary objects
load('2024-02-28_pseudobulk_killi_cell_types_AGING_GRZ_DEseq2_objects.RData'   )
load('2024-02-28_pseudobulk_killi_cell_types_AGING_GRZ_VST_data_objects.RData' )
load('2024-02-28_pseudobulk_killi_cell_types_AGING_ZMZ_DEseq2_objects.RData'   )
load('2024-02-28_pseudobulk_killi_cell_types_AGING_ZMZ_VST_data_objects.RData' )

######################################################
#######  A. Correlation analysis   ++++   GRZ  #######
######################################################

# get the common expressed background for comparison
my.exp.genes.grz <- Reduce(intersect,lapply(vst.cts.grz,rownames))

###
my.logFC.data.grz <- data.frame(matrix(NA,length(my.exp.genes.grz),length(deseq.res.list.genes.grz)))
colnames(my.logFC.data.grz) <- names(deseq.res.list.genes.grz)
rownames(my.logFC.data.grz) <- my.exp.genes.grz

my.FDR.data.grz <- data.frame(matrix(NA,length(my.exp.genes.grz),length(deseq.res.list.genes.grz)))
colnames(my.FDR.data.grz) <- names(deseq.res.list.genes.grz)
rownames(my.FDR.data.grz) <- my.exp.genes.grz

for (i in 1:length(deseq.res.list.genes.grz)) {
  for (j in 1:length(my.exp.genes.grz)) {
    my.idx <- rownames(deseq.res.list.genes.grz[[i]]) %in% my.exp.genes.grz[j]
    
    if (sum(my.idx) > 0) {
      
      my.logFC.data.grz[j,i] <- deseq.res.list.genes.grz[[i]]$log2FoldChange[my.idx]
      my.FDR.data.grz[j,i]   <- deseq.res.list.genes.grz[[i]]$padj[my.idx]
      
    }
  }
}

# get spearman rank correlation
my.cors.grz <- cor(my.logFC.data.grz, method = 'spearman', use = "complete.obs")

pdf(paste0(Sys.Date(),"_Killi_Brain_Aging_GRZ_correlation_of_FC_heatmap.pdf"))
Heatmap(my.cors.grz, border = T, rect_gp = gpar(col = "grey", lwd = 0.5),
        column_title = "log2FC correlation across cell types with Aging")
dev.off()


#########
source('Plot_Jaccard_function_AGING.R')
##
jacc.5 <- calc_jaccard(my.logFC.data.grz, my.FDR.data.grz, 0.05)

pdf(paste0(Sys.Date(),"_GRZ_Jaccard_Index_FDR5.pdf"))
Heatmap(jacc.5[[1]],
        col = colorRamp2(c(0,1), c("white","#CC3333"), transparency = 0, space = "LAB"),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5)  ,
        column_title = "Jaccard Aging Up")
Heatmap(jacc.5[[2]],
        col = colorRamp2(c(0,1), c("white","#333399"), transparency = 0, space = "LAB"),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
        column_title = "Jaccard Aging Down")
dev.off()


# genes reccurently up or down with aging
my.up.grz  <- my.logFC.data.grz > 0
my.sig.grz <- my.FDR.data.grz   < 0.05

up.genes.grz   <- data.frame(matrix(0,length(my.exp.genes.grz),length(deseq.res.list.genes.grz)))
colnames(up.genes.grz) <- names(deseq.res.list.genes.grz)
rownames(up.genes.grz) <- my.exp.genes.grz

down.genes.grz <- data.frame(matrix(0,length(my.exp.genes.grz),length(deseq.res.list.genes.grz)))
colnames(down.genes.grz) <- names(deseq.res.list.genes.grz)
rownames(down.genes.grz) <- my.exp.genes.grz

for (i in 1:length(deseq.res.list.genes.grz)) {
  for (j in 1:length(my.exp.genes.grz)) {
    
    up.genes.grz[j,i]   <- bitAnd(my.up.grz[j,i], my.sig.grz[j,i])
    down.genes.grz[j,i] <- bitAnd(!my.up.grz[j,i], my.sig.grz[j,i])
    
  }
}
up.genes.grz[is.na(up.genes.grz)]     <- 0
down.genes.grz[is.na(down.genes.grz)] <- 0

####
rownames(up.genes.grz)[apply(up.genes.grz,1,sum)>4] ### 16
#  [1] "LOC107374876" "LOC107374900" "LOC107374381" "LOC107376528" "helz2"        "irf9"         "stat1"        "znfx1"        "rnf213"       "LOC107387991"
#  [11] "LOC107389223" "parp9"        "ifih1"        "LOC107390931" "LOC107392596" "LOC107392842"

rownames(down.genes.grz)[apply(down.genes.grz,1,sum)>4]
# character(0)


######################################################
#######  B. Correlation analysis   ++++   ZMZ  #######
######################################################

# get the common expressed background for comparison
my.exp.genes.zmz <- Reduce(intersect,lapply(vst.cts.zmz,rownames))

###
my.logFC.data.zmz <- data.frame(matrix(NA,length(my.exp.genes.zmz),length(deseq.res.list.genes.zmz)))
colnames(my.logFC.data.zmz) <- names(deseq.res.list.genes.zmz)
rownames(my.logFC.data.zmz) <- my.exp.genes.zmz

my.FDR.data.zmz <- data.frame(matrix(NA,length(my.exp.genes.zmz),length(deseq.res.list.genes.zmz)))
colnames(my.FDR.data.zmz) <- names(deseq.res.list.genes.zmz)
rownames(my.FDR.data.zmz) <- my.exp.genes.zmz

for (i in 1:length(deseq.res.list.genes.zmz)) {
  for (j in 1:length(my.exp.genes.zmz)) {
    my.idx <- rownames(deseq.res.list.genes.zmz[[i]]) %in% my.exp.genes.zmz[j]
    
    if (sum(my.idx) > 0) {
      
      my.logFC.data.zmz[j,i] <- deseq.res.list.genes.zmz[[i]]$log2FoldChange[my.idx]
      my.FDR.data.zmz[j,i]   <- deseq.res.list.genes.zmz[[i]]$padj[my.idx]
      
    }
  }
}

# get spearman rank correlation
my.cors.zmz <- cor(my.logFC.data.zmz, method = 'spearman', use = "complete.obs")

pdf(paste0(Sys.Date(),"_Killi_Brain_Aging_ZMZ_correlation_of_FC_heatmap.pdf"))
Heatmap(my.cors.zmz, border = T, rect_gp = gpar(col = "grey", lwd = 0.5),
        column_title = "log2FC correlation across cell types with Aging")
dev.off()


#########
source('Plot_Jaccard_function_AGING.R')
##
jacc.5 <- calc_jaccard(my.logFC.data.zmz, my.FDR.data.zmz, 0.05)

pdf(paste0(Sys.Date(),"_ZMZ_Jaccard_Index_FDR5.pdf"))
Heatmap(jacc.5[[1]],
        col = colorRamp2(c(0,1), c("white","#CC3333"), transparency = 0, space = "LAB"),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5)  ,
        column_title = "Jaccard Aging Up")
Heatmap(jacc.5[[2]],
        col = colorRamp2(c(0,1), c("white","#333399"), transparency = 0, space = "LAB"),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
        column_title = "Jaccard Aging Down")
dev.off()


# genes reccurently up or down with aging
my.up.zmz  <- my.logFC.data.zmz > 0
my.sig.zmz <- my.FDR.data.zmz   < 0.05

up.genes.zmz   <- data.frame(matrix(0,length(my.exp.genes.zmz),length(deseq.res.list.genes.zmz)))
colnames(up.genes.zmz) <- names(deseq.res.list.genes.zmz)
rownames(up.genes.zmz) <- my.exp.genes.zmz

down.genes.zmz <- data.frame(matrix(0,length(my.exp.genes.zmz),length(deseq.res.list.genes.zmz)))
colnames(down.genes.zmz) <- names(deseq.res.list.genes.zmz)
rownames(down.genes.zmz) <- my.exp.genes.zmz

for (i in 1:length(deseq.res.list.genes.zmz)) {
  for (j in 1:length(my.exp.genes.zmz)) {
    
    up.genes.zmz[j,i]   <- bitAnd(my.up.zmz[j,i], my.sig.zmz[j,i])
    down.genes.zmz[j,i] <- bitAnd(!my.up.zmz[j,i], my.sig.zmz[j,i])
    
  }
}
up.genes.zmz[is.na(up.genes.zmz)]     <- 0
down.genes.zmz[is.na(down.genes.zmz)] <- 0

####
rownames(up.genes.zmz)[apply(up.genes.zmz,1,sum)>4] ### 32
# [1] "irf2"         "LOC107375315" "LOC107374900" "syn3"         "LOC107386289" "LOC107374381" "LOC107376528" "helz2"        "LOC107376321" "LOC107376289"
# [11] "LOC107377161" "ezh1"         "dhx58"        "LOC107382742" "irf9"         "stat1"        "znfx1"        "rnf213"       "LOC107387991" "LOC107389223"
# [21] "ifih1"        "LOC107390133" "LOC107390931" "LOC107391063" "LOC107392516" "LOC107392912" "LOC107392842" "LOC107396709" "LOC107372890" "LOC107373776"
# [31] "LOC107373896" "LOC107374505"

rownames(down.genes.zmz)[apply(down.genes.zmz,1,sum)>4]
# [1] "LOC107386615"


######################################################
#######   C. Overlap analysis      GRZ/ ZMZ    #######
######################################################

########## Compute overlap for each cell types
fisher.res <- data.frame(matrix(0,length(deseq.res.list.genes.grz), 5))
colnames(fisher.res) <- c("Cell_Type","Overlap_up","Overlap_up_enrich_Fisher","Overlap_dwn","Overlap_dwn_enrich_Fisher")
fisher.res$Cell_Type <- names(deseq.res.list.genes.grz)

celltype.common.data <- data.frame(col.names = c("Gene", "Cell_Type",
                                                 "baseMean.grz", "log2FoldChange.grz", "padj.grz",
                                                 "baseMean.zmz", "log2FoldChange.zmz", "padj.zmz"))

for (i in 1:length(deseq.res.list.genes.grz)) {
  
  if (i != 2) {  # take care of case where GRZ does not have sig genes (GABAergic neurons)
    
    # gene list
    grz.up   <- rownames(deseq.res.list.genes.grz[[i]])[bitAnd(deseq.res.list.genes.grz[[i]]$padj < 0.05, 
                                                               deseq.res.list.genes.grz[[i]]$log2FoldChange >0)>0]
    grz.dwn  <- rownames(deseq.res.list.genes.grz[[i]])[bitAnd(deseq.res.list.genes.grz[[i]]$padj < 0.05, 
                                                               deseq.res.list.genes.grz[[i]]$log2FoldChange <0)>0]
    zmz.up   <-rownames(deseq.res.list.genes.zmz[[i]])[bitAnd(deseq.res.list.genes.zmz[[i]]$padj < 0.05, 
                                                              deseq.res.list.genes.zmz[[i]]$log2FoldChange >0)>0]
    zmz.dwn  <- rownames(deseq.res.list.genes.zmz[[i]])[bitAnd(deseq.res.list.genes.zmz[[i]]$padj < 0.05, 
                                                               deseq.res.list.genes.zmz[[i]]$log2FoldChange <0)>0]
    genes.bckd  <- union(rownames(deseq.res.list.genes.grz[[i]]),rownames(deseq.res.list.genes.zmz[[i]]))
    
    ################ Upregulated genes
    my.criteria.up <- list("GRZ Up"    = grz.up,  
                           "ZMZ Up"    = zmz.up)  
    my.Venn.up <- Venn(my.criteria.up)
    
    
    pdf(paste0(Sys.Date(),"_VennDiagram_",fisher.res$Cell_Type[i],"_UpAging_Overlap_ZMZ_GRZ.pdf"))
    plot(my.Venn.up, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
    dev.off()
    
    # prepare Fisher hypergeometric test
    a <- length(my.Venn.up@IntersectionSets$`11`)
    b <- length(my.Venn.up@IntersectionSets$`10`)
    c <- length(my.Venn.up@IntersectionSets$`01`)
    d <- length(genes.bckd) - a - b - c
    
    my.fisher.up <- fisher.test(matrix(c(a,b,c,d),2,2))
    
    # put in result table
    fisher.res$Overlap_up[i]                <- a
    fisher.res$Overlap_up_enrich_Fisher[i]  <- my.fisher.up$p.value
    
    
    ################ Downregulated genes
    my.criteria.dwn <- list("GRZ down"    = grz.dwn,  
                            "ZMZ down"    = zmz.dwn)  
    my.Venn.dwn <- Venn(my.criteria.dwn)
    
    
    pdf(paste0(Sys.Date(),"_VennDiagram_",fisher.res$Cell_Type[i],"_DownAging_Overlap_ZMZ_GRZ.pdf"))
    plot(my.Venn.dwn, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
    dev.off()
    
    # prepare Fisher hypergeometric test
    a <- length(my.Venn.dwn@IntersectionSets$`11`)
    b <- length(my.Venn.dwn@IntersectionSets$`10`)
    c <- length(my.Venn.dwn@IntersectionSets$`01`)
    d <- length(genes.bckd) - a - b - c
    
    my.fisher.dwn <- fisher.test(matrix(c(a,b,c,d),2,2))
    
    # put in result table
    fisher.res$Overlap_dwn[i]                <- a
    fisher.res$Overlap_dwn_enrich_Fisher[i]  <- my.fisher.dwn$p.value
    
    
    ##########################################
    # extract stats/information for overlapping genes
    data.up.grz   <- deseq.res.list.genes.grz[[i]][my.Venn.up@IntersectionSets$`11`,]
    data.up.zmz   <- deseq.res.list.genes.zmz[[i]][my.Venn.up@IntersectionSets$`11`,]
    data.up.merge <- merge(data.up.grz, data.up.zmz, by = "row.names", suffixes = c(".grz",".zmz"))
    colnames(data.up.merge)[1] <- "Gene"
    
    data.dwn.grz   <- deseq.res.list.genes.grz[[i]][my.Venn.dwn@IntersectionSets$`11`,]
    data.dwn.zmz   <- deseq.res.list.genes.zmz[[i]][my.Venn.dwn@IntersectionSets$`11`,]
    data.dwn.merge <- merge(data.dwn.grz, data.dwn.zmz, by = "row.names", suffixes = c(".grz",".zmz"))
    colnames(data.dwn.merge)[1] <- "Gene"
    
    
    my.tmp <- rbind(data.up.merge, data.dwn.merge)
    my.tmp$Cell_Type <- fisher.res$Cell_Type[i]
    my.tmp <- my.tmp[,c("Gene","Cell_Type", "baseMean.grz", "log2FoldChange.grz", "padj.grz", "baseMean.zmz", "log2FoldChange.zmz", "padj.zmz")]
    
    if ( i == 1) {
      celltype.common.data <- my.tmp
    } else {
      celltype.common.data <- rbind(celltype.common.data, my.tmp)
    }
    
  }
}

# export significance results
# up overlap is always more than expected by chance, down is more variable
write.table(fisher.res, file = paste0(Sys.Date(),"_Fisher_Aging_FDR5_byCellType_Overlap_ZMZ_GRZ.txt"), sep = "\t", row.names = F, quote = F)

# add combined fisher FDR for help sorting
celltype.common.data$FisherCombinedPval <- celltype.common.data$padj.grz * celltype.common.data$padj.zmz
write.table(celltype.common.data, file = paste0(Sys.Date(),"_Aging_DEseq2_results_FDR5_byCellType_Overlap_ZMZ_GRZ.txt"), sep = "\t", row.names = F, quote = F)




########## Compute overlap for recurrent UPregulated
# no need for down, since GRZ doesn't have any
up.rec.grz  <- rownames(up.genes.grz)[apply(up.genes.grz,1,sum)>4]
up.rec.zmz  <- rownames(up.genes.zmz)[apply(up.genes.zmz,1,sum)>4]

my.criteria.rec <- list("GRZ recurrent Up"    = up.rec.grz,  
                        "ZMZ recurrent Up"    = up.rec.zmz)  
my.Venn.rec <- Venn(my.criteria.rec)

pdf(paste0(Sys.Date(),"_VennDiagram_Recurrent_UpAging_Overlap_ZMZ_GRZ.pdf"))
plot(my.Venn.rec, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
dev.off()

# all detected genes for background
bckd.all <- union(unlist(lapply(deseq.res.list.genes.grz,rownames)),unlist(lapply(deseq.res.list.genes.zmz,rownames)))

my.fisher.all <- fisher.test(matrix(c(13, 3, 19, length(bckd.all) - 13 - 3 - 19),2,2))
my.fisher.all$p.value # 3.798167e-35

my.Venn.rec@IntersectionSets$`11`
# [1] "LOC107374900" "LOC107374381" "LOC107376528" "helz2"        "irf9"         "stat1"        "znfx1"        "rnf213"       "LOC107387991" "LOC107389223"
# [11] "ifih1"        "LOC107390931" "LOC107392842"

write.table(my.Venn.rec@IntersectionSets$`11`, file = paste0(Sys.Date(),"_Recurrent_genes_UP_FDR5_overlap_ZMZ_GRZ.txt"), sep = "\t", 
            row.names = F, quote = F, col.names = F)


##########  Plot heatmaps for recurrent genes
my.rec.genes <- read.table('2024-02-28_Recurrent_genes_UP_FDR5_overlap_ZMZ_GRZ.txt')

my.grz.data <- matrix(0, length(my.rec.genes$V1), length(deseq.res.list.genes.grz))
colnames(my.grz.data) <- names(deseq.res.list.genes.grz)
rownames(my.grz.data) <- my.rec.genes$V1

# get log2FC for all the rec genes
for (i in 1:length(deseq.res.list.genes.grz)) {
  my.grz.data[,i] <- deseq.res.list.genes.grz[[i]][my.rec.genes$V1,]$log2FoldChange
}

my.zmz.data <- matrix(0, length(my.rec.genes$V1), length(deseq.res.list.genes.zmz))
colnames(my.zmz.data) <- names(deseq.res.list.genes.zmz)
rownames(my.zmz.data) <- my.rec.genes$V1

# get log2FC for all the rec genes
for (i in 1:length(deseq.res.list.genes.zmz)) {
  my.zmz.data[,i] <- deseq.res.list.genes.zmz[[i]][my.rec.genes$V1,]$log2FoldChange
}


my.consensus.sort.order <- sort(apply(cbind(my.grz.data,my.zmz.data),1,median), index.return = T, decreasing = T)

pdf(paste0(Sys.Date(),"_Recurrent_genes_UP_FDR5_overlap_ZMZ_GRZ_HEATMAP_in_GRZ.pdf"))
grz.heat <- Heatmap(my.grz.data[my.consensus.sort.order$ix,], 
                    col = colorRamp2(c(0,0.7), c("white", "#CC3333"), transparency = 0, space = "LAB"),
                    border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
                    cluster_rows = F,
                    cluster_columns = F,
                    column_title = "GRZ log2(FC)")
grz.heat
dev.off()

pdf(paste0(Sys.Date(),"_Recurrent_genes_UP_FDR5_overlap_ZMZ_GRZ_HEATMAP-in_ZMZ.pdf"))
zmz.heat <- Heatmap(my.zmz.data[my.consensus.sort.order$ix,], 
                    col = colorRamp2(c(0,0.7), c("white", "#CC3333"), transparency = 0, space = "LAB"),
                    border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
                    cluster_rows = F,
                    cluster_columns = F,
                    column_title = "ZMZ log2(FC)")
zmz.heat
dev.off()

pdf(paste0(Sys.Date(),"_Recurrent_genes_UP_FDR5_overlap_ZMZ_GRZ_HEATMAP_in_GRZ_and_ZMZ_PANELED.pdf"), width = 6, height = 5)
grz.heat + zmz.heat
dev.off()


##########  Plot Violins for recurrent genes
my.rec.genes <- read.table('2024-02-28_Recurrent_genes_UP_FDR5_overlap_ZMZ_GRZ.txt')
qc.celltypes <- names(deseq.res.list.genes.grz)

load('2024-02-16_Seurat_objects_SPLIT_PER_STRAIN_withSampleID_with_Manual_annotation_FINAL.RData')
killi.brain.grz <- subset(killi.brain.grz, subset = Cell_Identity %in% qc.celltypes)
killi.brain.zmz <- subset(killi.brain.zmz, subset = Cell_Identity %in% qc.celltypes)

killi.brain.grz$Group <- factor(killi.brain.grz$Group, levels = c("GRZ_Y_F","GRZ_M_F","GRZ_O_F","GRZ_Y_M","GRZ_M_M","GRZ_O_M"))
killi.brain.zmz$Group <- factor(killi.brain.zmz$Group, levels = c("ZMZ_Y_F","ZMZ_M_F","ZMZ_O_F","ZMZ_G_F","ZMZ_Y_M","ZMZ_M_M","ZMZ_O_M","ZMZ_G_M"))

killi.brain.grz$Age_Group <- factor(killi.brain.grz$Age_Group, levels = c("Y","M","O"))
killi.brain.zmz$Age_Group <- factor(killi.brain.zmz$Age_Group, levels = c("Y","M","O", "G"))

killi.brain.grz$Age_weeks <- as.numeric(killi.brain.grz$Age_weeks)
killi.brain.zmz$Age_weeks <- as.numeric(killi.brain.zmz$Age_weeks)

pdf(paste0(Sys.Date(),"_Seurat_violinPlot_expression_recurrent_Up_GRZ_byAge.pdf"), height = 8, width = 10)
VlnPlot(killi.brain.grz, 
        features = rev(my.rec.genes$V1), 
        split.by = "Age_Group", 
        cols = c("grey90" ,
                 "grey65",
                 "grey45"),
        pt.size = 0, assay = 'RNA', stack = T, flip = T)
dev.off()  

pdf(paste0(Sys.Date(),"_Seurat_violinPlot_expression_recurrent_Up_ZMZ_byAge.pdf"), height = 8, width = 10)
VlnPlot(killi.brain.zmz, 
        features = rev(my.rec.genes$V1), 
        split.by = "Age_Group", 
        cols = c("grey90" ,
                 "grey65",
                 "grey45",
                 "grey30"),
        pt.size = 0, assay = 'RNA', stack = T, flip = T)
dev.off()  


pdf(paste0(Sys.Date(),"_Seurat_violinPlot_expression_helz2_perBioGroup_GRZ.pdf"), height = 4, width = 8)
VlnPlot(killi.brain.grz, 
        features = "helz2", 
        split.by = "Group", 
        cols = c("deeppink" ,"deeppink3","deeppink4",
                 "deepskyblue","deepskyblue3","deepskyblue4"),
        pt.size = 0, assay = 'RNA', flip = T)
dev.off()  

pdf(paste0(Sys.Date(),"_Seurat_violinPlot_expression_helz2_perBioGroup_ZMZ.pdf"), height = 4, width = 9)
VlnPlot(killi.brain.zmz, 
        features = "helz2", 
        split.by = "Group", 
        cols = c("deeppink" ,"deeppink3","deeppink4","magenta4",
                 "deepskyblue","deepskyblue3","deepskyblue4","royalblue4"),
        pt.size = 0, assay = 'RNA', flip = T)
dev.off()  

# rnf213a / LOC107374381


### get TE frequency as metadata
killi.brain.grz[["percent.TE"]]   <- PercentageFeatureSet(object = killi.brain.grz, pattern = "^NotFur")
killi.brain.zmz[["percent.TE"]]   <- PercentageFeatureSet(object = killi.brain.zmz, pattern = "^NotFur")
killi.brain.grz[["percent.LINE"]] <- PercentageFeatureSet(object = killi.brain.grz, pattern = "LINE")
killi.brain.zmz[["percent.LINE"]] <- PercentageFeatureSet(object = killi.brain.zmz, pattern = "LINE")

pdf(paste0(Sys.Date(),"_Boxplot_TE_LINE_UMI_percentage.pdf"), width = 8, height = 8)
par(mfrow = c(2,2))
boxplot(percent.TE ~ Group, 
        data = killi.brain.grz@meta.data,
        las = 2, 
        col = c("deeppink" ,"deeppink3","deeppink4",
                 "deepskyblue","deepskyblue3","deepskyblue4"),
        ylab = "% TE UMIs",
        xlab = 0,
        outline = F,
        ylim = c(10,50)
)

boxplot(percent.TE ~ Group, 
        data = killi.brain.zmz@meta.data,
        las = 2, 
        col = c("deeppink" ,"deeppink3","deeppink4","magenta4",
                "deepskyblue","deepskyblue3","deepskyblue4","royalblue4"),
        ylab = "% TE UMIs",
        xlab = 0,
        outline = F,
        ylim = c(10,50)
)

boxplot(percent.LINE ~ Group, 
        data = killi.brain.grz@meta.data,
        las = 2, 
        col = c("deeppink" ,"deeppink3","deeppink4",
                "deepskyblue","deepskyblue3","deepskyblue4"),
        ylab = "% LINE UMIs",
        xlab = 0,
        outline = F,
        ylim = c(10,30)
)

boxplot(percent.LINE ~ Group, 
        data = killi.brain.zmz@meta.data,
        las = 2, 
        col = c("deeppink" ,"deeppink3","deeppink4","magenta4",
                "deepskyblue","deepskyblue3","deepskyblue4","royalblue4"),
        ylab = "% LINE UMIs",
        xlab = 0,
        outline = F,
        ylim = c(10,30)
)
dev.off()


te.grz.lm <- lm(percent.TE ~ Age_weeks + Sex, data = killi.brain.grz@meta.data)
summary(te.grz.lm)
#                Estimate Std. Error t value Pr(>|t|)    
#   (Intercept) 26.190728   0.043506 602.000  < 2e-16 ***
#   Age_weeks    0.092160   0.003604  25.573  < 2e-16 ***
#   SexM        -0.102112   0.029808  -3.426 0.000614 ***
  
te.zmz.lm <- lm(percent.TE ~ Age_weeks + Sex, data = killi.brain.zmz@meta.data)
summary(te.zmz.lm)
#               Estimate Std. Error t value Pr(>|t|)    
#   (Intercept) 27.164912   0.034737  782.02   <2e-16 ***
#   Age_weeks    0.024858   0.001994   12.46   <2e-16 ***
#   SexM         0.743903   0.028647   25.97   <2e-16 ***  

line.grz.lm <- lm(percent.LINE ~ Age_weeks + Sex, data = killi.brain.grz@meta.data)
summary(line.grz.lm)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 17.219289   0.030612  562.51   <2e-16 ***
# Age_weeks    0.058092   0.002536   22.91   <2e-16 ***
# SexM        -0.260794   0.020974  -12.43   <2e-16 ***

line.zmz.lm <- lm(percent.LINE ~ Age_weeks + Sex, data = killi.brain.zmz@meta.data)
summary(line.zmz.lm)
#               Estimate Std. Error t value Pr(>|t|)    
#  (Intercept) 17.855820   0.024458  730.06   <2e-16 ***
#  Age_weeks    0.015328   0.001404   10.92   <2e-16 ***
#  SexM         0.447401   0.020170   22.18   <2e-16 ***

##########################################################################################################################################


#######################
sink(file = paste(Sys.Date(),"_MuscatDEseq2_PB_DESeq2_GSEA_scRNAseq_BrainAtlas_session_Info.txt", sep =""))
sessionInfo()
sink()