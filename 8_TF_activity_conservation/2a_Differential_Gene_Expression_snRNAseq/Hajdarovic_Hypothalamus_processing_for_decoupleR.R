setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Mouse_Datasets_for_Comparison/Muscat_DESeq2')
options(stringsAsFactors = F)

#### Packages
library('Seurat')         # 
library(sctransform)      # 
library("singleCellTK")   # 
library("anndata")

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

# 2024-03-22
# Parse/annotate Hajdarovic Brain data for DecoupleR

################################################################################################################################################################
#### 0. read data and metadata, create new Seurat and SCE objects
hyp.data <- readRDS('./Processed_Objects/Hajdarovic_GSE188646_hypo.integrated.final.20210719.RDS')
hyp.data
# An object of class Seurat 
# 27135 features across 40064 samples within 2 assays 
# Active assay: integrated (2000 features, 2000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

table(hyp.data@meta.data$orig.ident,hyp.data@meta.data$group)
#         Oligodendrocyte Neuron Astrocyte  OPC Microglia/Macrophage Tanycyte Ependymocyte Pericyte/Endothelial Cell Vascular and Leptomeningeal Cell
# Aged_1             1246   4224       720  247                  353      133           16                        41                               31
# Aged_2             1030   5248       926  220                  288      142           32                        35                               36
# Aged_3              455   3091       400   93                   98       55           45                        56                               19
# Aged_4              606   3259       358  108                   96       64           10                        20                                7
# Young_1              87    567       373   22                   42       84           38                         8                                8
# Young_2             838   4255       898  208                  192      112           18                        30                               22
# Young_3             316   2824       351   90                   54       48           29                        18                               10
# Young_4             449   3472       468  105                   95       50           14                        36                               25

# bring RNA as main assay fro processing
DefaultAssay(hyp.data) <- "RNA"

# convert to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
hyp.data.sce <- as.SingleCellExperiment(hyp.data)
save(hyp.data.sce, file = paste(Sys.Date(),"Hajdarovic_GSE188646_SingleCellExperimnent_object.RData",sep = "_"))
###############################################################################################


###############################################################################################
# 1. Run muscat for pseudobulking and extraction of samples

######################################################
####### Data preparation   ++++   Hajdarovic   #######
######################################################

hyp.data.sce.cl <- prepSCE(hyp.data.sce, 
                            kid    = "group"          ,  # population assignments
                            gid    = "stim"           ,  # group IDs (ctrl/stim)
                            sid    = "orig.ident"     ,  # sample IDs (ctrl/stim.1234)
                            drop   = TRUE           )  # drop all other colData columns

# store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk  <- length(kids <- levels(hyp.data.sce.cl$cluster_id))
ns  <- length(sids <- levels(hyp.data.sce.cl$sample_id))
names(kids) <- kids; names(sids) <- sids

# Aggregation of single-cell to pseudobulk data
pb <- aggregateData(hyp.data.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# one list item per cell type
assayNames(pb)
# [1] "Oligodendrocyte"                  "Neuron"                           "Astrocyte"                        "OPC"                             
# [5] "Microglia/Macrophage"             "Tanycyte"                         "Ependymocyte"                     "Pericyte/Endothelial Cell"       
# [9] "Vascular and Leptomeningeal Cell"

# Number of cells in each sample and cell type
cell.per.samp.tab <- t(table(hyp.data.sce.cl$cluster_id, hyp.data.sce.cl$sample_id))
#          Oligodendrocyte Neuron Astrocyte  OPC Microglia/Macrophage Tanycyte Ependymocyte Pericyte/Endothelial Cell Vascular and Leptomeningeal Cell
# Aged_1             1246   4224       720  247                  353      133           16                        41                               31
# Aged_2             1030   5248       926  220                  288      142           32                        35                               36
# Aged_3              455   3091       400   93                   98       55           45                        56                               19
# Aged_4              606   3259       358  108                   96       64           10                        20                                7
# Young_1              87    567       373   22                   42       84           38                         8                                8
# Young_2             838   4255       898  208                  192      112           18                        30                               22
# Young_3             316   2824       351   90                   54       48           29                        18                               10
# Young_4             449   3472       468  105                   95       50           14                        36                               25

# extract pseudobulk information
counts.pb.tmp <- pb@assays@data

# get the genes with no reads in at least half the samples out, they mess up the algorithm
for (i in 1:length(counts.pb.tmp)) {
  my.good <- which(apply(counts.pb.tmp[[i]]>0, 1, sum) >= nrow(cell.per.samp.tab)/2) # see deseq2 vignette, need to remove too low genes
  counts.pb.tmp[[i]] <- counts.pb.tmp[[i]][my.good,]
}

# cell types with at least 10 cells in all samples (20 samples for GRZ)
celltype.qc <- colnames(cell.per.samp.tab)[colSums(cell.per.samp.tab  >= 10) == nrow(cell.per.samp.tab)]
celltype.qc
# [1] "Oligodendrocyte"      "Neuron"               "Astrocyte"            "OPC"                  "Microglia/Macrophage" "Tanycyte"            
# [7] "Ependymocyte"  

# extract pseudobulk information for samples that pass the cell number cutoff
counts.pb <- counts.pb.tmp[celltype.qc]

#### save counts
save(counts.pb, file = paste0(Sys.Date(),"_Hajdarovic_GSE188646_muscat_PB_objects_QC_Clean.RData"))
#############################################################################################

##############################################################################################
# 2. Use SVA to clean up batch effects on expression and DEseq2 for DE analysis

# clean up memory and reload only muscat PBs
load('2024-03-22_Hajdarovic_GSE188646_muscat_PB_objects_QC_Clean.RData')

# will run SVA to clean up noise

######################################################
#######   DEG analysis   ++++   Hajdarovic     #######
######################################################

# import metadata and order it

colnames(counts.pb[[1]])
# "Aged_1"  "Aged_2"  "Aged_3"  "Aged_4"  "Young_1" "Young_2" "Young_3" "Young_4"

my.meta           <- data.frame(matrix(0, 8,4))
colnames(my.meta) <- c("SampleID","Age_Group", "Age_month", "Sex")
my.meta$Age_Group <- factor(c(rep("Y",4),rep("O",4)) , levels = c("Y","O"))
my.meta$Age_month <- c(rep(3,4),rep(19,4))
my.meta$Sex       <- rep("F",8)
my.meta$SampleID  <- c("Young_1", "Young_2", "Young_3", "Young_4", "Aged_1", "Aged_2", "Aged_3", "Aged_4")
rownames(my.meta) <- my.meta$SampleID
my.meta

# reorder count tables in sensical order
for  (i in 1:length(counts.pb)) {
  counts.pb[[i]] <- counts.pb[[i]][,my.meta$SampleID]
}

names(counts.pb)[5] <- "Microglia"

counts.pb <- counts.pb[sort(names(counts.pb))]

# Create list object to receive SVA normalized counts
sva.cts        <- vector(mode = "list", length = length(counts.pb))
names(sva.cts) <- names(counts.pb)

# Create list object to receive VST normalized counts
vst.cts        <- vector(mode = "list", length = length(counts.pb))
names(vst.cts) <- names(counts.pb)

# Create list object to receive DESeq2 results
deseq.res.list        <- vector(mode = "list", length = length(counts.pb))
names(deseq.res.list) <- names(counts.pb)

# loop over pseudobulk data
for  (i in 1:length(counts.pb)) {
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_Pseudobulk_Hajdarovic_GSE188646_",names(counts.pb)[[i]])
  
  ###################################
  #######       Run SVA      #######
  
  # build design matrix
  sva.dataDesign = data.frame( row.names = my.meta$SampleID , 
                               age       = my.meta$Age_month)
  
  # Set null and alternative models (ignore batch)
  mod1    = model.matrix(~ age , data = sva.dataDesign)
  n.sv.be = num.sv(counts.pb[[i]], mod1, method="be") 
  
  # apply SVAseq algortihm
  my.svseq = svaseq(as.matrix(counts.pb[[i]]), mod1, n.sv=n.sv.be, constant = 0.1)
  
  # remove RIN and SV, preserve age and sex
  my.clean <- removeBatchEffect(log2(counts.pb[[i]] + 0.1), 
                                covariates = cbind(my.svseq$sv),
                                design     = mod1)
  
  # delog and round data for DEseq2 processing
  my.filtered.sva <- round(2^my.clean-0.1)
  
  # keep only robustly expressed genes
  sva.cts[[i]] <- my.filtered.sva
  
  # legend
  my.cols  <- rep("",nrow(my.meta))
  my.cols[my.meta$Age_Group %in% "Y"] <- "deeppink"
  my.cols[my.meta$Age_Group %in% "O"] <- "deeppink4"

  # get matrix using age as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = sva.cts[[i]],
                                colData   = my.meta,
                                design    = ~ Age_month)
  
  # run DESeq normalizations and export results
  dds.deseq <- DESeq(dds)
  
  # get DESeq2 normalized expression value
  vst.cts[[i]] <- getVarianceStabilizedData(dds.deseq)
  
  # MDS analysis
  mds.result <- cmdscale(1-cor(vst.cts[[i]],method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  x <- mds.result[, 1]
  y <- mds.result[, 2]
  
  pdf(paste0(my.outprefix,"_MDS_plot.pdf"))
  plot(x, y,
       xlab = "MDS dimension 1", ylab = "MDS dimension 2",
       main= paste0(names(counts.pb)[[i]]," (MDS)"),
       cex=3, col = my.cols, pch = 16,
       cex.lab = 1.25,
       cex.axis = 1.25, las = 1)
  dev.off()
  
  # extract gene significance by DEseq2
  res.age <- results(dds.deseq, name = "Age_month") # FC per month
  
  # exclude genes with NA FDR value
  res.age <- res.age[!is.na(res.age$padj),]
  
  # store results
  deseq.res.list[[i]]       <- data.frame(res.age)

  ### get sex dimorphic changes at FDR5
  genes.age <- rownames(res.age)[res.age$padj < 0.05]
  my.num.age <- length(genes.age)
  
  if (my.num.age > 2) {
    # heatmap drawing - only if there is at least 2 gene
    my.heatmap.out <- paste0(my.outprefix,"_AGING_Heatmap_FDR5_GENES.pdf")
    
    pdf(my.heatmap.out, onefile = F, height = 10, width = 10)
    my.heatmap.title <- paste0(names(counts.pb)[[i]], " aging significant (FDR<5%), ", my.num.age, " genes")
    pheatmap::pheatmap(vst.cts[[i]][genes.age,],
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
  
  # # output result tables of combined analysis to text files
  # my.out.ct.mat <- paste0(my.outprefix,"_AGING_VST_log2_counts_matrix.txt")
  # write.table(vst.cts[[i]], file = my.out.ct.mat , sep = "\t" , row.names = T, quote = F)
  # 
  # my.out.stats.age <- paste0(my.outprefix,"_AGING_all_genes_statistics.txt")
  # write.table(deseq.res.list[[i]], file = my.out.stats.age , sep = "\t" , row.names = T, quote = F)
  
}

# save R object with all DEseq2 results
my.rdata.age <- paste0(Sys.Date(),"_PB_Hajdarovic_CellType_GSE188646_DEseq2_objects.RData")
save(deseq.res.list, deseq.res.list, file = my.rdata.age)

my.vst.age <- paste0(Sys.Date(),"_PB_Hajdarovic_CellType_GSE188646_VST_data_objects.RData")
save(vst.cts, file = my.vst.age)
##############################################################################################

#######################
sink(file = paste0(Sys.Date(),"_Hajdarovic_GSE188646_processing_session_Info.txt"))
sessionInfo()
sink()
