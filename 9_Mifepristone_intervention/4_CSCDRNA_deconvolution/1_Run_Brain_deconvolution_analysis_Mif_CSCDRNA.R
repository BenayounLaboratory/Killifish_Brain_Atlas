setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/GR_signaling/Mifepristone/bulk_Brain_RNAseq/Deconvolution/')
options(stringsAsFactors = F)

# Loading necessary libraries
library(Seurat)
library(CSCDRNA)
library(Biobase)
library(beeswarm)
library(scales)

# 2025-07-16
# run CSCDRNA deconvolution
# https://github.com/empiricalbayes/CSCDRNA

# 2025-07-18
# increase stringency

###################################################################################################################
# 1. load Seurat object and  Build ExpressionSet with single-cell data
load('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Differential_Expression_Age/DGE_Analysis/2024-02-16_Seurat_objects_SPLIT_PER_STRAIN_withSampleID_with_Manual_annotation_FINAL.RData')
killi.brain.grz
# An object of class Seurat 
# 21160 features across 104665 samples within 1 assay 
# Active assay: RNA (21160 features, 5000 variable features)
# 3 layers present: counts, data, scale.data
# 3 dimensional reductions calculated: pca, umap, harmony

### extract necessary information
sc.counts.matrix    <- data.matrix(killi.brain.grz@assays$RNA@counts)
individual.labels   <- killi.brain.grz@meta.data$SampleID
cell.type.labels    <- killi.brain.grz@meta.data$Cell_Identity
sample.ids          <- colnames(sc.counts.matrix)

# make eset
sc.pheno  <- data.frame(check.names = FALSE, check.rows = FALSE, stringsAsFactors = FALSE, row.names = sample.ids, SubjectName = individual.labels,cellType = cell.type.labels)
sc.meta   <- data.frame(labelDescription = c("SubjectName","cellType"), row.names = c("SubjectName","cellType"))
sc.pdata  <- new("AnnotatedDataFrame",data  =  sc.pheno, varMetadata  =  sc.meta)
sc.eset   <- Biobase::ExpressionSet(assayData = sc.counts.matrix, phenoData = sc.pdata)
###################################################################################################################

###################################################################################################################
# 2. Build ExpressionSet with bulk brain RNAseq data 

#### 
bulk_mif.F.df   <- read.table("../DESeq2/2025-07-15_DEseq2_Bulk_GRZ_Brain_Aging_Mif_SVA_Females_SVA_corrected_counts_matrix.txt", sep = "\t", header = TRUE)
bulk_mif.F.eset <- Biobase::ExpressionSet(assayData = data.matrix(bulk_mif.F.df))

bulk_mif.M.df   <- read.table("../DESeq2/2025-07-15_DEseq2_Bulk_GRZ_Brain_Aging_Mif_SVA_Males_SVA_corrected_counts_matrix.txt", sep = "\t", header = TRUE)
bulk_mif.M.eset <- Biobase::ExpressionSet(assayData = data.matrix(bulk_mif.M.df))

###################################################################################################################

###################################################################################################################
# 3. format for CSCDRNA input
# individual.labels and cell.types should be in the same order as in sample.ids.

### 
F.mif.brain.analysis <- CSCD(bulk.eset  = bulk_mif.F.eset,
                             sc.eset    = sc.eset,
                             min.p      = 0.3,
                             markers    = NULL,
                             cell.types = "cellType",
                             subj.names = "SubjectName",
                             verbose    = TRUE)
# Using 941 genes in both bulk and single-cell expression.



M.mif.brain.analysis <- CSCD(bulk.eset  = bulk_mif.M.eset,
                             sc.eset    = sc.eset,
                             min.p      = 0.3,
                             markers    = NULL,
                             cell.types = "cellType",
                             subj.names = "SubjectName",
                             verbose    = TRUE)
# Using 941 genes in both bulk and single-cell expression.


save(F.mif.brain.analysis, M.mif.brain.analysis, file = paste0(Sys.Date(),"_CSCDRNA_deconvolution_results.RData"))
###########################################################################################

###################################################################################################################
# 4. extract microglia information

# individual.labels and cell.types should be in the same order as in sample.ids.
F.mif.decon.res <- data.frame(t(F.mif.brain.analysis$bulk.props))
M.mif.decon.res <- data.frame(t(M.mif.brain.analysis$bulk.props))

microglia.props <- list("YF_CTL" = F.mif.decon.res$Microglia[1:3],
                        "OF_CTL" = F.mif.decon.res$Microglia[4:6],
                        "OF_MIF" = F.mif.decon.res$Microglia[7:8],
                        "YM_CTL" = M.mif.decon.res$Microglia[1:4],
                        "OM_CTL" = M.mif.decon.res$Microglia[5:7],
                        "OM_MIF" = M.mif.decon.res$Microglia[8:10])

pdf(paste0(Sys.Date(),"_microglia_frequency_CSCDRNA_boxplot.pdf"))
boxplot(microglia.props,
        outline = F, col = c("deeppink","deeppink4","hotpink2","deepskyblue","deepskyblue4","lightskyblue"),
        ylim = c(0,0.03), las = 2, ylab = "Fraction microglia (CSCDRNA)")
beeswarm::beeswarm(microglia.props, pch = 16, add = T)
dev.off()

wilcox.test(microglia.props$YF_CTL, microglia.props$OF_CTL)    # p-value = 0.1
wilcox.test(microglia.props$OF_CTL, microglia.props$OF_MIF)    # p-value = 1
wilcox.test(microglia.props$YM_CTL, microglia.props$OM_CTL)    # p-value = 0.4
wilcox.test(microglia.props$OM_CTL, microglia.props$OM_MIF)    # p-value = 0.4

microglia.props.v2 <- list("Y_CTL" = c(F.mif.decon.res$Microglia[1:3], M.mif.decon.res$Microglia[1:4]  ),
                           "O_CTL" = c(F.mif.decon.res$Microglia[4:6],M.mif.decon.res$Microglia[5:7]  ),
                           "O_MIF" = c(F.mif.decon.res$Microglia[7:8],M.mif.decon.res$Microglia[8:10] ) )


boxplot(microglia.props.v2,
        outline = F, col = c("grey","grey25","grey50"),
        ylim = c(0,0.03), las = 2, ylab = "Fraction microglia (CSCDRNA)")
beeswarm::beeswarm(microglia.props.v2, pch = 16, add = T)

test.age <- wilcox.test(microglia.props.v2$Y_CTL, microglia.props.v2$O_CTL)    # p-value = 0.03497
test.mif <- wilcox.test(microglia.props.v2$O_CTL, microglia.props.v2$O_MIF)    # p-value = 0.1775
test.rev <- wilcox.test(microglia.props.v2$Y_CTL, microglia.props.v2$O_MIF)    # p-value = 0.6389


pdf(paste0(Sys.Date(),"_microglia_frequency_CSCDRNA_boxplot.pdf"), width = 3.5, height = 5)
boxplot(microglia.props.v2,
        outline = F, col = c("grey","grey25","grey50"),
        ylim = c(0,0.03), las = 2, ylab = "Fraction microglia (CSCDRNA)")
beeswarm(microglia.props.v2, pch = 16, add = T,
         pwcol = c(rep("deepskyblue" , 4),
                   rep("deeppink"    , 3),
                   rep("deeppink4"   , 3),
                   rep("deepskyblue4", 3),
                   rep("hotpink2"    , 2),
                   rep("lightskyblue", 3)))
text(1.5, 0.025, scientific(test.age$p.value,2))
text(2.5, 0.025, scientific(test.mif$p.value,2))
text(2  , 0.03, scientific(test.rev$p.value,2))

dev.off()



write.table(rbind(F.mif.decon.res,
                  M.mif.decon.res), file = paste0(Sys.Date(),"_CSCDRNA_KilliBrain_Aging_Mif_analysis.txt"), quote = FALSE, sep = "\t")
###################################################################################################################

#######################
sink(file = paste(Sys.Date(),"_Killi_brain_CSCDRNA_session_Info.txt", sep =""))
sessionInfo()
sink()
