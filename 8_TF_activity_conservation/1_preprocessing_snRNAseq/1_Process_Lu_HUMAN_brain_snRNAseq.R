setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Species_Comparison/Human_Datasets_for_comparison/snRNA/GSE212606_RAW_Human')
options(stringsAsFactors = F)
options (future.globals.maxSize = 32000 * 1024^2)

# General use packages
library('Seurat')
library(bitops)
library(sctransform)
library(Matrix)

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
# 2025-07-10
# Process data from Lu et al, 2023 
################################################################################################

#### make seurat
mtx <- readMM('GSM6657986_gene_count.txt') 
mtx[1:4,1:4] 

dim(mtx) 
# [1]  60706 118240

gene.annot <- read.csv('GSM6657986_gene_annotation.csv')
cell.annot <- read.csv('GSM6657986_cell_annotation.csv')

dim(gene.annot) # [1] 60706     3
dim(cell.annot) # [1] 118240      9

rownames(mtx) <- gene.annot$gene_short_name 
colnames(mtx) <- cell.annot$Cell_ID

Lu.data <- CreateSeuratObject(counts = mtx , meta.data = cell.annot)
Lu.data
# An object of class Seurat 
# 60706 features across 118240 samples within 1 assay 
# Active assay: RNA (60706 features, 0 variable features)
# 2 layers present: counts, data

Lu.data.av      <- AverageExpression(Lu.data, assays = "RNA")
# View(Lu.data.av$RNA)

### 
Lu.meta <- readxl::read_xlsx('mmc6.xlsx', skip = 1)

### filter for brains whose meta data we have
Lu.data.flt <- subset(Lu.data, subset = Individual_ID %in% Lu.meta$Human_ID)
Lu.data.flt
# An object of class Seurat 
# 60706 features across 30801 samples within 1 assay 
# Active assay: RNA (60706 features, 0 variable features)
# 2 layers present: counts, data

Lu.data.flt.av      <- AverageExpression(Lu.data.flt, assays = "RNA")
View(Lu.data.flt.av$RNA)

### add meta data
Lu.data.flt$Age <- NA
Lu.data.flt$Sex <- NA

for (i in 1:length(Lu.meta$Human_ID)) {
  subj.id <- Lu.meta$Human_ID[i]
  Lu.data.flt$Age[Lu.data.flt$Individual_ID %in% subj.id] <- Lu.meta$AGE[Lu.meta$Human_ID %in% subj.id]
  Lu.data.flt$Sex[Lu.data.flt$Individual_ID %in% subj.id] <- Lu.meta$SEX[Lu.meta$Human_ID %in% subj.id]
  
}

head(Lu.data.flt@meta.data)
#                                       orig.ident nCount_RNA nFeature_RNA                               Cell_ID UMI_count Gene_count      Region
# Hippocampus_1_01.AACCGATTGCAATCCGGTCA Hippocampus       1286          882 Hippocampus_1_01.AACCGATTGCAATCCGGTCA      1286        882 Hippocampus
# Hippocampus_1_01.AACCGATTGCAATTGAGAGA Hippocampus       3083         1526 Hippocampus_1_01.AACCGATTGCAATTGAGAGA      3083       1526 Hippocampus
# Hippocampus_1_01.AACCGATTGCATTAACTTAA Hippocampus        349          213 Hippocampus_1_01.AACCGATTGCATTAACTTAA       349        213 Hippocampus
# Hippocampus_1_01.AACGACGAACTTGATGCTAT Hippocampus       1094          699 Hippocampus_1_01.AACGACGAACTTGATGCTAT      1094        699 Hippocampus
# Hippocampus_1_01.AACGACGAACCTTCATTAGA Hippocampus       7188         2712 Hippocampus_1_01.AACGACGAACCTTCATTAGA      7188       2712 Hippocampus
# Hippocampus_1_01.AACGGCCTAGCTGGTTGGTT Hippocampus      18208         4872 Hippocampus_1_01.AACGGCCTAGCTGGTTGGTT     18208       4872 Hippocampus
#                                                           Cell_type Condition Individual_ID     UMAP_1    UMAP_2 Age    Sex
# Hippocampus_1_01.AACCGATTGCAATCCGGTCA                     Microglia        WT          1304 -14.183059 -1.140423  81   MALE
# Hippocampus_1_01.AACCGATTGCAATTGAGAGA                    Astrocytes        WT          5459  -6.105235 10.676081  70 FEMALE
# Hippocampus_1_01.AACCGATTGCATTAACTTAA Cortical projection neurons 1        WT          1306   5.070580 -1.496925  83 FEMALE
# Hippocampus_1_01.AACGACGAACTTGATGCTAT              Oligodendrocytes        WT          5356  -1.870239 -1.075306  94   MALE
# Hippocampus_1_01.AACGACGAACCTTCATTAGA Cortical projection neurons 1        WT          1247  10.575827 -5.146237  94   MALE
# Hippocampus_1_01.AACGGCCTAGCTGGTTGGTT Cortical projection neurons 1        WT          1247   5.039482 -4.871882  94   MALE

Lu.data.flt.av      <- AverageExpression(Lu.data.flt, assays = "RNA")
View(Lu.data.flt.av$RNA)


# save filtered/annotated object
save(Lu.data.flt, file = paste0(Sys.Date(),"_Lu_Brain_Aging_Seurat_object.RData"))
################################################################################################


#######################
sink(file = paste(Sys.Date(),"_Lu_Human_Brain_Data_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()