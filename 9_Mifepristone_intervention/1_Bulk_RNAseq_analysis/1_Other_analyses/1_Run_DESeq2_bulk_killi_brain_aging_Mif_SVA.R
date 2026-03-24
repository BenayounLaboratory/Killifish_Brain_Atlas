setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/GR_signaling/Mifepristone/bulk_Brain_RNAseq/DESeq2')
options(stringsAsFactors = F)

library("DESeq2")        #
library("sva")           #
library("limma")         #
library("pheatmap")      #
library("bitops")      #
library(phenoTest)

library(ggplot2) 
library(scales) 
theme_set(theme_bw())


# 2025-01-24
# Analyze pilot dataset of mid-life Mifepristone treatment

# 2025-07-15
# analyze the 2 cohorts, with trimmed/remapped count tables
# use SVA to limit technical noise

source('1b_Mif_RNAseq_analysis_functions_SVA.R')

###################################################################################################################
# 1. Read and clean up count matrices

#### FEMALES (Cohort 1)
f.mif.cts            <- read.table("../STAR/2025-07-15_Killifish_Brain_Aging_GRZ_Females_Mif_RNAseq_counts.txt", sep = "\t", header = TRUE)
f.mif.cts            <- f.mif.cts[,-(2:6)] # remove superfluous annotation columns
f.mif.cts$Geneid     <- gsub("gene-","", f.mif.cts$Geneid,fixed = TRUE)
rownames(f.mif.cts)  <- f.mif.cts$Geneid
colnames(f.mif.cts)  <- gsub("_STAR_Aligned.sortedByCoord.out.bam","",colnames(f.mif.cts))
colnames(f.mif.cts)  <- gsub("GRZ_","",colnames(f.mif.cts))
colnames(f.mif.cts)  <- gsub("_Brain_c1_","",colnames(f.mif.cts))

head(f.mif.cts)
#                    Geneid FY_CTL_rep1 FY_CTL_rep2 FY_CTL_rep3 FO_CTL_rep1 FO_CTL_rep2 FO_CTL_rep3 FO_MIF_rep1 FO_MIF_rep2
# LOC107382895 LOC107382895       268.5       212.5       379.5       301.0      371.00      214.00       430.0      212.00
# LOC107382813 LOC107382813       473.5       360.5       467.0       470.0      390.50      342.00       493.0      416.50


#### MALES (Cohort 2)
m.mif.cts <- read.table("../STAR/2025-07-15_Killifish_Brain_Aging_GRZ_Males_Mif_RNAseq_counts.txt", sep = "\t", header = TRUE)
m.mif.cts            <- m.mif.cts[,-(2:6)] # remove superfluous annotation columns
m.mif.cts$Geneid     <- gsub("gene-","", m.mif.cts$Geneid,fixed = TRUE)
rownames(m.mif.cts)  <- m.mif.cts$Geneid
colnames(m.mif.cts)  <- gsub("_STAR_Aligned.sortedByCoord.out.bam","",colnames(m.mif.cts))
colnames(m.mif.cts)  <- gsub("GRZ_","",colnames(m.mif.cts))
colnames(m.mif.cts)  <- gsub("_Brain_c2_","",colnames(m.mif.cts))

head(m.mif.cts)
#                   Geneid MY_CTL_rep1 MY_CTL_rep2 MY_CTL_rep3 MY_CTL_rep4 MO_CTL_rep1 MO_CTL_rep2 MO_CTL_rep3 MO_MIF_rep1 MO_MIF_rep2 MO_MIF_rep3
# LOC107382895 LOC107382895       128.0      153.50      156.50       274.0       195.0      101.50       142.0       185.5       94.50      183.00
# LOC107382813 LOC107382813       128.5      194.00      210.00       199.5       138.5       90.00       148.5       150.0      199.50       84.00
###################################################################################################################

###################################################################################################################
# 2. Run analyses
F.mif.res <- process_mif_data(f.mif.cts, sex = "Females")
  
M.mif.res <- process_mif_data(m.mif.cts, sex = "Males")

save(F.mif.res, M.mif.res, file = paste0(Sys.Date(), "_Brain_MIF_DEseq2_results.RData"))

###################################################################################################################

###################################################################################################################
# 3. Export results
options(java.parameters = "-Xmx16g" )
require(openxlsx)

names(F.mif.res)[1:2] <- c("F_Aging", "F_Mifepristone")
names(M.mif.res)[1:2] <- c("M_Aging", "M_Mifepristone")

write.xlsx(c(F.mif.res[1:2],M.mif.res[1:2]), rowNames = TRUE, file = paste0(Sys.Date(),"_GRZ_Mifepristone_Aging_DESeq2_Results.xlsx"))
###################################################################################################################


#######################
sink(file = paste(Sys.Date(),"_DEseq2_bulk_killi_BrainAging_Mifepristone_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()
