setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/GR_signaling/Mifepristone/bulk_Brain_RNAseq/DESeq2')
options(stringsAsFactors = F)

library("DESeq2")        #
library("sva")           #
library("limma")         #
library("pheatmap")      #
library("bitops")        #
library(phenoTest)

library(ggplot2) 
library(scales) 
theme_set(theme_bw())

library(beeswarm)

# 2025-07-16
# get microglia genes

#####################################################################################################################
#### 1. Run GSEA

load('2025-07-15_Brain_MIF_DEseq2_results.RData')
ls()
# [1] "F.mif.res" "M.mif.res"

#####################################################################################################################
# 1. microglia markers in female dataset
F.apoeb.cts <- F.mif.res$VST["LOC107379395",]

F.apoeb.cts.list <- list("YF_CTL" = F.apoeb.cts[1:3],
                         "OF_CTL" = F.apoeb.cts[4:6],
                         "OF_MIF" = F.apoeb.cts[7:8])

F.apoeb.res.1 <- F.mif.res$Aging["LOC107379395",]
F.apoeb.res.2 <- F.mif.res$Mif  ["LOC107379395",]

pdf(paste0(Sys.Date(),"_Apoeb_Microglia_markers_gene_exp_FEMALE.pdf"),height = 5, width=3.5)
boxplot(F.apoeb.cts.list, col = c("deeppink","deeppink4","hotpink2"), 
        main = "GRZ brain (F)", ylab = "apoeb/LOC107379395 DESeq2 log2 VST counts", ylim = c(11,14), las = 2)
beeswarm(F.apoeb.cts.list, add = T, pch = 16)
text(1.5,14, scientific(F.apoeb.res.1$padj,2))
text(2.5,14, scientific(F.apoeb.res.2$padj,2))
dev.off()

#
F.csf1r.cts <- F.mif.res$VST["csf1r",]

F.csf1r.cts.list <- list("YF_CTL" = F.csf1r.cts[1:3],
                         "OF_CTL" = F.csf1r.cts[4:6],
                         "OF_MIF" = F.csf1r.cts[7:8])

F.csf1r.res.1 <- F.mif.res$Aging["csf1r",]
F.csf1r.res.2 <- F.mif.res$Mif  ["csf1r",]

pdf(paste0(Sys.Date(),"_csf1r_Microglia_markers_gene_exp_FEMALE.pdf"),height = 5, width=3.5)
boxplot(F.csf1r.cts.list, col = c("deeppink","deeppink4","hotpink2"), 
        main = "GRZ brain (F)", ylab = "csf1r DESeq2 log2 VST counts", ylim = c(8,11), las = 2)
beeswarm(F.csf1r.cts.list, add = T, pch = 16)
text(1.5,11, scientific(F.csf1r.res.1$padj,2))
text(2.5,11, scientific(F.csf1r.res.2$padj,2))
dev.off()



##############################################################################################################
# 2. microglia markers in male dataset

M.apoeb.cts <- M.mif.res$VST["LOC107379395",]

M.apoeb.cts.list <- list("YM_CTL" = M.apoeb.cts[1:4],
                         "OM_CTL" = M.apoeb.cts[5:7],
                         "OM_MIF" = M.apoeb.cts[8:10])


M.apoeb.res.1 <- M.mif.res$Aging["LOC107379395",]
M.apoeb.res.2 <- M.mif.res$Mif  ["LOC107379395",]


pdf(paste0(Sys.Date(),"_ApoebMicroglia_markers_gene_exp_MALE.pdf"), height = 5, width=3.5)
boxplot(M.apoeb.cts.list, col = c("deepskyblue","deepskyblue4","lightskyblue"), 
        main = "GRZ brain (M)", ylab = "apoeb/LOC107379395 DESeq2 log2 VST counts", ylim = c(11,14), las = 2)
beeswarm(M.apoeb.cts.list, add = T, pch = 16)
text(1.5,14, scientific(M.apoeb.res.1$padj,2))
text(2.5,14, scientific(M.apoeb.res.2$padj,2))
dev.off()


# 
M.csf1r.cts <- M.mif.res$VST["csf1r",]

M.csf1r.cts.list <- list("YM_CTL" = M.csf1r.cts[1:4],
                         "OM_CTL" = M.csf1r.cts[5:7],
                         "OM_MIF" = M.csf1r.cts[8:10])


M.csf1r.res.1 <- M.mif.res$Aging["csf1r",]
M.csf1r.res.2 <- M.mif.res$Mif  ["csf1r",]


pdf(paste0(Sys.Date(),"_csf1r_Microglia_markers_gene_exp_MALE.pdf"), height = 5, width=3.5)
boxplot(M.csf1r.cts.list, col = c("deepskyblue","deepskyblue4","lightskyblue"), 
        main = "GRZ brain (M)", ylab = "csf1r DESeq2 log2 VST counts", ylim = c(8,11), las = 2)
beeswarm(M.csf1r.cts.list, add = T, pch = 16)
text(1.5,11, scientific(M.csf1r.res.1$padj,2))
text(2.5,11, scientific(M.csf1r.res.2$padj,2))
dev.off()


