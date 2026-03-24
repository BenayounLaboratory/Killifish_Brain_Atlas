setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/ATAC/Bulk_ATAC_Analysis_Folder/DESeq2')
options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('pvclust')
library('bitops')
library('sva')
library('limma')
library(RColorBrewer)
library(fields)

# 2024-05-01
# Killifish brain bulk ATAC aging
# analyze overlap of DA peaks with H3K4me3 ChIPseq


my.outprefix.zmz <- paste(Sys.Date(),"ZMZ_brain_Aging_ATAC_DESeq2_Analysis",sep="_")
my.outprefix.grz <- paste(Sys.Date(),"GRZ_brain_Aging_ATAC_DESeq2_Analysis",sep="_")

#####################################################
#######   H3K4me3  analysis   ++++   GRZ     #######
#####################################################
load('2024-04-09_GRZ_brain_Aging_ATAC_BOTH.RData')

# get H3K4me3 overlaps
my.peak.annot                <- read.csv('../Diffbind/GRZ_Aging_Brains_MACS2_MSPC_diffbind_peaks_H3K4me3_overlap.bed', sep = "\t", header = F)
my.peak.annot                <- my.peak.annot[,1:8]
colnames(my.peak.annot)      <- c("Chr", "Start", "End","Chr_K4", "Start_K4", "End_K4","PeakID_K4","Overlap_K4")
my.peak.annot$PeakName       <- paste(my.peak.annot$Chr,my.peak.annot$Start,my.peak.annot$End,sep = ":")
my.peak.annot$H3K4me3_Status <- ifelse(my.peak.annot$PeakID_K4 == ".", "Not_H3K4me3","H3K4me3")

##### Output AGE bed
res.age.grz$PeakName <- rownames(res.age.grz)
res.age.grz          <- data.frame(res.age.grz)
res.age.grz.annot    <- merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "PeakID_K4",  "H3K4me3_Status")],res.age.grz, by = "PeakName")


background.grz <- res.age.grz.annot
age.up.grz     <- res.age.grz.annot[bitAnd(res.age.grz.annot$padj < 0.05, res.age.grz.annot$log2FoldChange > 0 )>0,]
age.dwn.grz    <- res.age.grz.annot[bitAnd(res.age.grz.annot$padj < 0.05, res.age.grz.annot$log2FoldChange < 0 )>0,]

cts.background.grz <- aggregate(background.grz$H3K4me3_Status, by = list(background.grz$H3K4me3_Status), FUN = length)
cts.age.up.grz     <- aggregate(age.up.grz  $H3K4me3_Status  , by = list(age.up.grz$H3K4me3_Status    ), FUN = length)
cts.age.dwn.grz    <- aggregate(age.dwn.grz $H3K4me3_Status  , by = list(age.dwn.grz$H3K4me3_Status   ), FUN = length)

pdf(paste0(my.outprefix.grz,"_PieCharts_H3K4me3_FDR5.pdf"), height = 3, width = 9)
par(mfrow=c(1,3) )
pie(as.numeric(cts.background.grz$x), labels = cts.background.grz$Group.1, main = "All ATAC (GRZ)" , col = c("goldenrod1","palevioletred1") )
pie(as.numeric(cts.age.up.grz$x    ), labels = cts.age.up.grz$Group.1    , main = "ATAC Up (GRZ)"  , col = c("goldenrod1","palevioletred1")  )
pie(as.numeric(cts.age.dwn.grz$x   ), labels = cts.age.dwn.grz$Group.1   , main = "ATAC Down (GRZ)", col = c("goldenrod1","palevioletred1")  )
par(mfrow=c(1,1) )
dev.off()

chisq.up <- chisq.test(rbind(cts.age.up.grz$x,cts.background.grz$x))
chisq.up$p.value # 1.040293e-122

chisq.dwn <- chisq.test(rbind(cts.age.dwn.grz$x,cts.background.grz$x))
chisq.dwn$p.value #  1.1374e-170


#####################################################
#######   H3K4me3  analysis   ++++   ZMZ      #######
#####################################################
load('2024-04-09_ZMZ_brain_Aging_ATAC_BOTH.RData')

# get H3K4me3 overlaps
my.peak.annot                <- read.csv('../Diffbind/ZMZ_Aging_Brains_MACS2_MSPC_diffbind_peaks_H3K4me3_overlap.bed', sep = "\t", header = F)
my.peak.annot                <- my.peak.annot[,1:8]
colnames(my.peak.annot)      <- c("Chr", "Start", "End","Chr_K4", "Start_K4", "End_K4","PeakID_K4","Overlap_K4")
my.peak.annot$PeakName       <- paste(my.peak.annot$Chr,my.peak.annot$Start,my.peak.annot$End,sep = ":")
my.peak.annot$H3K4me3_Status <- ifelse(my.peak.annot$PeakID_K4 == ".", "Not_H3K4me3","H3K4me3")

##### Output AGE bed
res.age.zmz$PeakName <- rownames(res.age.zmz)
res.age.zmz          <- data.frame(res.age.zmz)
res.age.zmz.annot    <- merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "PeakID_K4",  "H3K4me3_Status")],res.age.zmz, by = "PeakName")


background.zmz <- res.age.zmz.annot
age.up.zmz     <- res.age.zmz.annot[bitAnd(res.age.zmz.annot$padj < 0.05, res.age.zmz.annot$log2FoldChange > 0 )>0,]
age.dwn.zmz    <- res.age.zmz.annot[bitAnd(res.age.zmz.annot$padj < 0.05, res.age.zmz.annot$log2FoldChange < 0 )>0,]

cts.background.zmz <- aggregate(background.zmz$H3K4me3_Status, by = list(background.zmz$H3K4me3_Status), FUN = length)
cts.age.up.zmz     <- aggregate(age.up.zmz  $H3K4me3_Status  , by = list(age.up.zmz$H3K4me3_Status    ), FUN = length)
cts.age.dwn.zmz    <- aggregate(age.dwn.zmz $H3K4me3_Status  , by = list(age.dwn.zmz$H3K4me3_Status   ), FUN = length)

pdf(paste0(my.outprefix.zmz,"_PieCharts_H3K4me3_FDR5.pdf"), height = 3, width = 9)
par(mfrow=c(1,3) )
pie(as.numeric(cts.background.zmz$x), labels = cts.background.zmz$Group.1, main = "All ATAC (ZMZ)" , col = c("goldenrod1","palevioletred1") )
pie(as.numeric(cts.age.up.zmz$x    ), labels = cts.age.up.zmz$Group.1    , main = "ATAC Up (ZMZ)"  , col = c("goldenrod1","palevioletred1")  )
pie(as.numeric(cts.age.dwn.zmz$x   ), labels = cts.age.dwn.zmz$Group.1   , main = "ATAC Down (ZMZ)", col = c("goldenrod1","palevioletred1")  )
par(mfrow=c(1,1) )
dev.off()

chisq.up <- chisq.test(rbind(cts.age.up.zmz$x,cts.background.zmz$x))
chisq.up$p.value # 0.2869419

chisq.dwn <- chisq.test(rbind(cts.age.dwn.zmz$x,cts.background.zmz$x))
chisq.dwn$p.value #  1.82126e-100
################################################################################################



#######################
sink(file = paste0(Sys.Date(),"_Killi_Brain_aging_bulkATAC_H4K4me3_comp_session_Info.txt"))
sessionInfo()
sink()
