setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/ATAC/Bulk_ATAC_Analysis_Folder/Diffbind')
options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library('DiffBind')
library("rtracklayer")

# 2024-04-06
# Diffbind on GRZ/ZMZ brain ATAC (MSPC consensus peaks)
# use MSPC peaks with 20/47 sample peak support
# need to remove the score of the MSPC peak file as "inf" makes diffbind error out

################################################################################
####################    Diffbind Analysis GRZ brain aging   ####################
################################################################################
grz.brain.aging <- dba(sampleSheet="ATAC_Brain_Aging_GRZ_samples.csv",skipLines=0,attributes=c(DBA_ID,DBA_CONDITION))
grz.brain.aging <- dba.count(grz.brain.aging)

pdf(paste0(Sys.Date(),"GRZ_Aging_Brains_heatmap_2cohorts.pdf"))
plot(grz.brain.aging, colScheme="Reds")
dev.off()

grz.brain.aging <- dba.contrast(grz.brain.aging, categories=DBA_CONDITION,minMembers=2)
grz.brain.aging <- dba.analyze(grz.brain.aging,method=DBA_ALL_METHODS, design = F)
grz.brain.aging
# 24 Samples, 73378 sites in matrix:
#   ID Condition Replicate     Reads FRiP
# 1  OF1        OF         1  47531090 0.42
# 2  OF2        OF         2  28360666 0.36
# 3  OF3        OF         3  40764983 0.36
# 4  OM1        OM         1  57187320 0.42
# 5  OM2        OM         2  52266024 0.44
# 6  OM3        OM         3  54067912 0.41
# 7  YF1        YF         1  64291734 0.33
# 8  YF2        YF         2  34906227 0.35
# 9  YF3        YF         3  56572297 0.42
# 10 YM1        YM         1  56020243 0.30
# 11 YM2        YM         2  44249788 0.33
# 12 YM3        YM         3  36929659 0.41
# 13 OF4        OF         4 102566129 0.37
# 14 OF5        OF         5  97010453 0.34
# 15 OF6        OF         6  93773974 0.34
# 16 OM4        OM         4  77353489 0.38
# 17 OM5        OM         5  59016597 0.35
# 18 OM6        OM         6  69434326 0.25
# 19 YF4        YF         4 110912138 0.32
# 20 YF5        YF         5  64098224 0.34
# 21 YF6        YF         6 102390755 0.40
# 22 YM4        YM         4  80324679 0.39
# 23 YM5        YM         5  81969794 0.38
# 24 YM6        YM         6 108067710 0.31
# 
# 6 Contrasts:
#   Group Samples Group2 Samples2 DB.edgeR
# 1    OF       6     OM        6       15
# 2    OF       6     YF        6        1
# 3    OF       6     YM        6        5
# 4    OM       6     YF        6       26
# 5    OM       6     YM        6        0
# 6    YM       6     YF        6        3

save(grz.brain.aging, file = paste0(Sys.Date(),"ATAC_GRZ_Aging_Brains_Diffbind_MACS2_MSPC.RData"))


my.norm.counts <- dba.peakset(grz.brain.aging, bRetrieve=TRUE)

write.table(data.frame(my.norm.counts), 
            file = paste0(Sys.Date(),"_GRZ_Aging_Brains_MACS2_MSPC_Normalized_count_matrix.txt"),
            quote = F, col.names = T, row.names = F, sep = "\t")

write.table(data.frame(my.norm.counts)[,c(1:3)], 
            file = paste0(Sys.Date(),"_GRZ_Aging_Brains_MACS2_MSPC_diffbind_peaks.bed"),
            quote = F, col.names = F, row.names = F, sep = "\t")


################################################################################
####################    Diffbind Analysis ZMZ brain aging   ####################
################################################################################
zmz.brain.aging <- dba(sampleSheet="ATAC_Brain_Aging_ZMZ_samples.csv",skipLines=0,attributes=c(DBA_ID,DBA_CONDITION))
zmz.brain.aging <- dba.count(zmz.brain.aging)

pdf(paste0(Sys.Date(),"ZMZ_Aging_Brains_heatmap_2cohorts.pdf"))
plot(zmz.brain.aging, colScheme="Reds")
dev.off()

zmz.brain.aging <- dba.contrast(zmz.brain.aging, categories=DBA_CONDITION,minMembers=2)
zmz.brain.aging <- dba.analyze(zmz.brain.aging,method=DBA_ALL_METHODS, design = F)
zmz.brain.aging
# 23 Samples, 71199 sites in matrix:
#   ID Condition Replicate    Reads FRiP
# 1  GF1        GF         1 93683515 0.44
# 2  GF2        GF         2 17919462 0.20
# 3  GM1        GM         1 28564685 0.12
# 4  GM2        GM         2 44343164 0.39
# 5  OF1        OF         1 82131059 0.33
# 6  OM1        OM         1 43520774 0.45
# 7  OM2        OM         2 20716650 0.18
# 8  YF1        YF         1 65012260 0.34
# 9  YF2        YF         2 15046039 0.15
# 10 YM1        YM         1 45382344 0.37
# 11 YM2        YM         2 48769751 0.38
# 12 GF3        GF         3 33150906 0.40
# 13 GF4        GF         4 58178669 0.47
# 14 GM3        GM         3 33575110 0.48
# 15 GM4        GM         4 27620796 0.48
# 16 OF3        OF         3 53916675 0.45
# 17 OF4        OF         4 49782289 0.44
# 18 OM3        OM         3 40514545 0.45
# 19 OM4        OM         4 48118535 0.49
# 20 YF3        YF         3 33418844 0.43
# 21 YF4        YF         4 33828695 0.43
# 22 YM3        YM         3 31944678 0.45
# 23 YM4        YM         4 41760287 0.41
# 
# Group Samples Group2 Samples2 DB.edgeR
# 1     GF       4     GM        4        3
# 2     GF       4     OF        3        3
# 3     GF       4     OM        4        5
# 4     GF       4     YF        4       10
# 5     GF       4     YM        4      283
# 6     GM       4     OF        3       12
# 7     GM       4     OM        4        1
# 8     GM       4     YF        4        8
# 9     GM       4     YM        4        4
# 10    OF       3     OM        4       10
# 11    OF       3     YF        4        4
# 12    OF       3     YM        4       47
# 13    OM       4     YF        4        9
# 14    OM       4     YM        4       12
# 15    YM       4     YF        4       11


save(zmz.brain.aging, file = paste0(Sys.Date(),"ATAC_ZMZ_Aging_Brains_Diffbind_MACS2_MSPC.RData"))


my.norm.counts <- dba.peakset(zmz.brain.aging,bRetrieve=TRUE)

write.table(data.frame(my.norm.counts), 
            file = paste0(Sys.Date(),"_ZMZ_Aging_Brains_MACS2_MSPC_Normalized_count_matrix.txt"),
            quote = F, col.names = T, row.names = F, sep = "\t")

write.table(data.frame(my.norm.counts)[,c(1:3)], 
            file = paste0(Sys.Date(),"_ZMZ_Aging_Brains_MACS2_MSPC_diffbind_peaks.bed"),
            quote = F, col.names = F, row.names = F, sep = "\t")


###############################################################################################################################################################
sink(file = paste(Sys.Date(),"_Killi_Brain_Aging_ATACseq_Diffbind_session_Info.txt", sep =""))
sessionInfo()
sink()
