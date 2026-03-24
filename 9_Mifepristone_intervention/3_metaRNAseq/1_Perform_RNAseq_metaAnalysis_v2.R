setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/GR_signaling/Mifepristone/bulk_Brain_RNAseq/metaRNAseq/')
options(stringsAsFactors = F)

library("DESeq2")        #
library('metaRNASeq')    #
library("pheatmap")      #

library(Vennerable)
library(bitops)

library(clusterProfiler)   #
library(ggplot2)           #
library(scales)            #

theme_set(theme_bw())

# 2025-07-28
# Perform meta analysis of mif/aging effects for females/males (Cohort1/Cohort2)
# include genes and TE to match base analysis

#####################################################################################################################
#### 1. Load DEseq2 results

load('../DESeq2/2025-07-15_Brain_MIF_DEseq2_results.RData')
# F.mif.res
# M.mif.res

### Filter for genes detected in both cohorts (since meta-analysis of other genes is nonsensical)
genes.F      <- intersect(rownames(F.mif.res$Aging), rownames(F.mif.res$Mif)) # 18005
genes.M      <- intersect(rownames(M.mif.res$Aging), rownames(M.mif.res$Mif)) # 22805
genes.common <- intersect(genes.F,genes.M) # 18001


# filter results accordingly
F.mif.res$Aging <- F.mif.res$Aging[genes.common, ]
M.mif.res$Aging <- M.mif.res$Aging[genes.common, ]

F.mif.res$Mif <- F.mif.res$Mif[genes.common, ]
M.mif.res$Mif <- M.mif.res$Mif[genes.common, ]
#####################################################################################################################

#####################################################################################################################
#### 2. AGING   meta analysis

# We recommand to store both p-value and Fold Change results in lists in order to perform
# meta-analysis and keep track of the potential conflicts (see section 5)
age.rawpval <- list("pval.F" = F.mif.res$Aging$pvalue, 
                    "pval.M" = M.mif.res$Aging$pvalue)
age.FC      <- list("FC.F"   = F.mif.res$Aging$log2FoldChange,
                    "FC.M"   = M.mif.res$Aging$log2FoldChange)

# Differentially expressed genes in each individual study can also be marked in a matrix DE:
age.adjpval <- list("adjpval.F" = F.mif.res$Aging$padj,
                    "adjpval.M" = M.mif.res$Aging$padj)

age.studies      <- c("FEMALES", "MALES")
age.DE           <- data.frame(mapply(age.adjpval, FUN=function(x) ifelse(x <= 0.05, 1, 0)))
colnames(age.DE) <- paste("DE", age.studies, sep=".")
rownames(age.DE) <- rownames(F.mif.res$Aging)

# DE returns a list with 1 for genes identified as differentially expressed and 0 otherwise
# Since the proposed p-value combination techniques rely on the assumption that p-values follow a uniform distribution under the null hypothesis, 
# it is necessary to check that the histograms of raw-pvalues reflect that assumption
# The peak near 0 corresponds to differentially expressed genes, no other peak should appear
pdf(paste0(Sys.Date(), "_rawp_histogram_aging.pdf"), height = 5, width = 8)
par(mfrow = c(1,2))
hist(age.rawpval[[1]], breaks=100, col="grey", main="Female aging", xlab="Raw p-value")
hist(age.rawpval[[2]], breaks=100, col="grey", main="Male aging", xlab="Raw p-value")
par(mfrow = c(1,1))
dev.off()

# Account for DESeq2 automatic independent filtering. 
age.filtered <- lapply(age.adjpval, FUN=function(pval) which(is.na(pval)))
age.rawpval[[1]][age.filtered[[1]]]=NA
age.rawpval[[2]][age.filtered[[2]]]=NA

# Check raw p-values are roughly uniformly distributed
pdf(paste0(Sys.Date(), "_rawp_histogram_aging_FILTERED.pdf"), height = 5, width = 8)
par(mfrow = c(1,2))
hist(age.rawpval[[1]], breaks=100, col="grey", main="Female aging", xlab="Raw p-value")
hist(age.rawpval[[2]], breaks=100, col="grey", main="Male aging", xlab="Raw p-value")
par(mfrow = c(1,1))
dev.off()

# Use of p-value combination techniques
# The p-value combination using the Fisher method may be performed with the fishercomb
# function, and the subsequent p-values obtained from the meta-analysis may be examined (Figure 3, left):
age.fishcomb    <- fishercomb(age.rawpval, BHth = 0.05)
age.invnormcomb <- invnorm(age.rawpval, nrep=c(6,7), BHth = 0.05) # F: 3Y, 3O; M: 4Y, 3O

pdf(paste0(Sys.Date(), "_META_pvalue_histogram_aging.pdf"), height = 5, width = 8)
par(mfrow = c(1,2))
hist(age.fishcomb$rawpval   , breaks=100, col="grey", main="Aging - Fisher method"        , xlab = "Raw p-values (meta-analysis)")
hist(age.invnormcomb$rawpval, breaks=100, col="grey", main="Aging - Inverse normal method", xlab = "Raw p-values (meta-analysis)")
par(mfrow = c(1,1))
dev.off()

# summarize the results of the individual differential analyses as well
# as the differential meta-analysis (using the Fisher and inverse normal methods) in a data.frame:
age.DEresults <- data.frame(age.DE,
                            "DE.fishercomb"= ifelse(age.fishcomb$adjpval<=0.05,1,0),
                            "DE.invnorm"   = ifelse(age.invnormcomb$adjpval<=0.05,1,0))
head(age.DEresults)
#              DE.FEMALES DE.MALES DE.fishercomb DE.invnorm
# LOC107382895          0        0             0          0
# LOC107382813          0        1             1          1
# fmr1                  0        1             1          1
# aff2                  0        0             0          0
# ids                   1        1             1          1
# dcps                  0        0             0          0

# Treatment of conflicts in differential expression
# We build a matrix signsFC gathering all signs of fold changes from individual studies.
age.signsFC      <- mapply(age.FC, FUN=function(x) sign(x))
age.sumsigns     <- apply(age.signsFC,1,sum)
age.commonsgnFC  <- ifelse(abs(age.sumsigns)==dim(age.signsFC)[2], sign(age.sumsigns),0)

# The vector commonsgnFC returns:
#    - a value of 1 if the gene has a positive log2 fold change in all studies, 
#    - -1 if the gene has a negative log2 fold change in all studies,
#    - 0 if contradictory log2 fold changes are observed across studies (i.e., positive in one and negative in the other). 
# By examining the elements of commonsgnFC, it is thus possible to identify genes displaying contradictory differential expression among studies.
age.unionDE     <- unique(c(age.fishcomb$DEindices,age.invnormcomb$DEindices))
age.FC.selecDE  <- data.frame(age.DEresults[age.unionDE,],
                              do.call(cbind,age.FC)[age.unionDE],
                              signFC=age.commonsgnFC[age.unionDE])
age.keepDE     <- age.FC.selecDE[which(abs(age.FC.selecDE$signFC)==1),]
age.conflictDE <- age.FC.selecDE[which(age.FC.selecDE$signFC == 0),]
dim(age.FC.selecDE) # 7605    6
dim(age.keepDE    ) # 5469    6
dim(age.conflictDE) # 2136    6

# IDD, IRR and Venn Diagram
age.fishcomb_de <- rownames(age.keepDE)[which(age.keepDE[,"DE.fishercomb"]==1)]
age.invnorm_de  <- rownames(age.keepDE)[which(age.keepDE[,"DE.invnorm"]==1)]
age.indstudy_de <- list(rownames(age.keepDE)[which(age.keepDE[,"DE.FEMALES"]==1)], rownames(age.keepDE)[which(age.keepDE[,"DE.MALES"]==1)])
IDD.IRR(age.fishcomb_de,age.indstudy_de)
# DE     IDD    Loss     IDR     IRR 
#5342.00  469.00    4.00    8.78    0.08 
IDD.IRR(age.invnorm_de ,age.indstudy_de)
# DE     IDD    Loss     IDR     IRR 
# 4958.00  575.00  494.00   11.60   10.13  

length(intersect(age.fishcomb_de,age.invnorm_de)) # 4479

age.meta.results <- cbind(age.DEresults,
                          "meta_sign"   = age.commonsgnFC,
                          "Fisher_FDR"  = age.fishcomb$adjpval   ,
                          "InvNorm_FDR" = age.invnormcomb$adjpval)
# sum(age.meta.results$InvNorm_FDR[abs(age.meta.results$meta_sign) == 1] < 0.05) # 4958

save(age.meta.results, file = paste0(Sys.Date(),"_AGING_meta_analysis_results.RData"))
#####################################################################################################################


#####################################################################################################################
#### 3. MIFEPRISTONE   meta analysis

# We recommand to store both p-value and Fold Change results in lists in order to perform
# meta-analysis and keep track of the potential conflicts (see section 5)
mif.rawpval <- list("pval.F" = F.mif.res$Mif$pvalue, 
                    "pval.M" = M.mif.res$Mif$pvalue)
mif.FC      <- list("FC.F"   = F.mif.res$Mif$log2FoldChange,
                    "FC.M"   = M.mif.res$Mif$log2FoldChange)

# Differentially expressed genes in each individual study can also be marked in a matrix DE:
mif.adjpval <- list("adjpval.F" = F.mif.res$Mif$padj,
                    "adjpval.M" = M.mif.res$Mif$padj)

mif.studies      <- c("FEMALES", "MALES")
mif.DE           <- data.frame(mapply(mif.adjpval, FUN=function(x) ifelse(x <= 0.05, 1, 0)))
colnames(mif.DE) <- paste("DE", mif.studies, sep=".")
rownames(mif.DE) <- rownames(F.mif.res$Mif)

# DE returns a list with 1 for genes identified as differentially expressed and 0 otherwise
# Since the proposed p-value combination techniques rely on the assumption that p-values follow a uniform distribution under the null hypothesis, 
# it is necessary to check that the histograms of raw-pvalues reflect that assumption
# The peak near 0 corresponds to differentially expressed genes, no other peak should appear
pdf(paste0(Sys.Date(), "_rawp_histogram_mifepristone.pdf"), height = 5, width = 8)
par(mfrow = c(1,2))
hist(mif.rawpval[[1]], breaks=100, col="grey", main="Female mifepristone", xlab="Raw p-value")
hist(mif.rawpval[[2]], breaks=100, col="grey", main="Male mifepristone", xlab="Raw p-value")
par(mfrow = c(1,1))
dev.off()

# Account for DESeq2 automatic independent filtering. 
mif.filtered <- lapply(mif.adjpval, FUN=function(pval) which(is.na(pval)))
mif.rawpval[[1]][mif.filtered[[1]]]=NA
mif.rawpval[[2]][mif.filtered[[2]]]=NA

# Check raw p-values are roughly uniformly distributed
pdf(paste0(Sys.Date(), "_rawp_histogram_mifepristone_FILTERED.pdf"), height = 5, width = 8)
par(mfrow = c(1,2))
hist(mif.rawpval[[1]], breaks=100, col="grey", main="Female mifepristone", xlab="Raw p-value")
hist(mif.rawpval[[2]], breaks=100, col="grey", main="Male mifepristone", xlab="Raw p-value")
par(mfrow = c(1,1))
dev.off()

# Use of p-value combination techniques
# The p-value combination using the Fisher method may be performed with the fishercomb
# function, and the subsequent p-values obtained from the meta-analysis may be examined (Figure 3, left):
mif.fishcomb    <- fishercomb(mif.rawpval, BHth = 0.05)
mif.invnormcomb <- invnorm(mif.rawpval, nrep=c(5,6), BHth = 0.05) # F: 3OC, 2OM; M: 3OC, 3OM

pdf(paste0(Sys.Date(), "_META_pvalue_histogram_mifepristone.pdf"), height = 5, width = 8)
par(mfrow = c(1,2))
hist(mif.fishcomb$rawpval   , breaks=100, col="grey", main="mifepristone - Fisher method"        , xlab = "Raw p-values (meta-analysis)")
hist(mif.invnormcomb$rawpval, breaks=100, col="grey", main="mifepristone - Inverse normal method", xlab = "Raw p-values (meta-analysis)")
par(mfrow = c(1,1))
dev.off()

# Finally, we suggest summarizing the results of the individual differential analyses as well
# as the differential meta-analysis (using the Fisher and inverse normal methods) in a data.frame:
mif.DEresults <- data.frame(mif.DE,
                            "DE.fishercomb"= ifelse(mif.fishcomb$adjpval<=0.05,1,0),
                            "DE.invnorm"   = ifelse(mif.invnormcomb$adjpval<=0.05,1,0))
head(mif.DEresults)
#              DE.FEMALES DE.MALES DE.fishercomb DE.invnorm
# LOC107382895          0        0             0          0
# LOC107382813          0        0             0          0
# fmr1                  0        0             0          0
# aff2                  0        0             0          0
# ids                   0        0             0          0
# dcps                  0        0             0          0

# Treatment of conflicts in differential expression
# We build a matrix signsFC gathering all signs of fold changes from individual studies.
mif.signsFC      <- mapply(mif.FC, FUN=function(x) sign(x))
mif.sumsigns     <- apply(mif.signsFC,1,sum)
mif.commonsgnFC  <- ifelse(abs(mif.sumsigns)==dim(mif.signsFC)[2], sign(mif.sumsigns),0)

# The vector commonsgnFC returns:
#    - a value of 1 if the gene has a positive log2 fold change in all studies, 
#    - -1 if the gene has a negative log2 fold change in all studies,
#    - 0 if contradictory log2 fold changes are observed across studies (i.e., positive in one and negative in the other). 
# By examining the elements of commonsgnFC, it is thus possible to identify genes displaying contradictory differential expression among studies.
mif.unionDE     <- unique(c(mif.fishcomb$DEindices,mif.invnormcomb$DEindices))
mif.FC.selecDE  <- data.frame(mif.DEresults[mif.unionDE,],
                              do.call(cbind,mif.FC)[mif.unionDE],
                              signFC=mif.commonsgnFC[mif.unionDE])
mif.keepDE     <- mif.FC.selecDE[which(abs(mif.FC.selecDE$signFC)==1),]
mif.conflictDE <- mif.FC.selecDE[which(mif.FC.selecDE$signFC == 0),]
dim(mif.FC.selecDE) # 841    6
dim(mif.keepDE    ) # 660    6
dim(mif.conflictDE) # 181    6

# IDD, IRR and Venn Diagram
mif.fishcomb_de <- rownames(mif.keepDE)[which(mif.keepDE[,"DE.fishercomb"]==1)]
mif.invnorm_de  <- rownames(mif.keepDE)[which(mif.keepDE[,"DE.invnorm"]==1)]
mif.indstudy_de <- list(rownames(mif.keepDE)[which(mif.keepDE[,"DE.FEMALES"]==1)], rownames(mif.keepDE)[which(mif.keepDE[,"DE.MALES"]==1)])
IDD.IRR(mif.fishcomb_de,mif.indstudy_de)
# DE     IDD    Loss     IDR     IRR 
# 623.00 236.00   0.00  37.88   0.00 
IDD.IRR(mif.invnorm_de ,mif.indstudy_de)
# DE     IDD    Loss     IDR     IRR 
# 550.00 259.00  96.00  47.09  24.81 

length(intersect(mif.fishcomb_de,mif.invnorm_de)) # 513

mif.meta.results <- cbind(mif.DEresults,
                          "meta_sign"   = mif.commonsgnFC,
                          "Fisher_FDR"  = mif.fishcomb$adjpval   ,
                          "InvNorm_FDR" = mif.invnormcomb$adjpval) 
# sum(mif.meta.results$InvNorm_FDR[abs(mif.meta.results$meta_sign) == 1] < 0.05) # 550

save(mif.meta.results, file = paste0(Sys.Date(),"_MIFEPRISTONE_meta_analysis_results.RData"))
#####################################################################################################################

#####################################################################################################################
### export to excel
options(java.parameters = "-Xmx16g" )
require(openxlsx)


write.xlsx(list("AGING" = age.meta.results, "MIF" = mif.meta.results), rowNames = TRUE, file = paste0(Sys.Date(),"_GRZ_Mifepristone_Aging_DESeq2_Results_METARNASEQ.xlsx"))
#####################################################################################################################

#####################################################################################################################
#### 4. Comparison
# mif.meta.results
# age.meta.results

####### a. Aging up/Mif down
AgeU.MifD <- list("Meta Aging Up" = rownames(age.meta.results)[bitAnd(age.meta.results$meta_sign > 0, age.meta.results$InvNorm_FDR < 0.05)>0],
                  "Meta Mif Down" = rownames(mif.meta.results)[bitAnd(mif.meta.results$meta_sign < 0, mif.meta.results$InvNorm_FDR < 0.05)>0])
lapply(AgeU.MifD,length)
#$`Meta Aging Up`
# [1] 2187
# 
# $`Meta Mif Down`
# [1] 276

my.Venn <- Venn(AgeU.MifD)

pdf(paste0(Sys.Date(),"_META_Aging_Up_MIf_Down.pdf"))
plot(my.Venn, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
dev.off()

test.AgeU.MifD <- fisher.test(matrix(c(134, 2053, 142, nrow(age.meta.results) - 134 - 2053 - 142),2,2))
test.AgeU.MifD$p.value
# 9.10024e-51

AUMD_g <- intersect(AgeU.MifD$`Meta Aging Up`, AgeU.MifD$`Meta Mif Down`)

####### b. Aging down/Mif up
AgeD.MifU <- list("Meta Aging Down" = rownames(age.meta.results)[bitAnd(age.meta.results$meta_sign < 0, age.meta.results$InvNorm_FDR < 0.05)>0],
                  "Meta Mif Up"     = rownames(mif.meta.results)[bitAnd(mif.meta.results$meta_sign > 0, mif.meta.results$InvNorm_FDR < 0.05)>0])
lapply(AgeD.MifU,length)
# $`Meta Aging Down`
# [1] 2771
# 
# $`Meta Mif Up`
# [1] 274

my.Venn <- Venn(AgeD.MifU)

pdf(paste0(Sys.Date(),"_META_Aging_Down_MIf_Up.pdf"))
plot(my.Venn, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
dev.off()

test.AgeD.MifU <- fisher.test(matrix(c(158, 2613, 116, nrow(age.meta.results) - 158 - 2613 - 116),2,2))
test.AgeD.MifU$p.value
# 5.468514e-59

ADMU_g <- intersect(AgeD.MifU$`Meta Aging Down`, AgeD.MifU$`Meta Mif Up`)

write.table(AgeD.MifU$`Meta Aging Down`, file = paste0(Sys.Date(),"_META_Aging_Down_FDR5.txt"), quote = F, row.names = F, col.names = F)
write.table(AgeD.MifU$`Meta Mif Up`    , file = paste0(Sys.Date(),"_META_Mif_Up_FDR5.txt")    , quote = F, row.names = F, col.names = F)

write.table(AgeU.MifD$`Meta Aging Up`  , file = paste0(Sys.Date(),"_META_Aging_Up_FDR5.txt")  , quote = F, row.names = F, col.names = F)
write.table(AgeU.MifD$`Meta Mif Down`  , file = paste0(Sys.Date(),"_META_Mif_Down_FDR5.txt")  , quote = F, row.names = F, col.names = F)

write.table(ADMU_g , file = paste0(Sys.Date(),"_META_Aging_Down_Mif_up_FDR5.txt")  , quote = F, row.names = F, col.names = F)
write.table(AUMD_g , file = paste0(Sys.Date(),"_META_Aging_Up_Mif_down_FDR5.txt")  , quote = F, row.names = F, col.names = F)

write.table(rownames(age.meta.results) , file = paste0(Sys.Date(),"_META_BACKGROUND.txt")  , quote = F, row.names = F, col.names = F)
################################################################################################################################################


#####################################################################################################################
#### 5. ORA analysis

# load processed genesets
load('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways_GeneSets/ENSEMBL/2025-03-18_Killified_GeneSetCollections_for_ClusterProfilerORA.RData')
# "Killi.ENS.GO_ALL" "Killi.ENS.GO_BP"  "Killi.ENS.GO_CC"  "Killi.ENS.GO_MF" 

load('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways_GeneSets/MsigDB/2024-02-23_Killified_GeneSetCollections_for_ClusterProfilerORA.RData')
# "Killi.m2.reactome" 

########################################################################
ADMU_g.ora.go <- enricher(gene             = ADMU_g,
                          pAdjustMethod    = "BH",
                          universe         = rownames(age.meta.results),
                          minGSSize        = 5,
                          maxGSSize        = 5000,
                          TERM2GENE        = Killi.ENS.GO_ALL)

ADMU_g.ora.re <- enricher(gene             = ADMU_g,
                          pAdjustMethod    = "BH",
                          universe         = rownames(age.meta.results),
                          minGSSize        = 5,
                          maxGSSize        = 5000,
                          TERM2GENE        = Killi.m2.reactome)

ADMU_g.res.go <- ADMU_g.ora.go@result
ADMU_g.res.go <- ADMU_g.res.go[ADMU_g.res.go$p.adjust < 0.05,]

ADMU_g.res.re <- ADMU_g.ora.re@result
ADMU_g.res.re <- ADMU_g.res.re[ADMU_g.res.re$p.adjust < 0.05,]



########################################################################
AUMD_g.ora.go <- enricher(gene             = AUMD_g,
                          pAdjustMethod    = "BH",
                          universe         = rownames(age.meta.results),
                          minGSSize        = 5,
                          maxGSSize        = 5000,
                          TERM2GENE        = Killi.ENS.GO_ALL)

AUMD_g.ora.re <- enricher(gene             = AUMD_g,
                          pAdjustMethod    = "BH",
                          universe         = rownames(age.meta.results),
                          minGSSize        = 5,
                          maxGSSize        = 5000,
                          TERM2GENE        = Killi.m2.reactome)

AUMD_g.res.go <- AUMD_g.ora.go@result
AUMD_g.res.go <- AUMD_g.res.go[AUMD_g.res.go$p.adjust < 0.05,]

AUMD_g.res.re <- AUMD_g.ora.re@result
AUMD_g.res.re <- AUMD_g.res.re[AUMD_g.res.re$p.adjust < 0.05,]


write.xlsx(list("AGE_UP_MIF_DOWN" = rbind(AUMD_g.res.go, AUMD_g.res.re),
                "AGE_DOWN_MIF_UP" = rbind(ADMU_g.res.go, ADMU_g.res.re)), 
           rowNames = TRUE, file = paste0(Sys.Date(),"_GRZ_Mifepristone_Aging_FUNCTIONAL_ENRICHMENT_Results_METARNASEQ.xlsx"))


pdf(paste0(Sys.Date(), "_AUMD_Top5_GO.pdf"), height = 4, width = 6)
dotplot(AUMD_g.ora.go, showCategory = 5)
dev.off()

pdf(paste0(Sys.Date(), "_AUMD_Top5_REACTOME.pdf"), height = 4, width = 6)
dotplot(AUMD_g.ora.re, showCategory = 5)
dev.off()

pdf(paste0(Sys.Date(), "_ADMU_Top5_GO.pdf"), height = 4, width = 6)
dotplot(ADMU_g.ora.go, showCategory = 5)
dev.off()

pdf(paste0(Sys.Date(), "_ADMU_Top5_REACTOME.pdf"), height = 4, width = 6)
dotplot(ADMU_g.ora.re, showCategory = 5)
dev.off()

#####################################################################################################################

# filter for plotting
flt.AUMD_g.res.go <- data.frame(AUMD_g.res.go[1:5,]) # Top 5
flt.AUMD_g.res.re <- data.frame(AUMD_g.res.re[1:5,]) # Top 5
flt.ADMU_g.res.go <- data.frame(ADMU_g.res.go[1:5,]) # Top 5
flt.ADMU_g.res.re <- data.frame(ADMU_g.res.re      ) # only 3

flt.AUMD_g.res.go$minusLog10FDR <- -log10(flt.AUMD_g.res.go$p.adjust)
flt.AUMD_g.res.re$minusLog10FDR <- -log10(flt.AUMD_g.res.re$p.adjust)
flt.ADMU_g.res.go$minusLog10FDR <- -log10(flt.ADMU_g.res.go$p.adjust)
flt.ADMU_g.res.re$minusLog10FDR <- -log10(flt.ADMU_g.res.re$p.adjust)

get_enrich <- function (cp.res) {
  cp.res$fg.ratio <- as.numeric(unlist(lapply(strsplit(cp.res$GeneRatio,"/"),'[[',1)))/as.numeric(unlist(lapply(strsplit(cp.res$GeneRatio,"/"),'[[',2)))
  cp.res$bg.ratio <- as.numeric(unlist(lapply(strsplit(cp.res$BgRatio,"/"),'[[',1)))/as.numeric(unlist(lapply(strsplit(cp.res$BgRatio,"/"),'[[',2)))
  cp.res$Enrich   <-  cp.res$fg.ratio/cp.res$bg.ratio
  return(cp.res)
}

flt.AUMD_g.res.go <- get_enrich(flt.AUMD_g.res.go)
flt.AUMD_g.res.re <- get_enrich(flt.AUMD_g.res.re)
flt.ADMU_g.res.go <- get_enrich(flt.ADMU_g.res.go)
flt.ADMU_g.res.re <- get_enrich(flt.ADMU_g.res.re)

flt.AUMD_g.res.go$condition <- "AUMD_g"
flt.AUMD_g.res.re$condition <- "AUMD_g"
flt.ADMU_g.res.go$condition <- "ADMU_g"
flt.ADMU_g.res.re$condition <- "ADMU_g"

# get merged datafame for ggplot
my.max <- 30
my.min <- 0
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
flt.AUMD_g.res.go$Description <- factor(flt.AUMD_g.res.go$Description, levels = rev(unique(flt.AUMD_g.res.go$Description)))
flt.AUMD_g.res.re$Description <- factor(flt.AUMD_g.res.re$Description, levels = rev(unique(flt.AUMD_g.res.re$Description)))
flt.ADMU_g.res.go$Description <- factor(flt.ADMU_g.res.go$Description, levels = rev(unique(flt.ADMU_g.res.go$Description)))
flt.ADMU_g.res.re$Description <- factor(flt.ADMU_g.res.re$Description, levels = rev(unique(flt.ADMU_g.res.re$Description)))

pdf(paste0(Sys.Date(),"_META_Age_UP_Mif_Down_Top5_GO.pdf"),height = 4, width=10)
my.plot.1 <- ggplot(flt.AUMD_g.res.go,aes(x=condition,y=Description,colour=Enrich,size=minusLog10FDR))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot.1 <- my.plot.1 + ggtitle("AUMD ORA") + labs(x = "", y = "GeneSet")
my.plot.1 <- my.plot.1 + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(my.min,my.max))
my.plot.1 <- my.plot.1 + scale_size_area(limits = c(1.5,15)) 
print(my.plot.1)
dev.off()

pdf(paste0(Sys.Date(),"_META_Age_UP_Mif_Down_Top5_REACTOME.pdf"), height = 4, width=15)
my.plot.2 <- ggplot(flt.AUMD_g.res.re,aes(x=condition,y=Description,colour=Enrich,size=minusLog10FDR))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot.2 <- my.plot.2 + ggtitle("AUMD ORA") + labs(x = "", y = "GeneSet")
my.plot.2 <- my.plot.2 + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(my.min,my.max))
my.plot.2 <- my.plot.2 + scale_size_area(limits = c(1.5,15)) 
print(my.plot.2)
dev.off()


pdf(paste0(Sys.Date(),"_META_Age_Down_Mif_Up_Top5_GO.pdf"),height = 4, width=10)
my.plot.3 <- ggplot(flt.ADMU_g.res.go,aes(x=condition,y=Description,colour=Enrich,size=minusLog10FDR))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot.3 <- my.plot.3 + ggtitle("ADMU ORA") + labs(x = "", y = "GeneSet")
my.plot.3 <- my.plot.3 + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(my.min,my.max))
my.plot.3 <- my.plot.3 + scale_size_area(limits = c(1.5,15)) 
print(my.plot.3)
dev.off()

pdf(paste0(Sys.Date(),"_META_Age_Down_Mif_Up_Top5_REACTOME.pdf"), height = 4, width=15)
my.plot.4 <- ggplot(flt.ADMU_g.res.re,aes(x=condition,y=Description,colour=Enrich,size=minusLog10FDR))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot.4 <- my.plot.4 + ggtitle("ADMU ORA") + labs(x = "", y = "GeneSet")
my.plot.4 <- my.plot.4 + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(my.min,my.max))
my.plot.4 <- my.plot.4 + scale_size_area(limits = c(1.5,15)) 
print(my.plot.4)
dev.off()


pdf(paste0(Sys.Date(),"_META_Age_Mif_Rescue_Top5_GO_REACTOME.pdf"), width = 25, height = 7)
gridExtra::grid.arrange(my.plot.1, my.plot.2, my.plot.3, my.plot.4, nrow = 2)
dev.off()

#####################################################################################################################

#######################
sink(file = paste(Sys.Date(),"_AGE_MIF_META_RNAseq_BrainAtlas_session_Info.txt", sep =""))
sessionInfo()
sink()



