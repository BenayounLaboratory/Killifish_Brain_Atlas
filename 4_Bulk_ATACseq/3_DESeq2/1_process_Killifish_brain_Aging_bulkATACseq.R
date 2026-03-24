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

# 2024-04-09
# Killifish brain bulk ATAC aging
# use Diffbind count on MSPC consensus peaks from MACS2

################################################################################################
# 0. preprocess Diffbind matrices for use with DEseq2


###############################################
####### Data preparation   ++++   GRZ   #######
###############################################

my.outprefix.grz <- paste(Sys.Date(),"GRZ_brain_Aging_ATAC_DESeq2_Analysis",sep="_")

####### A. Read diffbind normalized counts
my.brain.grz1 <- read.csv('../Diffbind/2024-04-08_GRZ_Aging_Brains_MACS2_MSPC_Normalized_count_matrix.txt', sep = "\t", header = T)
my.brain.grz1$PeakName <- paste(my.brain.grz1$seqnames,my.brain.grz1$start,my.brain.grz1$end,sep = ":")

my.brain.grz <- my.brain.grz1[,c("PeakName",
                                 "YF1" , "YF2", "YF3" , "YF4" , "YF5", "YF6" ,   
                                 "OF1" , "OF2", "OF3" , "OF4" , "OF5", "OF6" ,  
                                 "YM1" , "YM2", "YM3" , "YM4" , "YM5", "YM6" ,  
                                 "OM1" , "OM2", "OM3" , "OM4" , "OM5", "OM6"  )]

# round counts (DESeq needs integers)
my.brain.grz[,-1]      <- round(my.brain.grz[,-1])
rownames(my.brain.grz) <- my.brain.grz$PeakName

# get the peaks with no reads out
my.good.grz     <- rownames(my.brain.grz)[apply(my.brain.grz[,-1] > 0, 1, sum) >= 12 ] # see deseq2 vignette, need to remove too low peaks/only if in more than half the samples
my.filt.mat.grz <- my.brain.grz[my.good.grz,-1] # 73378 peaks

####### B. get metadata
my.grz.meta <- data.frame("sex"   = c(rep("F",12),rep("M",12)),
                          "age"   = c(rep(6,6),rep(16,6),rep(6,6),rep(16,6)),
                          "batch" = rep(c(rep("c1",3),rep("c2",3)), 4) )
rownames(my.grz.meta) <- colnames(my.brain.grz)[-1]

####### C. Run SVA
# Set null and alternative models (ignore batch)
mod1.grz     = model.matrix(~ sex + age + batch, data = my.grz.meta)
n.sv.be.grz  = num.sv(my.filt.mat.grz, mod1.grz, method="be") # 1

# apply SVAseq algortihm
my.svseq.grz = svaseq(as.matrix(my.filt.mat.grz), mod1.grz, n.sv=n.sv.be.grz, constant = 0.1)

# remove RIN and SV, preserve age and sex
my.clean.grz <- removeBatchEffect(log2(my.filt.mat.grz + 0.1), 
                                  batch      = my.grz.meta$batch, 
                                  covariates = cbind(my.svseq.grz$sv),
                                  design     = mod1.grz[,1:3])

# delog and round data for DEseq2 processing
my.filt.sva.grz <- round(2^my.clean.grz-0.1)

write.table(my.filt.sva.grz, file =  paste0(my.outprefix.grz,"_SVA_counts_matrix.txt") , sep = "\t" , row.names = T, quote = F)


###############################################
####### Data preparation   ++++   ZMZ   #######
###############################################

my.outprefix.zmz <- paste(Sys.Date(),"ZMZ_brain_Aging_ATAC_DESeq2_Analysis",sep="_")

####### A. Read diffbind normalized counts
my.brain.zmz1 <- read.csv('../Diffbind/2024-04-08_ZMZ_Aging_Brains_MACS2_MSPC_Normalized_count_matrix.txt', sep = "\t", header = T)
my.brain.zmz1$PeakName <- paste(my.brain.zmz1$seqnames,my.brain.zmz1$start,my.brain.zmz1$end,sep = ":")

my.brain.zmz <- my.brain.zmz1[,c("PeakName",
                                 "YF1" , "YF2", "YF3" , "YF4" ,   
                                 "OF1" , "OF3" , "OF4",        # OF2 does not exist
                                 "GF1" , "GF2", "GF3" , "GF4" ,  
                                 "YM1" , "YM2", "YM3" , "YM4" ,  
                                 "OM1" , "OM2", "OM3" , "OM4" ,
                                 "GM1" , "GM2", "GM3" , "GM4" )]

# round counts (DESeq needs integers)
my.brain.zmz[,-1]      <- round(my.brain.zmz[,-1])
rownames(my.brain.zmz) <- my.brain.zmz$PeakName

# get the peaks with no reads out
my.good.zmz     <- rownames(my.brain.zmz)[apply(my.brain.zmz[,-1] > 0, 1, sum) >= 12 ] # see deseq2 vignette, need to remove too low peaks/only if in more than half the samples
my.filt.mat.zmz <- my.brain.zmz[my.good.zmz,-1] # 71199 peaks

####### B. get metadata
my.zmz.meta <- data.frame("sex"   = c(rep("F",11),rep("M",12)),
                          "age"   = c(rep(6,4),rep(16,3),rep(26,4),rep(6,4),rep(16,4),rep(26,4)),
                          "batch" = c(rep("c2",2),rep("c3",2),rep("c2",1),rep("c3",2),rep("c2",2),rep("c3",2),rep("c2",2),rep("c3",2),rep("c2",2),rep("c3",2),rep("c2",2),rep("c3",2)))
rownames(my.zmz.meta) <- colnames(my.brain.zmz)[-1]

####### C. Run SVA
# Set null and alternative models (ignore batch)
mod1.zmz     = model.matrix(~ sex + age + batch, data = my.zmz.meta)
n.sv.be.zmz  = num.sv(my.filt.mat.zmz, mod1.zmz, method="be") # 1

# apply SVAseq algortihm
my.svseq.zmz = svaseq(as.matrix(my.filt.mat.zmz), mod1.zmz, n.sv=n.sv.be.zmz, constant = 0.1)

# remove RIN and SV, preserve age and sex
my.clean.zmz <- removeBatchEffect(log2(my.filt.mat.zmz + 0.1), 
                                  batch      = my.zmz.meta$batch, 
                                  covariates = cbind(my.svseq.zmz$sv),
                                  design     = mod1.zmz[,1:3])

# delog and round data for DEseq2 processing
my.filt.sva.zmz <- round(2^my.clean.zmz-0.1)

write.table(my.filt.sva.zmz, file =  paste0(my.outprefix.zmz,"_SVA_counts_matrix.txt") , sep = "\t" , row.names = T, quote = F)


save(my.filt.sva.grz,my.grz.meta,
     my.filt.sva.zmz,my.zmz.meta,
     file = paste0(Sys.Date(),"_MetaData_and_SVA_counts_matrices_KilliBrain_Aging_ATAC.RData") )
################################################################################################


################################################################################################
# 1. DESeq2 on cleaned data

# clean workspace and load clean counts and metadata
load('2024-04-09_MetaData_and_SVA_counts_matrices_KilliBrain_Aging_ATAC.RData')

my.outprefix.zmz <- paste(Sys.Date(),"ZMZ_brain_Aging_ATAC_DESeq2_Analysis",sep="_")
my.outprefix.grz <- paste(Sys.Date(),"GRZ_brain_Aging_ATAC_DESeq2_Analysis",sep="_")

###############################################
#######   DEG analysis   ++++   GRZ     #######
###############################################

# get matrix using age as a modeling covariate
dds.grz <- DESeqDataSetFromMatrix(countData = my.filt.sva.grz,
                                  colData   = my.grz.meta,
                                  design    = ~ age + sex)

# run DESeq normalizations and export results
dds.deseq.grz <- DESeq(dds.grz)

# plot dispersion
my.disp.out.grz <- paste(my.outprefix.grz,"dispersion_plot.pdf",sep="_")

pdf(my.disp.out.grz)
plotDispEsts(dds.deseq.grz)
dev.off()

# normalized expression value
tissue.cts.grz <- getVarianceStabilizedData(dds.deseq.grz)

# color-code
my.colors <- rep("deeppink",24)
my.colors[grep("OF",colnames(my.filt.sva.grz))] <- "deeppink4"
my.colors[grep("YM",colnames(my.filt.sva.grz))] <- "deepskyblue"
my.colors[grep("OM",colnames(my.filt.sva.grz))] <- "deepskyblue4"


# do MDS analysis
mds.result <- cmdscale(1-cor(tissue.cts.grz,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.mds.out <- paste(my.outprefix.grz,"MDS_plot.pdf",sep="_")
pdf(my.mds.out)
plot(x, y,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Brain ATACseq MDS (GRZ)",
     cex=3, pch = 16, col = my.colors,
     xlim = c(-0.1,0.05),
     ylim = c(-0.06,0.05),
     cex.lab = 1.5,
     cex.axis = 1.5,
     las = 1)
dev.off()


my.mds.out <- paste(my.outprefix.grz,"MDS_plot_with_Labels.pdf",sep="_")
pdf(my.mds.out)
plot(x, y,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Brain ATACseq MDS (GRZ)",
     cex=3, pch = 16, col = my.colors,
     xlim = c(-0.1,0.05),
     ylim = c(-0.06,0.05),
     cex.lab = 1.5,
     cex.axis = 1.5,
     las = 1)
text(x, y, colnames(tissue.cts.grz), col = "grey")
dev.off()

# calculate MDS1 distance to median old of the same sex
x.of.av <- median(x[7:12])
x.om.av <- median(x[19:24])
x.yf.av <- median(x[1:6])
x.ym.av <- median(x[13:18])

my.mds.out <- paste(my.outprefix.grz,"MDS_plot_with_Xmedians.pdf",sep="_")
pdf(my.mds.out)
plot(x, y,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Brain ATACseq MDS (GRZ)",
     cex=3, pch = 16, col = my.colors,
     xlim = c(-0.1,0.05),
     ylim = c(-0.06,0.05),
     cex.lab = 1.5,
     cex.axis = 1.5,
     las = 1)
abline(v = x.of.av, col = "deeppink4", lty = "dashed")
abline(v = x.om.av, col = "deepskyblue4", lty = "dashed")
abline(v = x.yf.av, col = "deeppink", lty = "dashed")
abline(v = x.ym.av, col = "deepskyblue", lty = "dashed")
dev.off()

# PCA analysis
my.pos.var <- apply(tissue.cts.grz,1,var) > 0
my.pca <- prcomp(t(tissue.cts.grz[my.pos.var,]),scale = TRUE)
x <- my.pca$x[,1]
y <- my.pca$x[,2]

my.summary <- summary(my.pca)

my.pca.out <- paste(my.outprefix.grz,"PCA_plot.pdf",sep="_")
pdf(my.pca.out)
plot(x,y,
     cex=3, pch = 16, col = my.colors,
     xlab = paste('PC1 (', round(100*my.summary$importance[,1][2],1),"%)", sep=""),
     ylab = paste('PC2 (', round(100*my.summary$importance[,2][2],1),"%)", sep=""),
     cex.lab = 1.5,
     cex.axis = 1.5)
dev.off()

# expression range
pdf(paste0(my.outprefix.grz,"_Normalized_counts_boxplot.pdf"))
boxplot(tissue.cts.grz,col=my.colors,cex=0.5,ylab="Log2 DESeq2 Normalized counts", las = 2, outline = F)
dev.off()

###############################################################################################
## a. model aging with sex as covariate  %%%%%%%%%%%%%%
res.age.grz <- results(dds.deseq.grz, name= "age")

### get the heatmap of aging changes at FDR5; exclude NA
res.age.grz <- res.age.grz[!is.na(res.age.grz$padj),]

genes.aging.grz <- rownames(res.age.grz)[res.age.grz$padj < 0.05]
my.num.aging.grz <- length(genes.aging.grz) # 18884

my.heatmap.out.grz <- paste(my.outprefix.grz,"AGING_Heatmap_FDR5.pdf", sep = "_")
pdf(my.heatmap.out.grz, onefile = F)
my.heatmap.title <- paste("Aging significant (FDR<5%), ", my.num.aging.grz, " peaks",sep="")
pheatmap(tissue.cts.grz[genes.aging.grz,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15)
dev.off()

save(res.age.grz, file = paste(Sys.Date(),"GRZ_brain_Aging_ATAC_BOTH.RData", sep ="_"))

###############################################################################################
## b. sex with age as covariate
res.sex.grz <- results(dds.deseq.grz, contrast = c("sex","F","M")) # FC in females over Males

### get the heatmap of sex dimorphic changes at FDR5; exclude NA
res.sex.grz <- res.sex.grz[!is.na(res.sex.grz$padj),]

genes.sex.grz <- rownames(res.sex.grz)[res.sex.grz$padj < 0.05]
my.num.sex.grz <- length(genes.sex.grz) # 151

my.heatmap.out.grz <- paste(my.outprefix.grz,"SEX_DIM_Heatmap_FDR5.pdf", sep = "_")
pdf(my.heatmap.out.grz, onefile = F)
my.heatmap.title <- paste("Sex significant (FDR<5%), ", my.num.sex.grz, " peaks",sep="")
pheatmap(tissue.cts.grz[genes.sex.grz,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15)
dev.off()

save(res.sex.grz, file = paste(Sys.Date(),"GRZ_brain_Aging_SEX_ATAC.RData", sep ="_"))

# output result tables of combined analysis to text files
my.out.ct.mat <- paste(my.outprefix.grz,"_log2_counts_matrix_DEseq2.txt", sep = "_")
write.table(tissue.cts.grz, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)

my.out.stats.age <- paste(my.outprefix.grz,"AGING_all_genes_statistics.txt", sep = "_")
my.out.stats.sex <- paste(my.outprefix.grz,"SEX_DIM_all_genes_statistics.txt", sep = "_")
write.table(res.age.grz, file = my.out.stats.age , sep = "\t" , row.names = T, quote=F)
write.table(res.sex.grz, file = my.out.stats.sex , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.age <- paste(my.outprefix.grz,"AGING_FDR5_genes_statistics.txt", sep = "_")
my.out.fdr5.sex <- paste(my.outprefix.grz,"SEX_DIM_FDR5_genes_statistics.txt", sep = "_")
write.table(res.age.grz[genes.aging.grz,], file = my.out.fdr5.age, sep = "\t" , row.names = T, quote=F)
write.table(res.sex.grz[genes.sex.grz,], file = my.out.fdr5.sex, sep = "\t" , row.names = T, quote=F)

################################################################################################
# c. annotate Peaks

# get HOMER annotations 
my.peak.annot <- read.csv('../Diffbind/HOMER_2024-04-08_GRZ_Aging_Brains_MACS2_MSPC_diffbind_peaks.xls', sep = "\t", header = T)
colnames(my.peak.annot)[1] <- "PeakID"
my.peak.annot$PeakName <- paste(my.peak.annot$Chr,my.peak.annot$Start-1,my.peak.annot$End,sep = ":")

# clean gene names and genomic annotations
my.peak.annot$Gene.Name <- gsub("gene-","",my.peak.annot$Nearest.PromoterID)
my.peak.annot$Genomic_Context <- unlist(lapply(strsplit(my.peak.annot$Annotation, " "), '[[',1))

unique(my.peak.annot$Genomic_Context)
# [1] "Intergenic"   "intron"       "promoter-TSS" "exon"         NA             "TTS"         

##### Log2 normalized counts
tissue.cts.grz.2 <- data.frame(cbind(rownames(tissue.cts.grz),tissue.cts.grz))
colnames(tissue.cts.grz.2)[1] <- "PeakName"
tissue.cts.grz.annot <- data.frame(merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Genomic_Context", "Distance.to.TSS", "Gene.Name")],tissue.cts.grz.2, by = "PeakName"))
write.table(tissue.cts.grz.annot, file = paste(my.outprefix.grz,"_log2_counts_matrix_DEseq2_PeakAnnot.txt", sep = "_") , sep = "\t" , row.names = T, quote=F)


##### Output SEX bed
res.sex.grz$PeakName <- rownames(res.sex.grz)
res.sex.grz <- data.frame(res.sex.grz)
res.sex.grz.annot    <- data.frame(merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Genomic_Context", "Distance.to.TSS", "Gene.Name" )],res.sex.grz, by = "PeakName"))

my.out.stats.sex.grz <- paste(my.outprefix.grz,"SEX_DIM_all_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.sex.grz.annot, file = my.out.stats.sex.grz , sep = "\t" , row.names = F, quote=F)

my.out.fdr5.sex.grz <- paste(my.outprefix.grz,"SEX_DIM_FDR5_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.sex.grz.annot[res.sex.grz.annot$padj < 0.05,], file = my.out.fdr5.sex.grz, sep = "\t" , row.names = F, quote=F)

my.out.bckgd.sex.grz.BED <- paste0(my.outprefix.grz,"_SEX_DIM_Background.bed")
write.table(res.sex.grz.annot[,c("Chr", "Start", "End","PeakName")], file = my.out.bckgd.sex.grz.BED , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.sex.grz.BED.up <- paste0(my.outprefix.grz,"_SEX_DIM_FDR5_UP.bed")
write.table(res.sex.grz.annot[bitAnd(res.sex.grz.annot$padj < 0.05,res.sex.grz.annot$log2FoldChange >0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.sex.grz.BED.up , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.sex.grz.BED.dwn <- paste0(my.outprefix.grz,"_SEX_DIM_FDR5_DWN.bed")
write.table(res.sex.grz.annot[bitAnd(res.sex.grz.annot$padj < 0.05,res.sex.grz.annot$log2FoldChange <0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.sex.grz.BED.dwn , sep = "\t" , row.names = F, col.names = F, quote=F)


##### Output AGE bed
res.age.grz$PeakName <- rownames(res.age.grz)
res.age.grz        <- data.frame(res.age.grz)
res.age.grz.annot <- merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Genomic_Context",  "Distance.to.TSS", "Gene.Name" )],res.age.grz, by = "PeakName")

my.out.stats.age.grz <- paste(my.outprefix.grz,"AGING_all_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.age.grz.annot, file = my.out.stats.age.grz , sep = "\t" , row.names = F, quote=F)

my.out.fdr5.age.grz <- paste(my.outprefix.grz,"AGING_FDR5_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.age.grz.annot[res.age.grz.annot$padj < 0.05,], file = my.out.fdr5.age.grz, sep = "\t" , row.names = F, quote=F)

my.out.bckgd.age.grz.BED <- paste0(my.outprefix.grz,"_AGING_Background.bed")
write.table(res.age.grz.annot[,c("Chr", "Start", "End","PeakName")], file = my.out.bckgd.age.grz.BED , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.age.grz.BED.up <- paste0(my.outprefix.grz,"_AGING_FDR5_UP_with_Age.bed")
write.table(res.age.grz.annot[bitAnd(res.age.grz.annot$padj < 0.05,res.age.grz.annot$log2FoldChange >0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.age.grz.BED.up , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.age.grz.BED.dwn <- paste0(my.outprefix.grz,"_AGING_FDR5_DWN_with_Age.bed")
write.table(res.age.grz.annot[bitAnd(res.age.grz.annot$padj < 0.05,res.age.grz.annot$log2FoldChange <0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.age.grz.BED.dwn , sep = "\t" , row.names = F, col.names = F, quote=F)


#### Plot some diagnostic plots with genomic location information
tss.dists.grz <- list("All_ATAC_Peaks"      = tissue.cts.grz.annot$Distance.to.TSS,
                      "Male_Biased_Peaks"   = res.sex.grz.annot[bitAnd(res.sex.grz.annot$padj < 0.05, res.sex.grz.annot$log2FoldChange < 0 )>0,]$Distance.to.TSS,
                      "Female_Biased_Peaks" = res.sex.grz.annot[bitAnd(res.sex.grz.annot$padj < 0.05, res.sex.grz.annot$log2FoldChange > 0 )>0,]$Distance.to.TSS,
                      "Age_Up_Peaks"        = res.age.grz.annot[bitAnd(res.age.grz.annot$padj < 0.05, res.age.grz.annot$log2FoldChange > 0 )>0,]$Distance.to.TSS,
                      "Age_Down_Peaks"      = res.age.grz.annot[bitAnd(res.age.grz.annot$padj < 0.05, res.age.grz.annot$log2FoldChange < 0 )>0,]$Distance.to.TSS)

pdf(paste0(my.outprefix.grz,"_violinPlot_distances_to_TSS.pdf"))
vioplot::vioplot(tss.dists.grz[c(1,4,5)], 
                 las = 2, col = c("grey","#CC3333","#333399"),
                 ylab = "Distance to closest TSS (bp)",
                 main = "ATAC seq peak distance to TSS (GRZ)")
abline(h = 0, col = "red", lty = "dashed")
dev.off()

background.grz <- tissue.cts.grz.annot
age.up.grz     <- res.age.grz.annot[bitAnd(res.age.grz.annot$padj < 0.05, res.age.grz.annot$log2FoldChange > 0 )>0,]
age.dwn.grz    <- res.age.grz.annot[bitAnd(res.age.grz.annot$padj < 0.05, res.age.grz.annot$log2FoldChange < 0 )>0,]

cts.background.grz <- aggregate(background.grz$Genomic_Context, by = list(background.grz$Genomic_Context), FUN = length)
cts.age.up.grz     <- aggregate(age.up.grz  $Genomic_Context  , by = list(age.up.grz$Genomic_Context    ), FUN = length)
cts.age.dwn.grz    <- aggregate(age.dwn.grz $Genomic_Context  , by = list(age.dwn.grz$Genomic_Context   ), FUN = length)

pdf(paste0(my.outprefix.grz,"_PieCharts_Genomic_context_FDR5.pdf"), height = 3, width = 9)
par(mfrow=c(1,3) )
pie(as.numeric(cts.background.grz$x), labels = cts.background.grz$Group.1, main = "All ATAC (GRZ)" , col = c("darkorchid1","palevioletred1","darkturquoise","goldenrod1","mediumvioletred") )
pie(as.numeric(cts.age.up.grz$x    ), labels = cts.age.up.grz$Group.1    , main = "ATAC Up (GRZ)"  , col = c("darkorchid1","palevioletred1","darkturquoise","goldenrod1","mediumvioletred")  )
pie(as.numeric(cts.age.dwn.grz$x   ), labels = cts.age.dwn.grz$Group.1   , main = "ATAC Down (GRZ)", col = c("darkorchid1","palevioletred1","darkturquoise","goldenrod1","mediumvioletred")  )
par(mfrow=c(1,1) )
dev.off()

chisq.up <- chisq.test(rbind(cts.age.up.grz$x,cts.background.grz$x))
chisq.up$p.value # [1] 3.241338e-52

chisq.dwn <- chisq.test(rbind(cts.age.dwn.grz$x,cts.background.grz$x))
chisq.dwn$p.value # [1] 1.450899e-169



###############################################
#######   DEG analysis   ++++   ZMZ     #######
###############################################

# get matrix using age as a modeling covariate
dds.zmz <- DESeqDataSetFromMatrix(countData = my.filt.sva.zmz,
                                  colData   = my.zmz.meta,
                                  design    = ~ age + sex)

# run DESeq normalizations and export results
dds.deseq.zmz <- DESeq(dds.zmz)

# plot dispersion
my.disp.out.zmz <- paste(my.outprefix.zmz,"dispersion_plot.pdf",sep="_")

pdf(my.disp.out.zmz)
plotDispEsts(dds.deseq.zmz)
dev.off()

# normalized expression value
tissue.cts.zmz <- getVarianceStabilizedData(dds.deseq.zmz)

# color-code
my.colors <- rep("deeppink",23)
my.colors[grep("OF",colnames(my.filt.sva.zmz))] <- "deeppink4"
my.colors[grep("GF",colnames(my.filt.sva.zmz))] <- "magenta4"
my.colors[grep("YM",colnames(my.filt.sva.zmz))] <- "deepskyblue"
my.colors[grep("OM",colnames(my.filt.sva.zmz))] <- "deepskyblue4"
my.colors[grep("GM",colnames(my.filt.sva.zmz))] <- "royalblue4"


# do MDS analysis
mds.result <- cmdscale(1-cor(tissue.cts.zmz,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.mds.out <- paste(my.outprefix.zmz,"MDS_plot.pdf",sep="_")
pdf(my.mds.out)
plot(x, y,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Brain ATACseq MDS (ZMZ)",
     cex=3, pch = 16, col = my.colors,
     xlim = c(-0.15,0.12),
     ylim = c(-0.07,0.05),
     cex.lab = 1.5,
     cex.axis = 1.5,
     las = 1)
# points(x,y,cex=3)
dev.off()


my.mds.out <- paste(my.outprefix.zmz,"MDS_plot_with_Labels.pdf",sep="_")
pdf(my.mds.out)
plot(x, y,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Brain ATACseq MDS (ZMZ)",
     cex=3, pch = 16, col = my.colors,
     xlim = c(-0.15,0.12),
     ylim = c(-0.07,0.05),
     cex.lab = 1.5,
     cex.axis = 1.5,
     las = 1)
text(x, y, colnames(tissue.cts.zmz), col = "grey")
dev.off()


# PCA analysis
my.pos.var <- apply(tissue.cts.zmz,1,var) > 0
my.pca <- prcomp(t(tissue.cts.zmz[my.pos.var,]),scale = TRUE)
x <- my.pca$x[,1]
y <- my.pca$x[,2]

my.summary <- summary(my.pca)

my.pca.out <- paste(my.outprefix.zmz,"PCA_plot.pdf",sep="_")
pdf(my.pca.out)
plot(x,y,
     cex=3, pch = 16, col = my.colors,
     xlab = paste('PC1 (', round(100*my.summary$importance[,1][2],1),"%)", sep=""),
     ylab = paste('PC2 (', round(100*my.summary$importance[,2][2],1),"%)", sep=""),
     cex.lab = 1.5,
     cex.axis = 1.5)
dev.off()

# expression range
pdf(paste0(my.outprefix.zmz,"_Normalized_counts_boxplot.pdf"))
boxplot(tissue.cts.zmz,col=my.colors,cex=0.5,ylab="Log2 DESeq2 Normalized counts", las = 2, outline = F)
dev.off()

###############################################################################################
## a. model aging with sex as covariate  %%%%%%%%%%%%%%
res.age.zmz <- results(dds.deseq.zmz, name= "age")

### get the heatmap of aging changes at FDR5; exclude NA
res.age.zmz <- res.age.zmz[!is.na(res.age.zmz$padj),]

genes.aging.zmz <- rownames(res.age.zmz)[res.age.zmz$padj < 0.05]
my.num.aging.zmz <- length(genes.aging.zmz) # 15304

my.heatmap.out.zmz <- paste(my.outprefix.zmz,"AGING_Heatmap_FDR5.pdf", sep = "_")
pdf(my.heatmap.out.zmz, onefile = F)
my.heatmap.title <- paste("Aging significant (FDR<5%), ", my.num.aging.zmz, " peaks",sep="")
pheatmap(tissue.cts.zmz[genes.aging.zmz,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15)
dev.off()

save(res.age.zmz, file = paste(Sys.Date(),"ZMZ_brain_Aging_ATAC_BOTH.RData", sep ="_"))

###############################################################################################
## b. sex with age as covariate
res.sex.zmz <- results(dds.deseq.zmz, contrast = c("sex","F","M")) # FC in females over Males

### get the heatmap of sex dimorphic changes at FDR5; exclude NA
res.sex.zmz <- res.sex.zmz[!is.na(res.sex.zmz$padj),]

genes.sex.zmz <- rownames(res.sex.zmz)[res.sex.zmz$padj < 0.05]
my.num.sex.zmz <- length(genes.sex.zmz) # 132

my.heatmap.out.zmz <- paste(my.outprefix.zmz,"SEX_DIM_Heatmap_FDR5.pdf", sep = "_")
pdf(my.heatmap.out.zmz, onefile = F)
my.heatmap.title <- paste("Sex significant (FDR<5%), ", my.num.sex.zmz, " peaks",sep="")
pheatmap(tissue.cts.zmz[genes.sex.zmz,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15)
dev.off()

save(res.sex.zmz, file = paste(Sys.Date(),"ZMZ_brain_Aging_SEX_ATAC.RData", sep ="_"))

# output result tables of combined analysis to text files
my.out.ct.mat <- paste(my.outprefix.zmz,"_log2_counts_matrix_DEseq2.txt", sep = "_")
write.table(tissue.cts.zmz, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)

my.out.stats.age <- paste(my.outprefix.zmz,"AGING_all_genes_statistics.txt", sep = "_")
my.out.stats.sex <- paste(my.outprefix.zmz,"SEX_DIM_all_genes_statistics.txt", sep = "_")
write.table(res.age.zmz, file = my.out.stats.age , sep = "\t" , row.names = T, quote=F)
write.table(res.sex.zmz, file = my.out.stats.sex , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.age <- paste(my.outprefix.zmz,"AGING_FDR5_genes_statistics.txt", sep = "_")
my.out.fdr5.sex <- paste(my.outprefix.zmz,"SEX_DIM_FDR5_genes_statistics.txt", sep = "_")
write.table(res.age.zmz[genes.aging.zmz,], file = my.out.fdr5.age, sep = "\t" , row.names = T, quote=F)
write.table(res.sex.zmz[genes.sex.zmz,], file = my.out.fdr5.sex, sep = "\t" , row.names = T, quote=F)

################################################################################################
# c. annotate Peaks

# get HOMER annotations (Peaks should be the same for GRZ and ZMZ)
my.peak.annot <- read.csv('../Diffbind/HOMER_2024-04-08_ZMZ_Aging_Brains_MACS2_MSPC_diffbind_peaks.xls', sep = "\t", header = T)
colnames(my.peak.annot)[1] <- "PeakID"
my.peak.annot$PeakName <- paste(my.peak.annot$Chr,my.peak.annot$Start-1,my.peak.annot$End,sep = ":")

# clean gene names and genomic annotations
my.peak.annot$Gene.Name <- gsub("gene-","",my.peak.annot$Nearest.PromoterID)
my.peak.annot$Genomic_Context <- unlist(lapply(strsplit(my.peak.annot$Annotation, " "), '[[',1))

unique(my.peak.annot$Genomic_Context)
# "exon"         "intron"       "promoter-TSS" "TTS"          "Intergenic"   NA

##### Log2 normalized counts
tissue.cts.zmz.2 <- data.frame(cbind(rownames(tissue.cts.zmz),tissue.cts.zmz))
colnames(tissue.cts.zmz.2)[1] <- "PeakName"
tissue.cts.zmz.annot <- data.frame(merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Genomic_Context", "Distance.to.TSS", "Gene.Name")],tissue.cts.zmz.2, by = "PeakName"))
write.table(tissue.cts.zmz.annot, file = paste(my.outprefix.zmz,"_log2_counts_matrix_DEseq2_PeakAnnot.txt", sep = "_") , sep = "\t" , row.names = T, quote=F)


##### Output SEX bed
res.sex.zmz$PeakName <- rownames(res.sex.zmz)
res.sex.zmz <- data.frame(res.sex.zmz)
res.sex.zmz.annot    <- data.frame(merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Genomic_Context", "Distance.to.TSS", "Gene.Name" )],res.sex.zmz, by = "PeakName"))

my.out.stats.sex.zmz <- paste(my.outprefix.zmz,"SEX_DIM_all_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.sex.zmz.annot, file = my.out.stats.sex.zmz , sep = "\t" , row.names = F, quote=F)

my.out.fdr5.sex.zmz <- paste(my.outprefix.zmz,"SEX_DIM_FDR5_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.sex.zmz.annot[res.sex.zmz.annot$padj < 0.05,], file = my.out.fdr5.sex.zmz, sep = "\t" , row.names = F, quote=F)

my.out.bckgd.sex.zmz.BED <- paste0(my.outprefix.zmz,"_SEX_DIM_Background.bed")
write.table(res.sex.zmz.annot[,c("Chr", "Start", "End","PeakName")], file = my.out.bckgd.sex.zmz.BED , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.sex.zmz.BED.up <- paste0(my.outprefix.zmz,"_SEX_DIM_FDR5_UP.bed")
write.table(res.sex.zmz.annot[bitAnd(res.sex.zmz.annot$padj < 0.05,res.sex.zmz.annot$log2FoldChange >0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.sex.zmz.BED.up , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.sex.zmz.BED.dwn <- paste0(my.outprefix.zmz,"_SEX_DIM_FDR5_DWN.bed")
write.table(res.sex.zmz.annot[bitAnd(res.sex.zmz.annot$padj < 0.05,res.sex.zmz.annot$log2FoldChange <0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.sex.zmz.BED.dwn , sep = "\t" , row.names = F, col.names = F, quote=F)


##### Output AGE bed
res.age.zmz$PeakName <- rownames(res.age.zmz)
res.age.zmz          <- data.frame(res.age.zmz)
res.age.zmz.annot    <- merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Genomic_Context",  "Distance.to.TSS", "Gene.Name" )],res.age.zmz, by = "PeakName")

my.out.stats.age.zmz <- paste(my.outprefix.zmz,"AGING_all_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.age.zmz.annot, file = my.out.stats.age.zmz , sep = "\t" , row.names = F, quote=F)

my.out.fdr5.age.zmz <- paste(my.outprefix.zmz,"AGING_FDR5_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.age.zmz.annot[res.age.zmz.annot$padj < 0.05,], file = my.out.fdr5.age.zmz, sep = "\t" , row.names = F, quote=F)

my.out.bckgd.age.zmz.BED <- paste0(my.outprefix.zmz,"_AGING_Background.bed")
write.table(res.age.zmz.annot[,c("Chr", "Start", "End","PeakName")], file = my.out.bckgd.age.zmz.BED , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.age.zmz.BED.up <- paste0(my.outprefix.zmz,"_AGING_FDR5_UP_with_Age.bed")
write.table(res.age.zmz.annot[bitAnd(res.age.zmz.annot$padj < 0.05,res.age.zmz.annot$log2FoldChange >0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.age.zmz.BED.up , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.age.zmz.BED.dwn <- paste0(my.outprefix.zmz,"_AGING_FDR5_DWN_with_Age.bed")
write.table(res.age.zmz.annot[bitAnd(res.age.zmz.annot$padj < 0.05,res.age.zmz.annot$log2FoldChange <0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.age.zmz.BED.dwn , sep = "\t" , row.names = F, col.names = F, quote=F)


#### Plot some diagnostic plots with genomic location information
tss.dists.zmz <- list("All_ATAC_Peaks"      = tissue.cts.zmz.annot$Distance.to.TSS,
                      "Male_Biased_Peaks"   = res.sex.zmz.annot[bitAnd(res.sex.zmz.annot$padj < 0.05, res.sex.zmz.annot$log2FoldChange < 0 )>0,]$Distance.to.TSS,
                      "Female_Biased_Peaks" = res.sex.zmz.annot[bitAnd(res.sex.zmz.annot$padj < 0.05, res.sex.zmz.annot$log2FoldChange > 0 )>0,]$Distance.to.TSS,
                      "Age_Up_Peaks"        = res.age.zmz.annot[bitAnd(res.age.zmz.annot$padj < 0.05, res.age.zmz.annot$log2FoldChange > 0 )>0,]$Distance.to.TSS,
                      "Age_Down_Peaks"      = res.age.zmz.annot[bitAnd(res.age.zmz.annot$padj < 0.05, res.age.zmz.annot$log2FoldChange < 0 )>0,]$Distance.to.TSS)

pdf(paste0(my.outprefix.zmz,"_violinPlot_distances_to_TSS.pdf"))
vioplot::vioplot(tss.dists.zmz[c(1,4,5)], 
                 las = 2, col = c("grey","#CC3333","#333399"),
                 ylab = "Distance to closest TSS (bp)",
                 main = "ATAC seq peak distance to TSS (ZMZ)")
abline(h = 0, col = "red", lty = "dashed")
dev.off()

background.zmz <- tissue.cts.zmz.annot
age.up.zmz     <- res.age.zmz.annot[bitAnd(res.age.zmz.annot$padj < 0.05, res.age.zmz.annot$log2FoldChange > 0 )>0,]
age.dwn.zmz    <- res.age.zmz.annot[bitAnd(res.age.zmz.annot$padj < 0.05, res.age.zmz.annot$log2FoldChange < 0 )>0,]

cts.background.zmz <- aggregate(background.zmz$Genomic_Context, by = list(background.zmz$Genomic_Context), FUN = length)
cts.age.up.zmz     <- aggregate(age.up.zmz  $Genomic_Context  , by = list(age.up.zmz$Genomic_Context    ), FUN = length)
cts.age.dwn.zmz    <- aggregate(age.dwn.zmz $Genomic_Context  , by = list(age.dwn.zmz$Genomic_Context   ), FUN = length)

pdf(paste0(my.outprefix.zmz,"_PieCharts_Genomic_context_FDR5.pdf"), height = 3, width = 9)
par(mfrow=c(1,3) )
pie(as.numeric(cts.background.zmz$x), labels = cts.background.zmz$Group.1, main = "All ATAC (ZMZ)" , col = c("darkorchid1","palevioletred1","darkturquoise","goldenrod1","mediumvioletred") )
pie(as.numeric(cts.age.up.zmz$x    ), labels = cts.age.up.zmz$Group.1    , main = "ATAC Up (ZMZ)"  , col = c("darkorchid1","palevioletred1","darkturquoise","goldenrod1","mediumvioletred") )
pie(as.numeric(cts.age.dwn.zmz$x   ), labels = cts.age.dwn.zmz$Group.1   , main = "ATAC Down (ZMZ)", col = c("darkorchid1","palevioletred1","darkturquoise","goldenrod1","mediumvioletred") )
par(mfrow=c(1,1) )
dev.off()

chisq.up <- chisq.test(rbind(cts.age.up.zmz$x,cts.background.zmz$x))
chisq.up$p.value # [1] 4.728198e-22

chisq.dwn <- chisq.test(rbind(cts.age.dwn.zmz$x,cts.background.zmz$x))
chisq.dwn$p.value # [1] 8.694539e-129
################################################################################################


#######################
sink(file = paste0(Sys.Date(),"_Killi_Brain_aging_bulkATAC_session_Info.txt"))
sessionInfo()
sink()
