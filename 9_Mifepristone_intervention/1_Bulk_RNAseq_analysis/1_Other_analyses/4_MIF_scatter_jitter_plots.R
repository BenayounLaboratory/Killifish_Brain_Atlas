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

# 2025-07-17
# scatter/jitters

#####################################################################################################################
#### 1. Run GSEA

load('2025-07-15_Brain_MIF_DEseq2_results.RData')
ls()
# [1] "F.mif.res" "M.mif.res"

age.merge <- merge(F.mif.res$Aging, M.mif.res$Aging, by = "row.names", suffixes = c(".F",".M"))

my.spear.cor <- cor.test(age.merge$log2FoldChange.F, age.merge$log2FoldChange.M, method = 'spearman')
my.rho       <- signif(my.spear.cor$estimate,3)

pdf(paste0(Sys.Date(),"_Bulk_Brain_AGING_F_vs_M_FC_scatterplot.pdf"))
smoothScatter(age.merge$log2FoldChange.F, age.merge$log2FoldChange.M,
              xlim = c(-4,4), ylim = c(-4,4),
              xlab = paste("log2(FC) in Females with aging"),
              ylab = paste("log2(FC) in Males with aging"  ),
              main = "Killifish Brain")
abline(0,1, col = "grey", lty = "dashed")
abline(h = 0, col = "red", lty = "dashed")
abline(v = 0, col = "red", lty = "dashed")
text(-3.9, 4  , paste("Rho = ",my.rho), pos = 4)
text(-3.9, 3.7, paste("p = ",scientific(my.spear.cor$p.value,2)), pos = 4)
dev.off()

#########
mif.merge <- merge(F.mif.res$Mif, M.mif.res$Mif, by = "row.names", suffixes = c(".F",".M"))

my.spear.cor <- cor.test(mif.merge$log2FoldChange.F, mif.merge$log2FoldChange.M, method = 'spearman')
my.rho       <- signif(my.spear.cor$estimate,3)

pdf(paste0(Sys.Date(),"_Bulk_Brain_MIF_F_vs_M_FC_scatterplot.pdf"))
smoothScatter(mif.merge$log2FoldChange.F, mif.merge$log2FoldChange.M,
              xlim = c(-4,4), ylim = c(-4,4),
              xlab = paste("log2(FC) in Females with Mif"),
              ylab = paste("log2(FC) in Males with Mif"  ),
              main = "Killifish Brain")
abline(0,1, col = "grey", lty = "dashed")
abline(h = 0, col = "red", lty = "dashed")
abline(v = 0, col = "red", lty = "dashed")
text(-3.9, 4  , paste("Rho = ",my.rho), pos = 4)
text(-3.9, 3.7, paste("p = ",scientific(my.spear.cor$p.value,2)), pos = 4)
dev.off()
############################################################################################################

############################################################################################################
names(F.mif.res)[1:2] <- c("F_Aging", "F_Mifepristone")
names(M.mif.res)[1:2] <- c("M_Aging", "M_Mifepristone")

bulk.list <- c(F.mif.res[1:2],M.mif.res[1:2])


######### Make jitter plot of DE TEs #########
## Order by pvalue:
age.results <- lapply(bulk.list,function(x) {x[order(x$padj),]})
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

pdf(paste0(Sys.Date(),"_KilliBrain_Aging_MIF_GRZ_stripplot_DESeq2_with_reg_colors_FDR5.pdf"), width = 6, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 4.5),
     ylim = c(-6, 6),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change"
)
abline(h = 0)
abline(h = seq(-6, 6, by = 2)[-4],
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
     at = 1:4,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-6, 6, by = 2))
box()
dev.off()

png(paste0(Sys.Date(),"_KilliBrain_Aging_MIF_GRZ_stripplot_DESeq2_with_reg_colors_FDR5.png"), width = 2400, height = 2200, res = 300, units = "px")
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 4.5),
     ylim = c(-6, 6),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change"
)
abline(h = 0)
abline(h = seq(-6, 6, by = 2)[-4],
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
     at = 1:4,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-6, 6, by = 2))
box()
dev.off()


#######################
sink(file = paste(Sys.Date(),"_MIF_jitter_BrainAtlas_session_Info.txt", sep =""))
sessionInfo()
sink()



