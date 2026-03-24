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

# 2025-06-27
# make jitter plot for bulk ATAC dataset

################################################################################################
load("/Volumes/BB_HQ_3/Brain_aging_Simons/ATAC/Bulk_ATAC_Analysis_Folder/DESeq2/2024-04-09_ZMZ_brain_Aging_ATAC_BOTH.RData")
load("/Volumes/BB_HQ_3/Brain_aging_Simons/ATAC/Bulk_ATAC_Analysis_Folder/DESeq2/2024-04-09_GRZ_brain_Aging_ATAC_BOTH.RData")

deseq2.bATAC.res <- list("Brain_GRZ" = res.age.grz,
                         "Brain_ZMZ" = res.age.zmz)


######### Make jitter plot of DA peaks #########
## Order by pvalue:
age.results <- lapply(deseq2.bATAC.res,function(x) {x[order(x$padj),]})
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

pdf(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_DESeq2_with_reg_colors_FDR5_GENES_ONLY.pdf"), width = 3, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 2.5),
     ylim = c(-0.25, 0.25),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change per week of life"
)
abline(h = 0)
abline(h = seq(-0.5, 0.5, by = 0.1)[-6],
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
     at = 1:2,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-0.5, 0.5, by = 0.1))
box()
dev.off()

png(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_DESeq2_with_reg_colors_FDR5_GENES_ONLY.png"), width = 1000, height = 2200, res = 300, units = "px")
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 2.5),
     ylim = c(-0.25, 0.25),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change per week of life"
)
abline(h = 0)
abline(h = seq(-0.5, 0.5, by = 0.1)[-6],
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
     at = 1:2,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-0.5, 0.5, by = 0.1))
box()
dev.off()
###############################################



################################################################################################

#######################
sink(file = paste0(Sys.Date(),"_Killi_Brain_aging_bulkATAC_session_Info.txt"))
sessionInfo()
sink()
