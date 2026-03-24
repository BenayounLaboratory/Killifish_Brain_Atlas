setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Differential_Expression_Age/DGE_Analysis')
options(stringsAsFactors = F)

#### Packages
library('Seurat')         # 
library(sctransform)      # 
library("singleCellTK")   # 

library('muscat')         # 
library('DESeq2')         # 
library('sva')            # 
library('limma')          # 

library(ggplot2)          # 
library(scales)           # 
library("bitops")         # 
library(Vennerable)       # 
library(data.table)       #

library(ComplexHeatmap)   #
library(circlize)         #


theme_set(theme_bw())   

# 2025-06-30
# Plot jitter for TEs only
##########################################################################################################################################

##########################################################################################################################################
### 1. prepare data
load('2024-02-28_pseudobulk_killi_cell_types_AGING_GRZ_DEseq2_objects.RData')
load('2024-02-28_pseudobulk_killi_cell_types_AGING_ZMZ_DEseq2_objects.RData')
deseq.res.list.grz
deseq.res.list.zmz

de.res <- deseq.res.list.grz[[1]]

extract_TE_data <- function(de.res) {
  TE.res <- de.res[grep("^NotFur", rownames(de.res)),]
  TE.res
}

deseq.res.list.TE.grz <- lapply(deseq.res.list.grz, extract_TE_data)
deseq.res.list.TE.zmz <- lapply(deseq.res.list.zmz, extract_TE_data)
##########################################################################################################################################

##########################################################################################################################################
### 2. PLOT


###############################################
#######   DEG analysis   ++++   GRZ     #######
###############################################

######### Make jitter plot of DE TEs #########
## Order by pvalue:
age.results <- lapply(deseq.res.list.TE.grz,function(x) {x[order(x$padj),]})
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

pdf(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_DESeq2_with_reg_colors_FDR5_TEs_ONLY.pdf"), width = 6, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 8.5),
     ylim = c(-0.75, 0.75),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change per week of life"
)
abline(h = 0)
abline(h = seq(-1, 1, by = 0.25)[-5],
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
     at = 1:8,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-1, 1, by = 0.5))
box()
dev.off()

png(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_DESeq2_with_reg_colors_FDR5_TEs_ONLY.png"), width = 2400, height = 2200, res = 300, units = "px")
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 8.5),
     ylim = c(-0.75, 0.75),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change per week of life"
)
abline(h = 0)
abline(h = seq(-1, 1, by = 0.25)[-5],
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
     at = 1:8,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-1, 1, by = 0.5))
box()
dev.off()





###############################################
#######   DEG analysis   ++++   ZMZ     #######
###############################################

######### Make jitter plot of DE TEs #########
## Order by pvalue:
age.results <- lapply(deseq.res.list.TE.zmz,function(x) {x[order(x$padj),]})
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

pdf(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_ZMZ_stripplot_DESeq2_with_reg_colors_FDR5_TEs_ONLY.pdf"), width = 6, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 8.5),
     ylim = c(-0.75, 0.75),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change per week of life"
)
abline(h = 0)
abline(h = seq(-1, 1, by = 0.25)[-5],
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
     at = 1:8,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-1, 1, by = 0.5))
box()
dev.off()

png(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_ZMZ_stripplot_DESeq2_with_reg_colors_FDR5_TEs_ONLY.png"), width = 2400, height = 2200, res = 300, units = "px")
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 8.5),
     ylim = c(-0.75, 0.75),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change per week of life"
)
abline(h = 0)
abline(h = seq(-1, 1, by = 0.25)[-5],
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
     at = 1:8,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-1, 1, by = 0.5))
box()
dev.off()
##########################################################################################################################################


##########################################################################################################################################
### 3. TE Overlap analysis      GRZ/ ZMZ   

########## Compute overlap for each cell types
fisher.res <- data.frame(matrix(0,length(deseq.res.list.TE.grz), 5))
colnames(fisher.res) <- c("Cell_Type","Overlap_up","Overlap_up_enrich_Fisher","Overlap_dwn","Overlap_dwn_enrich_Fisher")
fisher.res$Cell_Type <- names(deseq.res.list.TE.grz)

celltype.common.data <- data.frame(col.names = c("Gene", "Cell_Type",
                                                 "baseMean.grz", "log2FoldChange.grz", "padj.grz",
                                                 "baseMean.zmz", "log2FoldChange.zmz", "padj.zmz"))

for (i in 1:length(deseq.res.list.TE.grz)) {
  
  if (!(i %in% c(2,5))) {  # take care of case where GRZ does not have sig genes (GABAergic neurons/NSPCs)
    
    # gene list
    grz.up   <- rownames(deseq.res.list.TE.grz[[i]])[bitAnd(deseq.res.list.TE.grz[[i]]$padj < 0.05, 
                                                            deseq.res.list.TE.grz[[i]]$log2FoldChange >0)>0]
    grz.dwn  <- rownames(deseq.res.list.TE.grz[[i]])[bitAnd(deseq.res.list.TE.grz[[i]]$padj < 0.05, 
                                                            deseq.res.list.TE.grz[[i]]$log2FoldChange <0)>0]
    zmz.up   <-rownames(deseq.res.list.TE.zmz[[i]])[bitAnd(deseq.res.list.TE.zmz[[i]]$padj < 0.05, 
                                                           deseq.res.list.TE.zmz[[i]]$log2FoldChange >0)>0]
    zmz.dwn  <- rownames(deseq.res.list.TE.zmz[[i]])[bitAnd(deseq.res.list.TE.zmz[[i]]$padj < 0.05, 
                                                            deseq.res.list.TE.zmz[[i]]$log2FoldChange <0)>0]
    genes.bckd  <- union(rownames(deseq.res.list.TE.grz[[i]]),rownames(deseq.res.list.TE.zmz[[i]]))
    
    if( (length(grz.up) > 0) && (length(zmz.up) > 0)) {
      ################ Upregulated genes
      my.criteria.up <- list("GRZ Up"    = grz.up,  
                             "ZMZ Up"    = zmz.up)  
      my.Venn.up <- Venn(my.criteria.up)
      
      
      pdf(paste0(Sys.Date(),"_VennDiagram_",fisher.res$Cell_Type[i],"_UpAging_Overlap_ZMZ_GRZ_TEs.pdf"))
      plot(my.Venn.up, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
      dev.off()
      
      # prepare Fisher hypergeometric test
      a <- length(my.Venn.up@IntersectionSets$`11`)
      b <- length(my.Venn.up@IntersectionSets$`10`)
      c <- length(my.Venn.up@IntersectionSets$`01`)
      d <- length(genes.bckd) - a - b - c
      
      my.fisher.up <- fisher.test(matrix(c(a,b,c,d),2,2))
      
      # put in result table
      fisher.res$Overlap_up[i]                <- a
      fisher.res$Overlap_up_enrich_Fisher[i]  <- my.fisher.up$p.value
      
    }
    
    if( (length(grz.dwn) > 0) && (length(zmz.dwn) > 0)) {
      
      ################ Downregulated genes
      my.criteria.dwn <- list("GRZ down"    = grz.dwn,  
                              "ZMZ down"    = zmz.dwn)  
      my.Venn.dwn <- Venn(my.criteria.dwn)
      
      
      pdf(paste0(Sys.Date(),"_VennDiagram_",fisher.res$Cell_Type[i],"_DownAging_Overlap_ZMZ_GRZ_TEs.pdf"))
      plot(my.Venn.dwn, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
      dev.off()
      
      # prepare Fisher hypergeometric test
      a <- length(my.Venn.dwn@IntersectionSets$`11`)
      b <- length(my.Venn.dwn@IntersectionSets$`10`)
      c <- length(my.Venn.dwn@IntersectionSets$`01`)
      d <- length(genes.bckd) - a - b - c
      
      my.fisher.dwn <- fisher.test(matrix(c(a,b,c,d),2,2))
      
      # put in result table
      fisher.res$Overlap_dwn[i]                <- a
      fisher.res$Overlap_dwn_enrich_Fisher[i]  <- my.fisher.dwn$p.value
    }
    
    ##########################################
    # extract stats/information for overlapping genes
    data.up.grz   <- deseq.res.list.TE.grz[[i]][my.Venn.up@IntersectionSets$`11`,]
    data.up.zmz   <- deseq.res.list.TE.zmz[[i]][my.Venn.up@IntersectionSets$`11`,]
    data.up.merge <- merge(data.up.grz, data.up.zmz, by = "row.names", suffixes = c(".grz",".zmz"))
    colnames(data.up.merge)[1] <- "Gene"
    
    data.dwn.grz   <- deseq.res.list.TE.grz[[i]][my.Venn.dwn@IntersectionSets$`11`,]
    data.dwn.zmz   <- deseq.res.list.TE.zmz[[i]][my.Venn.dwn@IntersectionSets$`11`,]
    data.dwn.merge <- merge(data.dwn.grz, data.dwn.zmz, by = "row.names", suffixes = c(".grz",".zmz"))
    colnames(data.dwn.merge)[1] <- "Gene"
    
    
    my.tmp <- rbind(data.up.merge, data.dwn.merge)
    
    if (nrow(my.tmp) >0) {
      my.tmp$Cell_Type <- fisher.res$Cell_Type[i]
      my.tmp <- my.tmp[,c("Gene","Cell_Type", "baseMean.grz", "log2FoldChange.grz", "padj.grz", "baseMean.zmz", "log2FoldChange.zmz", "padj.zmz")]
      
    }

    if ( i == 1) {
      celltype.common.data <- my.tmp
    } else {
      celltype.common.data <- rbind(celltype.common.data, my.tmp)
    }
    
  }
}

# export significance results
# up overlap is always more than expected by chance, down is more variable
write.table(fisher.res, file = paste0(Sys.Date(),"_Fisher_Aging_FDR5_byCellType_Overlap_ZMZ_GRZ_TEs.txt"), sep = "\t", row.names = F, quote = F)

# add combined fisher FDR for help sorting
celltype.common.data$FisherCombinedPval <- celltype.common.data$padj.grz * celltype.common.data$padj.zmz
write.table(celltype.common.data, file = paste0(Sys.Date(),"_Aging_DEseq2_results_FDR5_byCellType_Overlap_ZMZ_GRZ_TEs.txt"), sep = "\t", row.names = F, quote = F)

table(celltype.common.data$Gene)
# NotFur1-Contig-18-DNA-TcMar NotFur1-GapFilledScaffold-1020-121196-LINE-L1 NotFur1-GapFilledScaffold-3517-21206-LINE-RTE 
# 3                                             1                                             1 
# NotFur1-GapFilledScaffold-375-227342-LINE-R2 NotFur1-GapFilledScaffold-3993-11210-LINE-RTE    NotFur1-GapFilledScaffold-4042-758-LINE-L2 
# 1                                             1                                             1 
# NotFur1-GapFilledScaffold-741-18446-LINE-RTE  NotFur1-GapFilledScaffold-9035-13448-LINE-L1                     NotFur1-LTR-104-LTR-Gypsy 
# 1                                             1                                             1 
# NotFur1-LTR-1133-LTR-Gypsy                       NotFur1-LTR-830-LTR-LTR         NotFur1-rnd1-family1213-LINE-RexBabar 
# 1                                             1                                             1 
# NotFur1-rnd5-family6448-DNA-hAThAT5              NotFur1-rnd6-family1840-LTR-ERV1                 NotFur1-Unknown-627-LINE-LINE 
# 1                                             1                                             1 
# NotFur1-Unknown-750-LINE-LINE 
# 1 


# #######################
# sink(file = paste(Sys.Date(),"_MuscatDEseq2_PB_DESeq2_GSEA_scRNAseq_BrainAtlas_session_Info.txt", sep =""))
# sessionInfo()
# sink()