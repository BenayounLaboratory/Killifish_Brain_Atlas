setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/TF_activity_analysis')
options(stringsAsFactors = F)

library("DESeq2")        #
library("decoupleR")     #
library("OmnipathR")     #

library(dplyr)
library(tibble)
library(ggplot2)
library(beeswarm)
library(ggrepel)

library(bitops)
library(ComplexHeatmap)

library(Vennerable)       # 

theme_set(theme_bw())   

# 2025-07-07
# make jitter

#####################################################################################################################
#### 0. Load necessary objects
load('2025-07-07_DecoupleR_results.RData'   )

ls()
# [1] "decoupler.list.grz" "decoupler.list.zmz"

head(decoupler.list.grz$Astrocytes_Radial_Glia)
# statistic  source       condition  score p_value
# <chr>      <chr>        <chr>      <dbl>   <dbl>
# 1 norm_fgsea LOC107372367 stat       1.45  0.0593 
# 2 norm_fgsea LOC107372434 stat      -1.60  0.00960
# 3 norm_fgsea LOC107372679 stat      -0.984 0.485  
# 4 norm_fgsea LOC107373090 stat      -0.779 0.593  
# 5 norm_fgsea LOC107373153 stat       1.09  0.304  
# 6 norm_fgsea LOC107373230 stat      -0.755 0.748  

lapply(decoupler.list.grz, summary)
#####################################################################################################################


##########################################################################################################################################
### 2. PLOT


###############################################
#######   DecoupleR analysis   ++++   GRZ     #######
###############################################

######### Make jitter plot of DE TEs #########
## Order by pvalue:
age.results <- lapply(decoupler.list.grz,function(x) {x[order(x$p_value),]})
n        <- sapply(age.results, nrow)
names(n) <- names(age.results)

cols <- list()
xlab <- character(length = length(age.results))
for(i in seq(along = age.results)){
  cols[[i]]      <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), n[i]) # grey60
  ind.sig.i      <- age.results[[i]]$p_value < 0.05
  ind.sig.i.up   <- bitAnd(age.results[[i]]$p_value < 0.05, age.results[[i]]$score >0)>0
  ind.sig.i.down <- bitAnd(age.results[[i]]$p_value < 0.05, age.results[[i]]$score <0)>0
  
  cols[[i]][ind.sig.i.up]   <- "#CC3333"
  cols[[i]][ind.sig.i.down] <- "#333399"
  xlab[i] <- paste(names(age.results)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}
names(cols) <- names(age.results)

pdf(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_DecoupleR_with_reg_colors_FDR5.pdf"), width = 6, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 8.5),
     ylim = c(-3, 3),
     axes = FALSE,
     xlab = "",
     ylab = "DecoupleR fgsea score"
)
abline(h = 0)
abline(h = seq(-3, 3, by = 1)[-4],
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
         y = rev(age.results[[i]]$score),
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
     at = seq(-3, 3, by = 1))
box()
dev.off()

png(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_DecoupleR_with_reg_colors_FDR5.png"), width = 2400, height = 2200, res = 300, units = "px")
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 8.5),
     ylim = c(-3, 3),
     axes = FALSE,
     xlab = "",
     ylab = "DecoupleR fgsea score"
)
abline(h = 0)
abline(h = seq(-3, 3, by = 1)[-4],
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
         y = rev(age.results[[i]]$score),
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
     at = seq(-3, 3, by = 1))
box()
dev.off()





###############################################
#######   DecoupleR analysis   ++++   ZMZ     #######
###############################################

######### Make jitter plot of DE TEs #########
## Order by pvalue:
age.results <- lapply(decoupler.list.zmz,function(x) {x[order(x$p_value),]})
n        <- sapply(age.results, nrow)
names(n) <- names(age.results)

cols <- list()
xlab <- character(length = length(age.results))
for(i in seq(along = age.results)){
  cols[[i]]      <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), n[i]) # grey60
  ind.sig.i      <- age.results[[i]]$p_value < 0.05
  ind.sig.i.up   <- bitAnd(age.results[[i]]$p_value < 0.05, age.results[[i]]$score >0)>0
  ind.sig.i.down <- bitAnd(age.results[[i]]$p_value < 0.05, age.results[[i]]$score <0)>0
  
  cols[[i]][ind.sig.i.up]   <- "#CC3333"
  cols[[i]][ind.sig.i.down] <- "#333399"
  xlab[i] <- paste(names(age.results)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}
names(cols) <- names(age.results)

pdf(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_ZMZ_stripplot_DecoupleR_with_reg_colors_FDR5.pdf"), width = 6, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 8.5),
     ylim = c(-3, 3),
     axes = FALSE,
     xlab = "",
     ylab = "DecoupleR fgsea score"
)
abline(h = 0)
abline(h = seq(-3, 3, by = 1)[-4],
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
         y = rev(age.results[[i]]$score),
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
     at = seq(-3, 3, by = 1))
box()
dev.off()

png(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_ZMZ_stripplot_DecoupleR_with_reg_colors_FDR5.png"), width = 2400, height = 2200, res = 300, units = "px")
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 8.5),
     ylim = c(-3, 3),
     axes = FALSE,
     xlab = "",
     ylab = "DecoupleR fgsea score"
)
abline(h = 0)
abline(h = seq(-3, 3, by = 1)[-4],
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
         y = rev(age.results[[i]]$score),
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
     at = seq(-3, 3, by = 1))
box()
dev.off()
##########################################################################################################################################


##########################################################################################################################################
### 3. REgulon Overlap analysis      GRZ/ ZMZ   

########## Compute overlap for each cell types
fisher.res <- data.frame(matrix(0,length(decoupler.list.grz), 5))
colnames(fisher.res) <- c("Cell_Type","Overlap_up","Overlap_up_enrich_Fisher","Overlap_dwn","Overlap_dwn_enrich_Fisher")
fisher.res$Cell_Type <- names(decoupler.list.grz)

celltype.common.data <- data.frame(col.names = c("Gene", "Cell_Type",
                                                 "score.grz", "p_value.grz",
                                                 "score.zmz", "p_value.zmz"))

for (i in 1:length(decoupler.list.grz)) {
  
  # regulon list
  grz.up   <- rownames(decoupler.list.grz[[i]])[bitAnd(decoupler.list.grz[[i]]$p_value < 0.05, 
                                                       decoupler.list.grz[[i]]$score >0)>0]
  grz.dwn  <- rownames(decoupler.list.grz[[i]])[bitAnd(decoupler.list.grz[[i]]$p_value < 0.05, 
                                                       decoupler.list.grz[[i]]$score <0)>0]
  zmz.up   <-rownames(decoupler.list.zmz[[i]])[bitAnd(decoupler.list.zmz[[i]]$p_value < 0.05, 
                                                      decoupler.list.zmz[[i]]$score >0)>0]
  zmz.dwn  <- rownames(decoupler.list.zmz[[i]])[bitAnd(decoupler.list.zmz[[i]]$p_value < 0.05, 
                                                       decoupler.list.zmz[[i]]$score <0)>0]
  reg.bckd  <- union(rownames(decoupler.list.grz[[i]]),rownames(decoupler.list.zmz[[i]]))
  
  if( (length(grz.up) > 0) && (length(zmz.up) > 0)) {
    ################ Upregulated genes
    my.criteria.up <- list("GRZ Up"    = grz.up,  
                           "ZMZ Up"    = zmz.up)  
    my.Venn.up <- Venn(my.criteria.up)
    
    
    pdf(paste0(Sys.Date(),"_VennDiagram_",fisher.res$Cell_Type[i],"_UpAging_Overlap_ZMZ_GRZ_Regulons.pdf"))
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
    
    
    pdf(paste0(Sys.Date(),"_VennDiagram_",fisher.res$Cell_Type[i],"_DownAging_Overlap_ZMZ_GRZ_Regulons.pdf"))
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
  data.up.grz   <- decoupler.list.grz[[i]][my.Venn.up@IntersectionSets$`11`,]
  data.up.zmz   <- decoupler.list.zmz[[i]][my.Venn.up@IntersectionSets$`11`,]
  data.up.merge <- merge(data.up.grz, data.up.zmz, by = "source", suffixes = c(".grz",".zmz"))
  colnames(data.up.merge)[1] <- "TF_regulon"
  
  data.dwn.grz   <- decoupler.list.grz[[i]][my.Venn.dwn@IntersectionSets$`11`,]
  data.dwn.zmz   <- decoupler.list.zmz[[i]][my.Venn.dwn@IntersectionSets$`11`,]
  data.dwn.merge <- merge(data.dwn.grz, data.dwn.zmz, by = "source", suffixes = c(".grz",".zmz"))
  colnames(data.dwn.merge)[1] <- "TF_regulon"
  
  
  my.tmp <- rbind(data.up.merge, data.dwn.merge)
  
  if (nrow(my.tmp) >0) {
    my.tmp$Cell_Type <- fisher.res$Cell_Type[i]
    my.tmp <- my.tmp[,c("TF_regulon","Cell_Type", "score.grz", "p_value.grz", "score.zmz", "p_value.zmz")]
    
  }
  
  if ( i == 1) {
    celltype.common.data <- my.tmp
  } else {
    celltype.common.data <- rbind(celltype.common.data, my.tmp)
  }
  
  
}

# export significance results
# up overlap is always more than expected by chance, down is more variable
write.table(fisher.res, file = paste0(Sys.Date(),"_DecoupleR_regulons_Fisher_Aging_FDR5_byCellType_Overlap_ZMZ_GRZ.txt"), sep = "\t", row.names = F, quote = F)

# add combined fisher FDR for help sorting
celltype.common.data$FisherCombinedPval <- celltype.common.data$p_value.grz * celltype.common.data$p_value.zmz
write.table(celltype.common.data, file = paste0(Sys.Date(),"_Aging_DecoupleR_regulons_FDR5_byCellType_Overlap_ZMZ_GRZ.txt"), sep = "\t", row.names = F, quote = F)
#####################################################################################################################

#######################
sink(file = paste(Sys.Date(),"_R_session_Info_decoupleR_JITTER_PseudoBulk_KilliBrain_Aging_Per_Cell_Type.txt", sep =""))
sessionInfo()
sink()

