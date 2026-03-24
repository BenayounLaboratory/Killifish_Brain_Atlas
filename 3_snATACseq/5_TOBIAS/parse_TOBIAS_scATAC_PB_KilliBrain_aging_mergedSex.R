setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/snATAC_Brain_Aging_Meta/Downstream_Analyses/TOBIAS/Merged_sex')
options(stringsAsFactors = F)

# Loading necessary libraries
library(pheatmap)
library(bitops)
library(grDevices)
library(ggplot2) 
library(scales) 
theme_set(theme_bw())

# 2024-05-09
# parse TOBIAS BINDetect results on pseudobulk single cell ATAC data, merged sex for increased power

###################################################################################################################
# 0. process TOBIAS results

# read bindtect analysis results
tobias.files <- c("Astrocytes_Radial_Glia_Y_vs_O_scATAC_BINDetect_JASPAR/bindetect_results.txt",
                  "Ependymal_cells_Y_vs_O_scATAC_BINDetect_JASPAR/bindetect_results.txt",
                  "GABAergic_neurons_Y_vs_O_scATAC_BINDetect_JASPAR/bindetect_results.txt",
                  "Granule_Excitatory_Neurons_Y_vs_O_scATAC_BINDetect_JASPAR/bindetect_results.txt",
                  "Microglia_Y_vs_O_scATAC_BINDetect_JASPAR/bindetect_results.txt",
                  "NSPCs_Y_vs_O_scATAC_BINDetect_JASPAR/bindetect_results.txt",
                  "Oligodendrocytes_Y_vs_O_scATAC_BINDetect_JASPAR/bindetect_results.txt",
                  "OPCs_Y_vs_O_scATAC_BINDetect_JASPAR/bindetect_results.txt",
                  "PV_interneurons_Y_vs_O_scATAC_BINDetect_JASPAR/bindetect_results.txt")

read.tpb <- function(my.file) {
  tob.res <- read.csv(my.file, header = T, sep = "\t")
  return(tob.res)
}

tobias.res.list <- lapply(tobias.files,read.tpb)
names(tobias.res.list) <- gsub("_Y_vs_O_scATAC_BINDetect_JASPAR/bindetect_results.txt","", tobias.files)

# run FDR adjustment
tob.padj <- function (tob.df) {
  tob.df$FDR <- p.adjust(tob.df[,grep("pvalue", colnames(tob.df))], method = "fdr")
  return(tob.df)
}
tobias.res.list <- lapply(tobias.res.list,tob.padj)

# use neg of the score so positive means up with aging
adj.tob.sign <- function (tob.df) {
  tob.df$Aging_Change <- -tob.df[,grep("_change", colnames(tob.df))]
  return(tob.df)
}
tobias.res.list <- lapply(tobias.res.list,adj.tob.sign)

# remove superfluous columns
select.cols.tob <- function (tob.df) {
  tob.df <- tob.df[,c("output_prefix","name","motif_id","total_tfbs","Aging_Change","FDR")]
  return(tob.df)
}
tobias.res.list <- lapply(tobias.res.list,select.cols.tob)

# get summary table out
scTobias.res.grz <- data.frame(matrix(0,0,7))
colnames(scTobias.res.grz) <- c("output_prefix","name","motif_id","total_tfbs","Aging_Change","FDR","Cell_Type")

for (i in 1:length(tobias.res.list)) {
  my.sig.data <- tobias.res.list[[i]][tobias.res.list[[i]]$FDR < 0.05,]
  my.sig.data$Cell_Type <- names(tobias.res.list)[i]
  scTobias.res.grz <- rbind(scTobias.res.grz,my.sig.data)
}
write.table(scTobias.res.grz, file = paste0(Sys.Date(),"_TOBIAS_scATAC_PB_TF_footprint_by_CellType_GRZ_FDR5.txt"), sep = "\t", quote = F, row.names = F)
write.table(scTobias.res.grz, file = paste0(Sys.Date(),"_TOBIAS_scATAC_PB_TF_footprint_by_CellType_GRZ_FDR5.xls"), sep = "\t", quote = F, row.names = F)

###################################################################################################################

###################################################################################################################
# 1. plot jitter

######### Make jitter plot of differential footprints #########
## Order by pvalue:
tobias.results <- lapply(tobias.res.list,function(x) {x[order(x$FDR),]})
n        <- sapply(tobias.results, nrow)
names(n) <- names(tobias.results)

cols <- list()
xlab <- character(length = length(tobias.results))
for(i in seq(along = tobias.results)){
  cols[[i]]      <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), n[i]) # grey60
  ind.sig.i      <- tobias.results[[i]]$FDR < 0.05
  ind.sig.i.up   <- bitAnd(tobias.results[[i]]$FDR < 0.05, tobias.results[[i]]$Aging_Change >0)>0
  ind.sig.i.down <- bitAnd(tobias.results[[i]]$FDR < 0.05, tobias.results[[i]]$Aging_Change <0)>0
  
  cols[[i]][ind.sig.i.up]   <- "#CC3333"
  cols[[i]][ind.sig.i.down] <- "#333399"
  xlab[i] <- paste(names(tobias.results)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}
names(cols) <- names(tobias.results)

pdf(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_TOBIAS_with_reg_colors_FDR5.pdf"), width = 6, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 9.5),
     ylim = c(-0.4, 0.4),
     axes = FALSE,
     xlab = "",
     ylab = "TOBIAS Footprint Score with Aging"
)
abline(h = 0)
abline(h = seq(-0.4, 0.4, by = 0.2)[-3],
       lty = "dotted",
       col = "grey")
for(i in 1:length(tobias.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(tobias.results[[i]])), amount = 0.2),
         y = rev(tobias.results[[i]]$Aging_Change),
         pch = 16,
         col = rev(cols[[i]]),
         bg = rev(cols[[i]]),
         cex = 0.75)
}
axis(1,
     at = 1:9,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-0.4, 0.4, by = 0.2))
box()
dev.off()

png(paste0(Sys.Date(),"_KilliBrain_Aging_perCellType_GRZ_stripplot_TOBIAS_with_reg_colors_FDR5.png"), width = 2600, height = 2200, res = 300, units = "px")
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 9.5),
     ylim = c(-0.4, 0.4),
     axes = FALSE,
     xlab = "",
     ylab = "TOBIAS Footprint Score with Aging"
)
abline(h = 0)
abline(h = seq(-0.4, 0.4, by = 0.2)[-3],
       lty = "dotted",
       col = "grey")
for(i in 1:length(tobias.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(tobias.results[[i]])), amount = 0.2),
         y = rev(tobias.results[[i]]$Aging_Change),
         pch = 16,
         col = rev(cols[[i]]),
         bg = rev(cols[[i]]),
         cex = 0.75)
}
axis(1,
     at = 1:9,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-0.4, 0.4, by = 0.2))
box()
dev.off()
###################################################################################################################

###################################################################################################################
# 4. bubble plot output

nr3c1.res <- scTobias.res.grz[scTobias.res.grz$name %in% "NR3C1",]

nr3c1.res$minusLog10FDR <- -log10(nr3c1.res$FDR)

#### 
my.max <- max(nr3c1.res$Aging_Change)
my.min <- min(nr3c1.res$Aging_Change)
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
nr3c1.res$Cell_Type  <- factor(nr3c1.res$Cell_Type, levels = nr3c1.res$Cell_Type)
nr3c1.res$TF_footprint   <- "NR3C1"

pdf(paste0(Sys.Date(),"TOBIAS_NR3C1_footprints_AgingScore_by_CellandSex_FDR5.pdf"), height = 3, width=5)
my.plot <- ggplot(nr3c1.res,aes(x=Cell_Type,y=TF_footprint,colour=Aging_Change,size=minusLog10FDR))+ theme_bw() + geom_point(shape = 16)
my.plot <- my.plot + ggtitle("NR3C1_MA0113.3 Footprints") + labs(x = "Killifish Brain Aging", y = "Aging DA ATAC footprints (FDR < 5%)")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled)
my.plot <- my.plot + scale_x_discrete(guide = guide_axis(angle = 45)) + scale_size(range = c(3, 10)) + scale_size_area(limits=c(0,80))
print(my.plot)
dev.off()

