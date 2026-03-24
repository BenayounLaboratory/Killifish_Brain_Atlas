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


# 2025-06-23
#  plot top genes selected in RNAscope
# include PB heatmap
##### - LOC107374900	optineurin
##### - helz2	helicase with zinc finger domain 2
##### - stat1	signal transducer and activator of transcription 1
##### - znfx1	zinc finger, NFX1-type containing 1
##### - LOC107387991	irf3- interferon regulatory factor 3
##### - ifih1	interferon induced with helicase C domain 1

# 2025-06-30
# replot scatters with bigger point diameter
##########################################################################################################################################

##########################################################################################################################################
### 1. load and prepare data

#### in alphabetical order
# rnascope.targets <- c("LOC107374900","helz2","stat1","znfx1","LOC107387991","ifih1")
rnascope.targets <- c("helz2",
                      "ifih1",
                      "LOC107387991", ## irf3
                      "LOC107374900", ## opt
                      "stat1",
                      "znfx1")

# load deseq2
load('2024-02-28_pseudobulk_killi_cell_types_AGING_GRZ_DEseq2_objects.RData')
load('2024-02-28_pseudobulk_killi_cell_types_AGING_ZMZ_DEseq2_objects.RData')
load('2024-02-28_pseudobulk_killi_cell_types_AGING_GRZ_VST_data_objects.RData')
load('2024-02-28_pseudobulk_killi_cell_types_AGING_ZMZ_VST_data_objects.RData')

ls()
# [1] "deseq.res.list.genes.grz" "deseq.res.list.genes.zmz" "deseq.res.list.grz"       "deseq.res.list.zmz"       "vst.cts.grz"             
# [6] "vst.cts.zmz"    


##########################################################################################################################################
### 2. Speaman rank correlation analysis

### Prepare output objects
deseq.merged.list <- vector(mode = "list", length(deseq.res.list.genes.grz))
names(deseq.merged.list) <- names(deseq.res.list.genes.grz)

cor.test.res           <- data.frame(matrix(0,length(deseq.res.list.genes.grz),2), row.names = names(deseq.merged.list))
colnames(cor.test.res) <- c("Rho", "Pval")

for (i in 1: length((deseq.merged.list))) {
  
  my.cell.type <- names(deseq.merged.list)[i]
  
  # merge DEseq2 result tables
  deseq.merged.list[[i]] <- merge(deseq.res.list.genes.grz[[i]], deseq.res.list.genes.zmz[[i]], by = "row.names", suffixes = c(".grz",".zmz"))
  
  # correlation test
  test.res <- cor.test(deseq.merged.list[[i]]$log2FoldChange.grz, deseq.merged.list[[i]]$log2FoldChange.zmz, method = "spearman")
  cor.test.res[i,1] <- test.res$estimate
  cor.test.res[i,2] <- test.res$p.value
  
  # scatter plot
  pdf(paste0(Sys.Date(),"_", my.cell.type, "aging_GRZ_ZMZ_scatter_with_RNAscope_targets.pdf"), height = 5, width = 5)
  smoothScatter(deseq.merged.list[[i]]$log2FoldChange.grz,
                deseq.merged.list[[i]]$log2FoldChange.zmz,
                xlab = "log2FC per week (GRZ)",
                ylab = "log2FC per week (ZMZ1001)",
                main = my.cell.type,
                xlim = c(-0.7,0.7),
                ylim = c(-0.4,0.4))
  abline(h = 0, col = "grey", lty = "dashed")
  abline(v = 0, col = "grey", lty = "dashed")
  points(deseq.merged.list[[i]]$log2FoldChange.grz[deseq.merged.list[[i]]$Row.names %in% rnascope.targets],
         deseq.merged.list[[i]]$log2FoldChange.zmz[deseq.merged.list[[i]]$Row.names %in% rnascope.targets],
         col = "red", cex = 2)
  text(-0.7, 0.4 , paste0("Rho = ", signif(test.res$estimate,2)),pos = 4)
  text(-0.7, 0.37, paste0("p = "  , signif(test.res$p.value,2)),pos = 4)
  dev.off()
  
}

write.table(cor.test.res, file = paste0(Sys.Date(),"_Aging_DEseq2_SpearmanRho_byCellType_ZMZ_GRZ.txt"), sep = "\t", row.names = T, quote = F)

######## 
cor.test.res$minuslog10p <- -log10(cor.test.res$Pval)
my.cors <- cbind(rownames(cor.test.res),cor.test.res)
colnames(my.cors)[1] <- "CellType"
my.cors$condition <- rep("GRZvZMZ",nrow(my.cors))

my.max <- 0.3
my.min <- 0
my.values <- c(0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
my.cors$condition <- factor(my.cors$condition, levels = unique(my.cors$condition))
my.cors$CellType   <- factor(my.cors$CellType, levels = rev(unique(my.cors$CellType)))

pdf(paste0(Sys.Date(),"_Spearman_CellTypes.pdf"),height = 5, width=5)
my.plot <- ggplot(my.cors,aes(x=condition,y=CellType,colour=Rho,size=minuslog10p))+ theme_bw()+ geom_point(shape = 16)
my.plot <- ggplot(my.cors,aes(x=condition,y=CellType,colour=Rho,size=minuslog10p))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("cor") + labs(x = "", y = "CellType")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(my.min,my.max))
print(my.plot)
dev.off()


##########################################################################################################################################
### 3. RNAscope expression analysis

for (i in 1:length(vst.cts.grz)) {
  
  my.cell.type <- names(vst.cts.grz)[i]
  
  grz.av     <- data.frame(row.names = rnascope.targets)
  grz.av$YF  <- apply(vst.cts.grz[[i]][rnascope.targets,c(1:4)],1,median)
  grz.av$MF  <- apply(vst.cts.grz[[i]][rnascope.targets,c(9:11)],1,median)
  grz.av$OF  <- apply(vst.cts.grz[[i]][rnascope.targets,c(15:17)],1,median)
  grz.av$YM  <- apply(vst.cts.grz[[i]][rnascope.targets,c(5:8)],1,median)
  grz.av$MM  <- apply(vst.cts.grz[[i]][rnascope.targets,c(12:14)],1,median)
  grz.av$OM  <- apply(vst.cts.grz[[i]][rnascope.targets,c(18:20)],1,median)
  
  

  zmz.av     <- data.frame(row.names = rnascope.targets)
  zmz.av$YF  <- apply(vst.cts.zmz[[i]][rnascope.targets,c(1:3)],1,median)
  zmz.av$MF  <- apply(vst.cts.zmz[[i]][rnascope.targets,c(7:9)],1,median)
  zmz.av$OF  <- apply(vst.cts.zmz[[i]][rnascope.targets,c(13:15)],1,median)
  zmz.av$GF  <- apply(vst.cts.zmz[[i]][rnascope.targets,c(19:21)],1,median)
  zmz.av$YM  <- apply(vst.cts.zmz[[i]][rnascope.targets,c(4:6)],1,median)
  zmz.av$MM  <- apply(vst.cts.zmz[[i]][rnascope.targets,c(10:12)],1,median)
  zmz.av$OM  <- apply(vst.cts.zmz[[i]][rnascope.targets,c(16:18)],1,median)
  zmz.av$GM  <- apply(vst.cts.zmz[[i]][rnascope.targets,c(22:24)],1,median)
  
  
  grz.heat <- Heatmap(t(scale(t(grz.av))), 
                      col = colorRamp2(c(-2,0, 2), c("darkblue", "white", "#CC3333"), transparency = 0, space = "LAB"),
                      border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
                      cluster_rows = F,
                      cluster_columns = F,
                      column_title = paste0(my.cell.type, "(GRZ)"),
                      width  = ncol(grz.av)*unit(7, "mm"),
                      height = nrow(grz.av)*unit(7, "mm"))
  
  zmz.heat <- Heatmap(t(scale(t(zmz.av))), 
                      col = colorRamp2(c(-2,0, 2), c("darkblue", "white", "#CC3333"), transparency = 0, space = "LAB"),
                      border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
                      cluster_rows = F,
                      cluster_columns = F,
                      column_title = paste0(my.cell.type, "(ZMZ1001)"),
                      width  = ncol(zmz.av)*unit(7, "mm"),
                      height = nrow(zmz.av)*unit(7, "mm"))
  
  
  pdf(paste0(Sys.Date(),"_", my.cell.type, "_aging_GRZ_ZMZ_RNAscope_median_heatmap.pdf"), width = 12, height = 5)
  print(grz.heat + zmz.heat)
  dev.off()
  
}
############################################################################################

#######################
sink(file = paste(Sys.Date(),"_Correlation_RNAscope_Analysis_BrainAtlas_session_Info.txt", sep =""))
sessionInfo()
sink()