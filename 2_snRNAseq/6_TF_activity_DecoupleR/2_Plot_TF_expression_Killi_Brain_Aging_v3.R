setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/TF_activity_analysis')
options(stringsAsFactors = F)

library("DESeq2")        #

library(dplyr)
library(ggplot2)
library(ggrepel)

library(bitops)
library(ComplexHeatmap)
library(circlize)         #

theme_set(theme_bw())   

# 2025-06-12
# Plot log2FC of candidate TFs


#####################################################################################################################
#### 0. Load up annotated DEseq2 objects

# Load necessary objects
load('../Differential_Expression_Age/DGE_Analysis/2024-02-28_pseudobulk_killi_cell_types_AGING_GRZ_DEseq2_objects.RData'   )
load('../Differential_Expression_Age/DGE_Analysis/2024-02-28_pseudobulk_killi_cell_types_AGING_ZMZ_DEseq2_objects.RData'   )

ls()
# [1] "deseq.res.list.genes.grz" "deseq.res.list.genes.zmz" "deseq.res.list.grz"       "deseq.res.list.zmz"         
#####################################################################################################################


#####################################################################################################################
#### 1. extract logFC of candidate TFs

TFg    <- c("LOC107375066", "helz2")
TFname <- c("nr3c1"       , "helz2") 

TFg    <- c("LOC107375066",
            "stat1",
            "LOC107376749",
            "LOC107384049", 
            "LOC107385678", 
            "LOC107387991", 
            "nfkb1",
            "nfya",
            "egr1",
            "pax6",
            "sox2",
            "irf7"
            )
TFname <- c("nr3c1",
            "stat1",
            "rela",
            "bhlhe40", 
            "stat2", 
            "irf3",
            "nfkb1",
            "nfya",
            "egr1",
            "pax6",
            "sox2",
            "irf7"
            ) 


# make a table for plotting
my.cell.types <- sort(names(deseq.res.list.genes.grz))

my.plot.table.grz <- data.frame(matrix(0,length(TFg),length(my.cell.types)))
my.plot.table.zmz <- data.frame(matrix(0,length(TFg),length(my.cell.types)))
colnames(my.plot.table.grz) <- my.cell.types
colnames(my.plot.table.zmz) <- my.cell.types
rownames(my.plot.table.grz) <- TFname
rownames(my.plot.table.zmz) <- TFname

my.sig.table.grz <- data.frame(matrix(0,length(TFg),length(my.cell.types)))
my.sig.table.zmz <- data.frame(matrix(0,length(TFg),length(my.cell.types)))
colnames(my.sig.table.grz) <- my.cell.types
colnames(my.sig.table.zmz) <- my.cell.types
rownames(my.sig.table.grz) <- TFname
rownames(my.sig.table.zmz) <- TFname



for (i in 1:length(my.cell.types)) {
  
  # deseq2 data in GRZ and ZMZ
  grz.res <- deseq.res.list.genes.grz[[i]]
  zmz.res <- deseq.res.list.genes.zmz[[i]]
  
  for (j in 1:length(TFg)) {
    
    # 
    my.tf <- TFg[j]
    
    my.plot.table.grz[j,i]  <- grz.res[my.tf,]$log2FoldChange
    my.plot.table.zmz[j,i]  <- zmz.res[my.tf,]$log2FoldChange
    
    my.sig.table.grz[j,i]  <- grz.res[my.tf,]$padj
    my.sig.table.zmz[j,i]  <- zmz.res[my.tf,]$padj
  }
  
}


my.plot.table.grz

### get order

grz.heat <- Heatmap(my.plot.table.grz, border = T, rect_gp = gpar(col = "grey", lwd = 0.5),
                    cluster_rows = F, cluster_columns = F,
                    column_title = "DESeq2 log2FC (GRZ)",
                    col = colorRamp2(c(-0.5,0,0.5), c("darkblue","white", "firebrick4"), transparency = 0, space = "LAB"),
                    width  = ncol(my.plot.table.grz)*unit(7, "mm"),
                    height = nrow(my.plot.table.grz)*unit(7, "mm"))
zmz.heat <- Heatmap(my.plot.table.zmz, border = T, rect_gp = gpar(col = "grey", lwd = 0.5),
                    cluster_rows = F, cluster_columns = F,
                    column_title = "DESeq2 log2FC (ZMZ)",
                    col = colorRamp2(c(-0.5,0,0.5), c("darkblue","white", "firebrick4"), transparency = 0, space = "LAB"),
                    width  = ncol(my.plot.table.zmz)*unit(7, "mm"),
                    height = nrow(my.plot.table.zmz)*unit(7, "mm"))


pdf(paste0(Sys.Date(),"_DESeq2_log2FC_FDR_TopTFs_decoupleR_ZMZ_GRZ_PANELED.pdf"), width = 8, height = 10)
grz.heat + zmz.heat
dev.off()

write.table(my.sig.table.grz, file = paste0(Sys.Date(),"_DESeq2_FDR_TopTFs_decoupleR_GRZ.txt"), sep = "\t", quote = F)
write.table(my.sig.table.zmz, file = paste0(Sys.Date(),"_DESeq2_FDR_TopTFs_decoupleR_ZMZ.txt"), sep = "\t", quote = F)
write.table(my.plot.table.grz, file = paste0(Sys.Date(),"_DESeq2_log2FC_TopTFs_FDR_TopTFs_decoupleR_GRZ.txt"), sep = "\t", quote = F)
write.table(my.plot.table.zmz, file = paste0(Sys.Date(),"_DESeq2_log2FC_TopTFs_FDR_TopTFs_decoupleR_ZMZ.txt"), sep = "\t", quote = F)


#####################################################################################################################

#######################
sink(file = paste(Sys.Date(),"_R_session_Info_decoupleR_PseudoBulk_KilliBrain_Aging_Per_Cell_Type.txt", sep =""))
sessionInfo()
sink()

