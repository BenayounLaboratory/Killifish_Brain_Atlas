setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Mouse_Datasets_for_Comparison/TF_Activity_analysis')
options(stringsAsFactors = F)

library("DESeq2")        #
library("decoupleR")     #
library("OmnipathR")     #

library(dplyr)
library(tibble)
library(ggplot2)
library(beeswarm)
library(ggrepel)


library(ComplexHeatmap)
library(circlize)

theme_set(theme_bw())   

# 2024-03-25
# Use decoupleR bioconductor package (and mouse TF target DB)
# to infer differential TF activity in brain aging data
# to compare with Killi data
# Use pseudoblulk analysis 
# https://www.bioconductor.org/packages/release/bioc/vignettes/decoupleR/inst/doc/tf_bk.html

# 2024-03-28
# Include Hahn and Allen datasets to complement analyses

# 2024-04-12
# Include Ximerakis and Ogrodnik data

# 2024-11-08
# revise color scale to reflect signifiance level

#####################################################################################################################
#### 0. Load up annotated DEseq2/VST objects

# Load necessary objects
load('../Muscat_DESeq2/Allen_PFC/2024-03-27_PB_Hajdarovic_CellType_GSE188646_DEseq2_objects.RData'   )
pfc.deseq.res.list  <- deseq.res.list ; rm(deseq.res.list)

load('../Muscat_DESeq2/Buckley_SVZ/2024-03-25_PB_Buckley_CellType_DEseq2_objects.RData'   )
svz.deseq.res.list  <- deseq.res.list  ; rm(deseq.res.list)

load('../Muscat_DESeq2/Hajdarovic_Hypothalamus/2024-03-22_PB_Hajdarovic_CellType_GSE188646_DEseq2_objects.RData'   )
hyp.deseq.res.list  <- deseq.res.list  ; rm(deseq.res.list)

load('../Muscat_DESeq2/Hahn_Hippocampus_CaudatePutamen/2024-03-28_PB_Hahn_GSE212576_Caudate_Putamen_DEseq2_objects.RData'   )
cp.deseq.res.list  <- deseq.res.list  ; rm(deseq.res.list)

load('../Muscat_DESeq2/Hahn_Hippocampus_CaudatePutamen/2024-03-28_PB_Hahn_GSE212576_Hippocampus_DEseq2_objects.RData'   )
hippo.deseq.res.list  <- deseq.res.list  ; rm(deseq.res.list)

load('../Muscat_DESeq2/Ximerakis_Brain/2024-04-11_PB_SCP263_Ximerakis_DEseq2_objects.RData'   )
brain.deseq.res.list  <- deseq.res.list  ; rm(deseq.res.list)

#There is PISD uppercase as a gene which creates a duplicate, we will remove that line
for (i in 1:length(brain.deseq.res.list)) {
  brain.deseq.res.list[[i]] <- brain.deseq.res.list[[i]][!(rownames(brain.deseq.res.list[[i]]) %in% "PISD"),]
}

load('../Muscat_DESeq2/Ogrodnik_GSE161340/2024-04-12_PB_Ogrodnik_Hippocampus_DEseq2_objects.RData'   )
hippo_o.deseq.res.list  <- deseq.res.list  ; rm(deseq.res.list)

ls()
# [1] "brain.deseq.res.list" "cp.deseq.res.list"    "hippo.deseq.res.list" "hippo_o.deseq.res.list" "hyp.deseq.res.list"   "pfc.deseq.res.list"   "svz.deseq.res.list"
#####################################################################################################################


#####################################################################################################################
#### 1. Process TF network data from mouse model
mm.net <- get_collectri(organism='mouse', split_complexes=FALSE)
mm.net
# # A tibble: 42,595 × 3
#    source target   mor
#    <chr>  <chr>  <dbl>
# 1 MYC    TERT       1
# 2 SPI1   BGLAP      1
# 3 SMAD3  JUN        1
# 4 SMAD4  JUN        1
# 5 STAT5A IL2        1
# 6 STAT5B IL2        1
# 7 RELA   FAS        1
# 8 WT1    NR0B1      1
# 9 NR0B2  CASP1      1
# 10 SP1    ALDOA      1

# Will need to make gene names match mouse nomenclature
#####################################################################################################################


#####################################################################################################################
#### 2. Calculate TF network data

n_tfs <- 10

Top_TFs_res.list        <- vector(mode = "list", length = 7)
names(Top_TFs_res.list) <- c("Brain_Ximerakis","Caudate_Putamen_Hahn", "Hippocampus_Hahn", "Hippocampus_Ogrdnik", "Hypothalamus_Hajdarovic", "Prefrontal_Cortex_Allen","SVZ_Buckley")

my.deseq2.list <- list(brain.deseq.res.list, cp.deseq.res.list, hippo.deseq.res.list, hippo_o.deseq.res.list, hyp.deseq.res.list, pfc.deseq.res.list,svz.deseq.res.list )
names(my.deseq2.list) <- c("Brain_Ximerakis","Caudate_Putamen_Hahn", "Hippocampus_Hahn", "Hippocampus_Ogrdnik", "Hypothalamus_Hajdarovic", "Prefrontal_Cortex_Allen","SVZ_Buckley")


# Run decoupleR for each tissue
for (k in 1:length(my.deseq2.list)) {
  # infer pathway activities from the t-values of the DEGs with aging
  # for each cell type
  
  Top_TFs_res.list[[k]] <- data.frame(matrix(0,0,7))
  colnames(Top_TFs_res.list[[k]]) <- c("statistic", "source","condition","score","p_value","Reg_Rank","Cell_Type")
  
  # Loop over DEseq2 results
  for (i in 1:length(my.deseq2.list[[k]])) {
    
    #network has mosue gene names as upper case
    rownames(my.deseq2.list[[k]][[i]]) <- toupper(rownames(my.deseq2.list[[k]][[i]]))
    
    # Run fgsea scoring
    contrast_acts <- run_fgsea(mat     = my.deseq2.list[[k]][[i]][, 'stat', drop=FALSE],
                               net     = mm.net                                  , 
                               minsize = 5                                             )
    # contrast_acts
    # A tibble: 1,206 × 5
    # statistic  source condition  score p_value
    # <chr>      <chr>  <chr>      <dbl>   <dbl>
    # 1 fgsea      ABL1   stat      -0.474  0.311 
    # 2 norm_fgsea ABL1   stat      -1.11   0.311 
    # 3 fgsea      AHR    stat      -0.232  0.362 
    # 4 norm_fgsea AHR    stat      -1.02   0.362 
    # 5 fgsea      AIRE   stat       0.459  0.253 

    # We select the norm_fgsea activities and then we show changes in activity with aging
    # Filter norm_fgsea
    tf.scores         <- contrast_acts[contrast_acts$statistic == 'norm_fgsea',]
    tf.sig            <- tf.scores[tf.scores$p_value < 0.1,]
    tf.sig$Reg_Rank   <- NA
    tf.sig$Cell_Type  <- names(my.deseq2.list[[k]])[i]
    
    # Only keep TF regulons who TF is expressed in cell type
    tf.sig <- tf.sig[tf.sig$source %in% row.names(my.deseq2.list[[k]][[i]]),]
    
    # Add rank information
    msk <- tf.sig$score > 0
    tf.sig$Reg_Rank[msk ]  <- round(rank(-tf.sig[msk, 'score']))
    tf.sig$Reg_Rank[!msk]  <- round(rank(-abs(tf.sig[!msk, 'score'])))
    
    Top_TFs_res.list[[k]]  <- rbind(Top_TFs_res.list[[k]] ,tf.sig)
    
    # Filter top significant TFs in both signs
    tfs <- tf.sig %>%
      arrange(Reg_Rank) %>%
      head(n_tfs) %>%
      pull(source)
    
    f_contrast_acts <- tf.sig %>% filter(source %in% tfs)
    
    # Plot top each direction
    tf.reg.plot <- ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + geom_bar(aes(fill = score), stat = "identity") 
    tf.reg.plot <- tf.reg.plot + scale_fill_gradient2(low = "#333399", mid = "whitesmoke", high = "#CC3333", midpoint = 0, limits = c(-2.5,2.5)) 
    tf.reg.plot <- tf.reg.plot  + theme(axis.title = element_text(face = "bold", size = 12),
                                        axis.text.x = element_text(angle = 45, hjust = 1, size =10),
                                        axis.text.y = element_text(size =10, face= "bold"),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank() ) 
    tf.reg.plot <- tf.reg.plot + xlab("Enriched TF regulons") + ylab("decoupleR score") + theme(plot.title = element_text(hjust = 0.5))
    tf.reg.plot <- tf.reg.plot + ylim(c(-2.5,2.5)) + ggtitle(paste(names(my.deseq2.list[[k]])[i],names(my.deseq2.list)[k])) + coord_flip()
    
    pdf(paste0(Sys.Date(),"_Barplot_decoupleR_fgsea_Top_", n_tfs , "_TF_regulons_", names(my.deseq2.list[[k]])[i], "_",names(my.deseq2.list)[k],"Aging_FDR10.pdf"), width = 4, height = 5)
    print(tf.reg.plot)
    dev.off()
    
  }
  
  write.table(Top_TFs_res.list[[k]], file = paste0(Sys.Date(),"_decoupleR_fgsea_TFregulons_by_CellType_", names(my.deseq2.list)[k], "_FDR10.txt"), sep = "\t", quote = F, row.names = F)
  
  # Reccurrent TF in more than half the cell types
  concat <- function (x) {
    paste0(x, collapse = ",")
  }
  
  Top_TFs_res.list[[k]]$sign <- sign(Top_TFs_res.list[[k]]$score)
  my.num.comb   <- aggregate(Top_TFs_res.list[[k]]$Cell_Type, by = list("TF_regulon" = Top_TFs_res.list[[k]]$source), FUN = length)
  colnames(my.num.comb)[2] <- "N_CellType"
  my.signs.comb <- aggregate(Top_TFs_res.list[[k]]$sign, by = list("TF_regulon" =Top_TFs_res.list[[k]]$source), FUN = concat)
  colnames(my.signs.comb)[2] <- "Signs"
  my.sum.comb <- merge(my.num.comb,my.signs.comb)
  
  my.sum.comb[my.sum.comb$N_CellType > 1,]
  # recurrent perturbed regulons
  #     TF_regulon N_CellType     Signs
  # 9         ATF4          3     1,1,1
  # 16       CEBPA          3     1,1,1
  # 36        EGR1          3    1,1,-1

  pdf(paste0(Sys.Date(),"_decoupleR_fgsea_TFregulons_CellTypeSharing_Boxplot_", names(my.deseq2.list)[k], ".pdf"), height = 5, width = 3.5)
  boxplot( my.sum.comb$N_CellType, col = "purple", main = "DecoupleR Regulon Sharing (QC cell types)", ylab = "Number of cell types", outline = F, ylim = c(0,10))
  beeswarm(my.sum.comb$N_CellType[my.sum.comb$N_CellType > 1], method = "compactswarm",corral = "gutter", add = T, pch = 16)
  points(1, my.sum.comb[my.sum.comb$TF_regulon %in% "NR3C1",]$N_CellType,pch = 16, col = "red")
  text(1.2, my.sum.comb[my.sum.comb$TF_regulon %in% "NR3C1",]$N_CellType,"Nr3c1", pos = 4, cex = 0.5, col = "red")
  dev.off()
  
}

save(Top_TFs_res.list, file = paste0(Sys.Date(),"_decoupleR_fgsea_TFregulons_by_CellType_MouseBrain_FDR10.RData"))

options(java.parameters = "-Xmx16g" )
require(openxlsx)

write.xlsx(Top_TFs_res.list, rowNames = F, file = paste0(Sys.Date(),"_Mouse_Brain_Aging_DecoupleR_Results_FDR10.xlsx"))
########################################################################################################################

########################################################################################################################
# Plot summary
load('2024-04-12_decoupleR_fgsea_TFregulons_by_CellType_MouseBrain_FDR10.RData')

# TFs of interest
my.tfs <- c("NR3C1","STAT1","RELA","BHLHE40","STAT2", "IRF3","NFKB1","NFYA","EGR1","PAX6","SOX2","IRF7")

plot.list        <- vector(mode = "list",length(Top_TFs_res.list))
names(plot.list) <- names(Top_TFs_res.list)
  
for (i in 1:length(Top_TFs_res.list)) {
  
  # output data
  my.cell.types              <- sort(unique(Top_TFs_res.list[[i]]$Cell_Type))
  my.plot.table.mm           <- data.frame(matrix(0,length(my.tfs),length(my.cell.types)))
  colnames(my.plot.table.mm) <- my.cell.types
  rownames(my.plot.table.mm) <- my.tfs

  
  for (j in 1:length(my.tfs)) {
    
    my.cur.tf <- Top_TFs_res.list[[i]][Top_TFs_res.list[[i]]$source %in% my.tfs[j],]
    
    if(nrow(my.cur.tf) >0) {
      
      for (k in 1:nrow(my.cur.tf)) {
        # populate based on sign, significance and Cell Type
        
        if ( (my.cur.tf[k,]$sign > 0) && (my.cur.tf[k,]$p_value < 0.05)) {
          my.plot.table.mm[j,my.cur.tf[k,]$Cell_Type] <- 1
        } else if ( (my.cur.tf[k,]$sign > 0) && (my.cur.tf[k,]$p_value < 0.1)) {
          my.plot.table.mm[j,my.cur.tf[k,]$Cell_Type] <- 0.5
        } else if ( (my.cur.tf[k,]$sign < 0) && (my.cur.tf[k,]$p_value < 0.05)) {
          my.plot.table.mm[j,my.cur.tf[k,]$Cell_Type] <- -1
        } else if ( (my.cur.tf[k,]$sign < 0) && (my.cur.tf[k,]$p_value < 0.1)) {
          my.plot.table.mm[j,my.cur.tf[k,]$Cell_Type] <- -0.5
        }
       
      }
    }
  }
  
  pdf(paste0(Sys.Date(),"_KilliRecurrent_decoupleR_Regulons_byCellTypes_",names(Top_TFs_res.list)[i],"_FDR10.pdf"), width = 6, height = 8)
  cur.heat <- Heatmap(my.plot.table.mm, border = T, rect_gp = gpar(col = "grey65", lwd = 0.5),
                      cluster_rows = F, cluster_columns = F,
                      column_title = names(Top_TFs_res.list)[i],
                      width  = ncol(my.plot.table.mm)*unit(7, "mm"),
                      height = nrow(my.plot.table.mm)*unit(7, "mm"),
                      col = colorRamp2(c(-1, 0, 1), c("blue", "grey95", "red")),
                      column_names_rot = 45)
  print(cur.heat)
  dev.off()
  
  plot.list[[i]] <- cur.heat
}

my.plot <- plot.list[[1]] + plot.list[[2]] + plot.list[[3]] + plot.list[[4]] + plot.list[[5]] + plot.list[[6]]+ plot.list[[7]]

pdf(paste0(Sys.Date(),"_KilliRecurrent_decoupleR_Regulons_byCellTypes_ALL_DATA_FDR10_color.pdf"), width = 20, height = 8)
print(my.plot)
dev.off()
#####################################################################################################################

#######################
sink(file = paste(Sys.Date(),"_R_session_Info_decoupleR_MouseBrain_Aging_Per_Cell_Type.txt", sep =""))
sessionInfo()
sink()

