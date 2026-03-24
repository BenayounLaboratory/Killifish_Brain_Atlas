setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/Misc/GTEX_Human_Brain/DecoupleR')
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

# 2025-03-07
# Use decoupleR bioconductor package (and human TF target DB)
# to infer differential TF activity in GTEx brain aging data
# to compare with Killi data

#####################################################################################################################
#### 0. Load up annotated DEseq2/VST objects

# Load necessary objects
load('../DESeq2/2025-03-04GTEX_Brain_AGING_DEseq2_objects.RData'   )

ls()
# [1] "deseq.res.list.gtex"
#####################################################################################################################


#####################################################################################################################
#### 1. Process TF network data from mouse model
hs.net <- get_collectri(organism='human', split_complexes=FALSE)
hs.net
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
#####################################################################################################################


#####################################################################################################################
#### 2. Calculate TF network data

n_tfs <- 10

Top_TFs_res <- data.frame(matrix(0,0,7))
colnames(Top_TFs_res) <- c("statistic", "source","condition","score","p_value","Reg_Rank","Region")

# Loop over DEseq2 results
for (i in 1:length(deseq.res.list.gtex)) {
  
  # Run fgsea scoring
  contrast_acts <- run_fgsea(mat     = deseq.res.list.gtex[[i]][, 'stat', drop=FALSE],
                             net     = hs.net                                  , 
                             minsize = 5                                             )
  # contrast_acts
  # # A tibble: 1,472 × 5
  # statistic  source condition score  p_value
  # <chr>      <chr>  <chr>     <dbl>    <dbl>
  #   1 fgsea      ABL1   stat      0.653 0.0132  
  # 2 norm_fgsea ABL1   stat      1.82  0.0132  
  # 3 fgsea      AEBP1  stat      0.548 0.380   
  # 4 norm_fgsea AEBP1  stat      1.19  0.380     
  
  # We select the norm_fgsea activities and then we show changes in activity with aging
  # Filter norm_fgsea
  tf.scores         <- contrast_acts[contrast_acts$statistic == 'norm_fgsea',]
  tf.sig            <- tf.scores[tf.scores$p_value < 0.1,]
  tf.sig$Reg_Rank   <- NA
  tf.sig$Region  <- names(deseq.res.list.gtex)[i]
  
  # Only keep TF regulons who TF is expressed in cell type
  tf.sig <- tf.sig[tf.sig$source %in% row.names(deseq.res.list.gtex[[i]]),]
  
  # Add rank information
  msk <- tf.sig$score > 0
  tf.sig$Reg_Rank[msk ]  <- round(rank(-tf.sig[msk, 'score']))
  tf.sig$Reg_Rank[!msk]  <- round(rank(-abs(tf.sig[!msk, 'score'])))
  
  Top_TFs_res  <- rbind(Top_TFs_res ,tf.sig)
  
  # Filter top significant TFs in both signs
  tfs <- tf.sig %>%
    arrange(Reg_Rank) %>%
    head(n_tfs) %>%
    pull(source)
  
  f_contrast_acts <- tf.sig %>% filter(source %in% tfs)
  
  # Plot top each direction
  tf.reg.plot <- ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + geom_bar(aes(fill = score), stat = "identity") 
  tf.reg.plot <- tf.reg.plot + scale_fill_gradient2(low = "#333399", mid = "whitesmoke", high = "#CC3333", midpoint = 0, limits = c(-5,5)) 
  tf.reg.plot <- tf.reg.plot  + theme(axis.title = element_text(face = "bold", size = 12),
                                      axis.text.x = element_text(angle = 45, hjust = 1, size =10),
                                      axis.text.y = element_text(size =10, face= "bold"),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank() ) 
  tf.reg.plot <- tf.reg.plot + xlab("Enriched TF regulons") + ylab("decoupleR score") + theme(plot.title = element_text(hjust = 0.5))
  tf.reg.plot <- tf.reg.plot + ylim(c(-5,5)) + ggtitle(names(deseq.res.list.gtex)[i]) + coord_flip()
  # tf.reg.plot
  
  pdf(paste0(Sys.Date(),"_Barplot_decoupleR_fgsea_Top_", n_tfs , "_TF_regulons_", names(deseq.res.list.gtex)[i], "_Aging_FDR10.pdf"), width = 4, height = 5)
  print(tf.reg.plot)
  dev.off()
  
}

write.table(Top_TFs_res, file = paste0(Sys.Date(),"_decoupleR_fgsea_TFregulons_by_Region_GTEX_FDR10.txt"), sep = "\t", quote = F, row.names = F)

# Reccurrent TF in more than half the cell types
concat <- function (x) {
  paste0(x, collapse = ",")
}

Top_TFs_res$sign <- sign(Top_TFs_res$score)
my.num.comb   <- aggregate(Top_TFs_res$Region, by = list("TF_regulon" = Top_TFs_res$source), FUN = length)
colnames(my.num.comb)[2] <- "N_CellType"
my.signs.comb <- aggregate(Top_TFs_res$sign, by = list("TF_regulon" =Top_TFs_res$source), FUN = concat)
colnames(my.signs.comb)[2] <- "Signs"
my.sum.comb <- merge(my.num.comb,my.signs.comb)

my.sum.comb[my.sum.comb$N_CellType > 12,]
# recurrent perturbed regulons
#        TF_regulon N_CellType                                  Signs
# 40       CEBPA         13              1,1,1,1,1,1,1,1,1,1,1,1,1
# 47       CIITA         13              1,1,1,1,1,1,1,1,1,1,1,1,1
# 214       IRF1         13              1,1,1,1,1,1,1,1,1,1,1,1,1
# 217       IRF4         13              1,1,1,1,1,1,1,1,1,1,1,1,1
# 221       IRF8         13              1,1,1,1,1,1,1,1,1,1,1,1,1
# 325      NFKB1         13              1,1,1,1,1,1,1,1,1,1,1,1,1
# 356      NR3C1         13              1,1,1,1,1,1,1,1,1,1,1,1,1   *****
# 363       NRF1         13 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
# 408      PPARG         13              1,1,1,1,1,1,1,1,1,1,1,1,1
# 422        REL         13              1,1,1,1,1,1,1,1,1,1,1,1,1
# 423       RELA         13              1,1,1,1,1,1,1,1,1,1,1,1,1
# 428       RFX5         13              1,1,1,1,1,1,1,1,1,1,1,1,1
# 483       SPI1         13              1,1,1,1,1,1,1,1,1,1,1,1,1
# 490      STAT1         13              1,1,1,1,1,1,1,1,1,1,1,1,1
# 594     ZNF384         13              1,1,1,1,1,1,1,1,1,1,1,1,1

pdf(paste0(Sys.Date(),"_decoupleR_fgsea_TFregulons_Region_Sharing_Boxplot_GTEx.pdf"), height = 5, width = 3.5)
boxplot( my.sum.comb$N_CellType, col = "purple", main = "DecoupleR Regulon Sharing (QC cell types)", ylab = "Number of region", outline = F, ylim = c(0,15))
beeswarm(my.sum.comb$N_CellType[my.sum.comb$N_CellType > 10], method = "compactswarm",corral = "gutter", add = T, pch = 16)
points(1, my.sum.comb[my.sum.comb$TF_regulon %in% "NR3C1",]$N_CellType,pch = 16, col = "red")
text(1.2, my.sum.comb[my.sum.comb$TF_regulon %in% "NR3C1",]$N_CellType,"Nr3c1", pos = 4, cex = 0.5, col = "red")
dev.off()

save(Top_TFs_res, file = paste0(Sys.Date(),"_decoupleR_fgsea_TFregulons_by_Region_GTEx_FDR10.RData"))

options(java.parameters = "-Xmx16g" )
require(openxlsx)

write.xlsx(Top_TFs_res, rowNames = F, file = paste0(Sys.Date(),"_Human_Brain_Aging_GTEX_DecoupleR_Results_FDR10.xlsx"))
########################################################################################################################

########################################################################################################################
# Plot summary
load('2025-03-07_decoupleR_fgsea_TFregulons_by_Region_GTEx_FDR10.RData')

# TFs of interest
my.tfs <- c("NR3C1","STAT1","RELA","BHLHE40","STAT2", "IRF3","NFKB1","NFYA","EGR1","PAX6","SOX2","IRF7")

# output data
my.regions                 <- sort(unique(Top_TFs_res$Region))
my.plot.table.hs           <- data.frame(matrix(0,length(my.tfs),length(my.regions)))
colnames(my.plot.table.hs) <- my.regions
rownames(my.plot.table.hs) <- my.tfs


for (j in 1:length(my.tfs)) {
  
  my.cur.tf <- Top_TFs_res[Top_TFs_res$source %in% my.tfs[j],]
  
  if(nrow(my.cur.tf) >0) {
    
    for (k in 1:nrow(my.cur.tf)) {
      # populate based on sign, significance and Cell Type
      
      if ( (my.cur.tf[k,]$sign > 0) && (my.cur.tf[k,]$p_value < 0.05)) {
        my.plot.table.hs[j,my.cur.tf[k,]$Region] <- 1
      } else if ( (my.cur.tf[k,]$sign > 0) && (my.cur.tf[k,]$p_value < 0.1)) {
        my.plot.table.hs[j,my.cur.tf[k,]$Region] <- 0.5
      } else if ( (my.cur.tf[k,]$sign < 0) && (my.cur.tf[k,]$p_value < 0.05)) {
        my.plot.table.hs[j,my.cur.tf[k,]$Region] <- -1
      } else if ( (my.cur.tf[k,]$sign < 0) && (my.cur.tf[k,]$p_value < 0.1)) {
        my.plot.table.hs[j,my.cur.tf[k,]$Region] <- -0.5
      }
      
    }
  }
}

my.plot.table.hs <- as.matrix(my.plot.table.hs)

pdf(paste0(Sys.Date(),"_KilliRecurrent_decoupleR_Regulons_by_Region_Human_Brain_Aging_GTEx_FDR10.pdf"), width = 6, height = 8)
cur.heat <- Heatmap(my.plot.table.hs, border = T, rect_gp = gpar(col = "grey65", lwd = 0.5),
                    cluster_rows = F, cluster_columns = F,
                    column_title = "GTEx Brain Aging",
                    name = "Score",
                    width  = ncol(my.plot.table.hs)*unit(7, "mm"),
                    height = nrow(my.plot.table.hs)*unit(7, "mm"),
                    col = colorRamp2(c(-1, 0, 1), c("blue", "grey95", "red")),
                    column_names_rot = 45)
print(cur.heat)
dev.off()
#####################################################################################################################

#######################
sink(file = paste(Sys.Date(),"_R_session_Info_decoupleR_MouseBrain_Aging_Per_Region.txt", sep =""))
sessionInfo()
sink()

