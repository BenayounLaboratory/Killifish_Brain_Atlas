setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/Integration_snRNA_snATAC/mitch')
options(stringsAsFactors = F)

library("mitch")

library(ggplot2)          # 
library(scales)           # 
library("bitops")         # 
library(Vennerable)       # 
library(data.table)       #

library(bitops)

library(ComplexHeatmap)   #
library(circlize)         #

theme_set(theme_bw())   

# 2025-06-25
# plot top heatmap for each cell type

# 2025-06-30
# plot top 6 per cell type

################################################################################
# 1. load mitch results across cell types

# load mitch results
load('2025-06-23_GRZ_ZMZ_mitch_snRNA_snATAC_Brain_Aging_Results.RData')

# mitch.data.grz
# mitch.res.grz.goall
# mitch.res.grz.react


get_enrich_res <- function (mitch.res) {
  mitch.res$enrichment_result
}

mitch.goall <- lapply(mitch.res.grz.goall,get_enrich_res)
mitch.react <- lapply(mitch.res.grz.react,get_enrich_res)

###############################################################################################

################################################################################
# 2. load mitch results across cell types

# head(mitch.react$Astrocytes_Radial_Glia)

# filter mitch results to keep consistent signs only
filter_consistent <- function (mitch_enrich) {
  # check if all the same sign
  mitch_enrich$Agreement <- ifelse(abs(apply(sign(mitch_enrich[,c("s.RNA_GRZ","s.RNA_ZMZ","s.ATAC_GRZ")]),1,sum)) == 3, "CONSISTENT", "INCONSISTENT")
  
  #filter for consistency across omes and significance at FDR < 5%
  mitch_enrich.filt <- mitch_enrich[bitAnd(mitch_enrich$Agreement %in% "CONSISTENT", mitch_enrich$p.adjustMANOVA < 0.05)>0, ]
  
  # sort on significance
  mitch_enrich.filt <- mitch_enrich.filt[order(mitch_enrich.filt$p.adjustMANOVA),]
  rownames(mitch_enrich.filt) <- mitch_enrich.filt$set
  
  return(mitch_enrich.filt)
}

# filter and sort results
mitch.goall.flt <- lapply(mitch.goall,filter_consistent)
mitch.react.flt <- lapply(mitch.react,filter_consistent)

#### plots
go.all.plots <- vector(mode = "list", length = length(mitch.goall.flt))

for (i in 1:length(mitch.goall.flt)) {
  
  # get top sig and sort on effect
  my.top.6 <- as.matrix(mitch.goall.flt[[i]][1:6,c("s.RNA_GRZ","s.RNA_ZMZ","s.ATAC_GRZ")])
  med_eff <- apply(my.top.6,1,median)
  my.top.6 <- my.top.6[sort(med_eff, decreasing = T, index.return = T)$ix,]
  
  go.all.plots[[i]] <- Heatmap(my.top.6, 
                               col = colorRamp2(c(-0.75,0, 0.75), c("darkblue", "white", "#CC3333"), transparency = 0, space = "LAB"),
                               border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
                               cluster_rows = F,
                               cluster_columns = F,
                               column_title = names(mitch.goall.flt)[i],
                               width  = 3*unit(7, "mm"),
                               height = nrow(my.top.6)*unit(7, "mm") )
  
  
  pdf(paste0(Sys.Date(),"_",names(mitch.goall.flt)[i], "_Top6_GOALL_Mitch_FDR5.pdf"), width = 10, height = 5)
  print(go.all.plots[[i]] )
  dev.off()
  
}

###
react.plots <- vector(mode = "list", length = length(mitch.react.flt))

for (i in 1:length(mitch.react.flt)) {
  
  # get top sig and sort on effect
  my.top.6 <- as.matrix(mitch.react.flt[[i]][1:6,c("s.RNA_GRZ","s.RNA_ZMZ","s.ATAC_GRZ")])
  med_eff <- apply(my.top.6,1,median)
  my.top.6 <- my.top.6[sort(med_eff, decreasing = T, index.return = T)$ix,]
  
  react.plots[[i]] <- Heatmap(my.top.6, 
                               col = colorRamp2(c(-0.75,0, 0.75), c("darkblue", "white", "#CC3333"), transparency = 0, space = "LAB"),
                               border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
                               cluster_rows = F,
                               cluster_columns = F,
                               column_title = names(mitch.react.flt)[i],
                               width  = 3*unit(7, "mm"),
                               height = nrow(my.top.6)*unit(7, "mm") )
  
  
  pdf(paste0(Sys.Date(),"_",names(mitch.goall.flt)[i], "_Top6_REACTOME_Mitch_FDR5.pdf"), width = 10, height = 5)
  print(react.plots[[i]] )
  dev.off()
  
}

###############################################################################################

#######################
sink(file = paste(Sys.Date(),"_TOP_Mitch_Plotting_scRNAseq_BrainAtlas_session_Info.txt", sep =""))
sessionInfo()
sink()
