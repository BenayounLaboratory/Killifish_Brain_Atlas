setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/GR_signaling/Mifepristone/bulk_Brain_RNAseq/mitch')
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

# 2025-07-17
# plot top heatmap for each comparison

################################################################################
# 1. load mitch results across cell types

# load mitch results
load('2025-07-16_GRZ_MIF_2Cohorts_mitch_bRNA_Brain_Aging_Results.RData')

ls()

# "mitch.data.mif"      "mitch.res.MIF.goall" "mitch.res.MIF.react"


get_enrich_res <- function (mitch.res) {
  mitch.res$enrichment_result
}

mitch.goall <- get_enrich_res(mitch.res.MIF.goall)
mitch.react <- get_enrich_res(mitch.res.MIF.react)

###############################################################################################

################################################################################
# 2. load mitch results across cell types

# head(mitch.react)

# filter mitch results to keep consistent signs only
filter_consistent <- function (mitch_enrich) {
  # check if all the same sign
  mitch_enrich$Agreement_AGE <- ifelse(abs(apply(sign(mitch_enrich[,c("s.AGE_F","s.AGE_M")]),1,sum)) == 2, "CONSISTENT", "INCONSISTENT")
  mitch_enrich$Agreement_MIF <- ifelse(abs(apply(sign(mitch_enrich[,c("s.MIF_F","s.MIF_M")]),1,sum)) == 2, "CONSISTENT", "INCONSISTENT")
  
  #filter for consistency across groups and significance at FDR < 5%
  mitch_enrich.sig  <- mitch_enrich[mitch_enrich$p.adjustMANOVA < 0.05,]
  mitch_enrich.filt <- mitch_enrich[bitAnd(mitch_enrich.sig$Agreement_AGE %in% "CONSISTENT", mitch_enrich.sig$Agreement_MIF %in% "CONSISTENT")>0, ]
  
  # sort on significance
  mitch_enrich.filt <- mitch_enrich.filt[order(mitch_enrich.filt$p.adjustMANOVA),]
  rownames(mitch_enrich.filt) <- mitch_enrich.filt$set
  
  return(mitch_enrich.filt)
}

# filter and sort results
mitch.goall.flt <- filter_consistent(mitch.goall)
mitch.react.flt <- filter_consistent(mitch.react)

#### plots
# get top sig and sort on effect
my.top.6 <- as.matrix(mitch.goall.flt[1:6,c("s.AGE_F","s.AGE_M","s.MIF_F","s.MIF_M")])
med_eff  <- apply(my.top.6,1,median)
my.top.6 <- my.top.6[sort(med_eff, decreasing = T, index.return = T)$ix,]

go.all.plot <- Heatmap(my.top.6, 
                       col = colorRamp2(c(-0.75,0, 0.75), c("darkblue", "white", "#CC3333"), transparency = 0, space = "LAB"),
                       border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
                       cluster_rows = F,
                       cluster_columns = F,
                       column_title = "GRZ BRAIN",
                       width  = 3*unit(7, "mm"),
                       height = nrow(my.top.6)*unit(7, "mm") )
go.all.plot

pdf(paste0(Sys.Date(),"_GRZ_Brain_MIF_Top6_GOALL_Mitch_FDR5.pdf"), width = 10, height = 5)
print(go.all.plot )
dev.off()


# get top sig and sort on effect
my.top.6 <- as.matrix(mitch.react.flt[1:6,c("s.AGE_F","s.AGE_M","s.MIF_F","s.MIF_M")])
med_eff  <- apply(my.top.6,1,median)
my.top.6 <- my.top.6[sort(med_eff, decreasing = T, index.return = T)$ix,]

react.plot <- Heatmap(my.top.6, 
                       col = colorRamp2(c(-0.75,0, 0.75), c("darkblue", "white", "#CC3333"), transparency = 0, space = "LAB"),
                       border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
                       cluster_rows = F,
                       cluster_columns = F,
                       column_title = "GRZ BRAIN",
                       width  = 3*unit(7, "mm"),
                       height = nrow(my.top.6)*unit(7, "mm") )
react.plot

pdf(paste0(Sys.Date(),"_GRZ_Brain_MIF_Top6_REACTOME_Mitch_FDR5.pdf"), width = 10, height = 5)
print(react.plot )
dev.off()
###############################################################################################

#######################
sink(file = paste(Sys.Date(),"_TOP_Mitch_Plotting_scRNAseq_BrainAtlas_session_Info.txt", sep =""))
sessionInfo()
sink()
