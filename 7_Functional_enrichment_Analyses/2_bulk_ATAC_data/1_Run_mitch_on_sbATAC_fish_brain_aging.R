setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/ATAC/Bulk_ATAC_Analysis_Folder/Enrichment/mitch')
options(stringsAsFactors = F)

library("mitch")

library(bitops)
library(ggplot2)          # 
library(scales)           # 
library(ComplexHeatmap)   #
library(circlize)         #

theme_set(theme_bw())   


# 2025-06-26
# run with GRZ/ZMZ bATAC

################################################################################
# 1. load DEseq2 results across cell types and 'omes'

## grab annoated DESeq2 results
res.age.grz <- read.table("../../DESeq2/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_all_genes_statistics_PeakAnnot.txt", sep = "\t",header = T)
res.age.zmz <- read.table("../../DESeq2/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_all_genes_statistics_PeakAnnot.txt", sep = "\t",header = T)

deseq.res.list <- list("GRZ" = res.age.grz,
                       "ZMZ" = res.age.zmz)


# parse to keep only the peak closest to gene
parsed.deseq.res.list.annot        <- vector(mode = "list", length = length(deseq.res.list))
names(parsed.deseq.res.list.annot) <- names(deseq.res.list)

for (i in 1:length(parsed.deseq.res.list.annot)) {
  
  # grab genes
  my.genes <- unique(deseq.res.list$Gene.Name)
  my.keep <- c()
  
  # keep only closest peak for each gene
  for (j in 1:length(my.genes)) {
    my.rows    <- which(deseq.res.list$Gene.Name %in% my.genes[j])
    my.closest <- which.min(abs(deseq.res.list$Distance.to.TSS[my.rows]))
    my.keep <- c(my.keep,my.rows[my.closest])
  }
  
  # select closest
  parsed.deseq.res.list.annot <- deseq.res.list[my.keep,]
  
  # filter anything further than 10kb to the TSS
  parsed.deseq.res.list.annot <- parsed.deseq.res.list.annot[abs(parsed.deseq.res.list.annot$Distance.to.TSS) < 10000,]
  
}

rownames(parsed.deseq.res.list.annot$GRZ) <- parsed.deseq.res.list.annot$GRZ$Gene.Name
rownames(parsed.deseq.res.list.annot$ZMZ) <- parsed.deseq.res.list.annot$ZMZ$Gene.Name

save(parsed.deseq.res.list.annot, file = paste0(Sys.Date(),"_Annotated_Parsed_bATAC_DESeq2_aging_results.RData"))

# load prepped gene sets
load('../../../../scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways_GeneSets/ENSEMBL/2025-01-17_Killified_ENSEMBL_ZebraGO_GeneSetCollections_for_Phenotest_GSEA.RData')
load('../../../../scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways_GeneSets/MsigDB/2024-02-23_Killified_GeneSetCollections_for_Phenotest_GSEA.RData')
################################################################################

################################################################################
# 2. Import data for mitch

# create structure for the mitch data import
mitch.data.atac <- vector(mode = "list", length = 2 )
names(mitch.data.atac) <- c("GRZ","ZMZ")

# reformat ATAC to help matching
parsed.deseq.res.list.annot$GRZ <- parsed.deseq.res.list.annot$GRZ[,c("baseMean", "log2FoldChange","lfcSE","stat","pvalue","padj")]
parsed.deseq.res.list.annot$ZMZ <- parsed.deseq.res.list.annot$ZMZ[,c("baseMean", "log2FoldChange","lfcSE","stat","pvalue","padj")]

# create required list
mitch.list <- list("bATAC_GRZ"  = parsed.deseq.res.list.annot$GRZ,
                   "bATAC_ZMZ"  = parsed.deseq.res.list.annot$ZMZ)

# generate mitch input
mitch.data.atac <- mitch_import(mitch.list, DEtype="DESeq2")

# prioritisation by significance
mitch.res.grz.goall <- mitch_calc(mitch.data.atac, Killi.ENS.GO_ALL , priority = "significance", minsetsize = 10, resrows = 20, cores = 1)
mitch.res.grz.react <- mitch_calc(mitch.data.atac, Killi.m2.reactome, priority = "significance", minsetsize = 10, resrows = 20, cores = 1)

# head(mitch.res.grz.react$enrichment_result)

mitch_report(mitch.res.grz.goall , paste(Sys.Date(),"bulk_ATAC_GO_ALL_mitch_report.html", sep = "_"))
mitch_plots(mitch.res.grz.goall  , outfile = paste(Sys.Date(),"bulk_ATAC_GO_ALL_mitch_charts.pdf", sep = "_"))

mitch_report(mitch.res.grz.react , paste(Sys.Date(),"bulk_ATAC_REACTOME_mitch_report.html", sep = "_"))
mitch_plots(mitch.res.grz.react  , outfile = paste(Sys.Date(),"bulk_ATAC_REACTOME_mitch_charts.pdf", sep = "_"))

save(mitch.data.atac,
     mitch.res.grz.goall, 
     mitch.res.grz.react,
     file = paste0(Sys.Date(),"_GRZ_ZMZ_mitch_bulkATAC_Brain_Aging_Results.RData"))
###############################################################################################


###############################################################################################
# 3. Export to xlsx
options(java.parameters = "-Xmx16g" )
require(openxlsx)

##
GOALL.out <-  paste0(Sys.Date(),"_GRZ_GOALL_mitch_snRNA_snATAC_Brain_Aging_Results_FDR5.xlsx")
write.xlsx(mitch.res.grz.goall$enrichment_result[mitch.res.grz.goall$enrichment_result$p.adjustMANOVA < 0.05,], rowNames = F, file = GOALL.out)

##
REACT.out <-  paste0(Sys.Date(),"_GRZ_REACTOME_mitch_snRNA_snATAC_Brain_Aging_Results_FDR5.xlsx")
write.xlsx(mitch.res.grz.react$enrichment_result[mitch.res.grz.react$enrichment_result$p.adjustMANOVA < 0.05,], rowNames = F, file = REACT.out)
###############################################################################################


###############################################################################################
# 4. load mitch results across cell types

# head(mitch.res.grz.react$enrichment_result)

# filter mitch results to keep consistent signs only
filter_consistent <- function (mitch_enrich) {
  # check if all the same sign
  mitch_enrich$Agreement <- ifelse(abs(apply(sign(mitch_enrich[,c("s.bATAC_GRZ","s.bATAC_ZMZ")]),1,sum)) == 2, "CONSISTENT", "INCONSISTENT")
  
  #filter for consistency across omes and significance at FDR < 5%
  mitch_enrich.filt <- mitch_enrich[bitAnd(mitch_enrich$Agreement %in% "CONSISTENT", mitch_enrich$p.adjustMANOVA < 0.05)>0, ]
  
  # sort on significance
  mitch_enrich.filt <- mitch_enrich.filt[order(mitch_enrich.filt$p.adjustMANOVA),]
  rownames(mitch_enrich.filt) <- mitch_enrich.filt$set
  
  return(mitch_enrich.filt)
}

# filter and sort results
mitch.goall.flt <- filter_consistent(mitch.res.grz.goall$enrichment_result)
mitch.react.flt <- filter_consistent(mitch.res.grz.react$enrichment_result)

#### plots
# get top sig and sort on effect
my.top.10 <- as.matrix(mitch.goall.flt[1:10,c("s.bATAC_GRZ","s.bATAC_ZMZ")])
med_eff <- apply(my.top.10,1,median)
my.top.10 <- my.top.10[sort(med_eff, decreasing = T, index.return = T)$ix,]

go.all.plot <- Heatmap(my.top.10, 
                       col = colorRamp2(c(-0.75,0, 0.75), c("darkblue", "white", "#CC3333"), transparency = 0, space = "LAB"),
                       border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
                       cluster_rows = F,
                       cluster_columns = F,
                       column_title = "Brain",
                       width  = 2*unit(7, "mm"),
                       height = nrow(my.top.10)*unit(7, "mm") )


pdf(paste0(Sys.Date(),"_Brain_bATAC_Top10_GOALL_Mitch_FDR5.pdf"), width = 10, height = 5)
go.all.plot
dev.off()



# get top sig and sort on effect
my.top.10 <- as.matrix(mitch.react.flt[1:10,c("s.bATAC_GRZ","s.bATAC_ZMZ")])
med_eff <- apply(my.top.10,1,median)
my.top.10 <- my.top.10[sort(med_eff, decreasing = T, index.return = T)$ix,]

react.plot <- Heatmap(my.top.10, 
                       col = colorRamp2(c(-0.75,0, 0.75), c("darkblue", "white", "#CC3333"), transparency = 0, space = "LAB"),
                       border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
                       cluster_rows = F,
                       cluster_columns = F,
                       column_title = "Brain",
                       width  = 2*unit(7, "mm"),
                       height = nrow(my.top.10)*unit(7, "mm") )


pdf(paste0(Sys.Date(),"_Brain_bATAC_Top10_REACTOME_Mitch_FDR5.pdf"), width = 10, height = 5)
react.plot
dev.off()
###############################################################################################


#######################
sink(file = paste(Sys.Date(),"_mitch_bATAC_scRNAseq_BrainAtlas_session_Info.txt", sep =""))
sessionInfo()
sink()
