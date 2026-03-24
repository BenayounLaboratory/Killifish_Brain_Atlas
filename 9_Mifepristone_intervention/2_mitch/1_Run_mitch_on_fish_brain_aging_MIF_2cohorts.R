setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/GR_signaling/Mifepristone/bulk_Brain_RNAseq/mitch')
options(stringsAsFactors = F)

library("mitch")

# 2025-07-16
# Integrate data from aging and MIF (female and male)


################################################################################
# 1. load DEseq2 results across cell types and 'omes'

load('../DESeq2/2025-07-15_Brain_MIF_DEseq2_results.RData')


# load prepped gene sets
load('../../../../scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways_GeneSets/ENSEMBL/2025-01-17_Killified_ENSEMBL_ZebraGO_GeneSetCollections_for_Phenotest_GSEA.RData')
load('../../../../scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways_GeneSets/MsigDB/2024-02-23_Killified_GeneSetCollections_for_Phenotest_GSEA.RData')


################################################################################
# 2. Import data for mitch

# create required list
mif.mitch.list <- list("AGE_F"  = F.mif.res$Aging,
                       "AGE_M"  = M.mif.res$Aging,
                       "MIF_F"  = F.mif.res$Mif,
                       "MIF_M"  = M.mif.res$Mif)

# generate mitch input
mitch.data.mif <- mitch_import(mif.mitch.list, DEtype="DESeq2")

# prioritisation by significance
mitch.res.MIF.goall <- mitch_calc(mitch.data.mif, Killi.ENS.GO_ALL , priority = "significance", minsetsize = 10, resrows = 20, cores = 1)
mitch.res.MIF.react <- mitch_calc(mitch.data.mif, Killi.m2.reactome, priority = "significance", minsetsize = 10, resrows = 20, cores = 1)

# head(mitch.res.grz.react[[i]]$enrichment_result)

mitch_report(mitch.res.MIF.goall , paste(Sys.Date(),"MIF_2Cohorts_GO_ALL_mitch_report.html", sep = "_"))
mitch_plots(mitch.res.MIF.goall  , outfile = paste(Sys.Date(),"MIF_2Cohorts_GO_ALL_mitch_charts.pdf", sep = "_"))

mitch_report(mitch.res.MIF.react , paste(Sys.Date(),"MIF_2Cohorts_REACTOME_mitch_report.html", sep = "_"))
mitch_plots(mitch.res.MIF.react  , outfile = paste(Sys.Date(),"MIF_2Cohorts_REACTOME_mitch_charts.pdf", sep = "_"))


save(mitch.data.mif,
     mitch.res.MIF.goall, 
     mitch.res.MIF.react,
     file = paste0(Sys.Date(),"_GRZ_MIF_2Cohorts_mitch_bRNA_Brain_Aging_Results.RData"))
###############################################################################################


###############################################################################################
# 3. Export to xlsx
options(java.parameters = "-Xmx16g" )
require(openxlsx)

##
GOALL.out <-  paste0(Sys.Date(),"_MIF_2Cohorts_GOALL_mitch_bulkRNA_Brain_Aging_Results_FDR5.xlsx")
write.xlsx(mitch.res.MIF.goall$enrichment_result[mitch.res.MIF.goall$enrichment_result$p.adjustMANOVA < 0.05,], rowNames = F, file = GOALL.out)

##
REACT.out <-  paste0(Sys.Date(),"_MIF_2Cohorts_REACTOME_mitch_bulkRNA_Brain_Aging_Results_FDR5.xlsx")
write.xlsx(mitch.res.MIF.react$enrichment_result[mitch.res.MIF.react$enrichment_result$p.adjustMANOVA < 0.05,], rowNames = F, file = REACT.out)
###############################################################################################

#######################
sink(file = paste(Sys.Date(),"_Mitch_MIF_2Cohorts_DESeq2_GSEA_scRNAseq_BrainAtlas_session_Info.txt", sep =""))
sessionInfo()
sink()
