setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/Integration_snRNA_snATAC/mitch')
options(stringsAsFactors = F)

library("mitch")

# 2024-07-11
# Integrate data from snRNA and snATAC

# 2025-06-12
# rerun with updated GO terms

# 2025-06-23
# rerun with GRZ/ZMZ snRNA, GRZ anATAC for each cell type 

################################################################################
# 1. load DEseq2 results across cell types and 'omes'

# load GRZ aging PB snRNA analysis
load('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Differential_Expression_Age/DGE_Analysis/2024-02-28_pseudobulk_killi_cell_types_AGING_GRZ_DEseq2_objects.RData')
deseq.res.list.genes.grz

# rename for ease
rna.GRZ.deseq2 <- deseq.res.list.genes.grz; rm(deseq.res.list.genes.grz)

head(rna.GRZ.deseq2$Microglia)
#               baseMean log2FoldChange      lfcSE        stat      pvalue       padj
# LOC107382895 0.9133813   -0.135766141 0.12181910 -1.11448980 0.265069120 0.46973867
# LOC107382813 1.3085666   -0.053706546 0.09692577 -0.55409975 0.579510584 0.74964186
# fmr1         3.3062475   -0.018191648 0.05770851 -0.31523338 0.752584461 0.86329171
# aff2         9.4305222    0.136981151 0.05005389  2.73667354 0.006206385 0.03863104
# ids          4.3389202   -0.056173761 0.06102009 -0.92057816 0.357270710 0.56321768
# dcps         1.2668458    0.008089169 0.09022930  0.08965125 0.928564358 0.96288751


# load ZMZ aging PB snRNA analysis
load('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Differential_Expression_Age/DGE_Analysis/2024-02-28_pseudobulk_killi_cell_types_AGING_ZMZ_DEseq2_objects.RData')
deseq.res.list.genes.zmz

# rename for ease
rna.ZMZ.deseq2 <- deseq.res.list.genes.zmz; rm(deseq.res.list.genes.zmz)

head(rna.ZMZ.deseq2$Microglia)
#               baseMean log2FoldChange      lfcSE       stat    pvalue      padj
# LOC107382813  1.854668   -0.017391532 0.04073958 -0.4268952 0.6694557 0.8630634
# fmr1          3.933734   -0.015606817 0.02525646 -0.6179338 0.5366190 0.7883124
# aff2         12.405948   -0.029053223 0.01860607 -1.5614921 0.1184077 0.3773088
# ids           4.681040   -0.030354460 0.02350228 -1.2915539 0.1965117 0.4869655
# dcps          2.104653    0.025420087 0.03392973  0.7491980 0.4537379 0.7283313
# LOC107382040  1.706961   -0.009504424 0.04087406 -0.2325295 0.8161268 0.9309700


# load parsed GRZ aging PB snATAC
load('/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/snATAC_Brain_Aging_Meta/Downstream_Analyses/GSEA/2024-06-24_Annotated_Parsed_snATAC_DESeq2_aging_results.RData')
parsed.deseq.res.list.grz.annot

# rename for ease
atac.deseq2 <- parsed.deseq.res.list.grz.annot; rm(parsed.deseq.res.list.grz.annot)

head(atac.deseq2$Microglia)
# Row.names baseMean log2FoldChange      lfcSE        stat    pvalue      padj Gene.Name Distance.to.TSS Genomic_Context
# 1 NC_029649.1-10373388-10374425 24.13939   -0.022899362 0.02796488 -0.81886135 0.4128655 0.9107354     hvcn1             444            exon
# 2 NC_029649.1-10411040-10411935 20.14214    0.053199557 0.03264233  1.62977200 0.1031497 0.6889504     tctn1             373            exon
# 4   NC_029649.1-1056907-1057846 14.50231   -0.002294323 0.03822776 -0.06001719 0.9521419 0.9934442   fam199x             190            exon
# 5   NC_029649.1-1082161-1083046 19.31279   -0.047959990 0.03132300 -1.53114310 0.1257340 0.7203879      pigg              90    promoter-TSS
# 6 NC_029649.1-10872368-10873250 11.02426    0.024149315 0.04392091  0.54983634 0.5824316 0.9496952     uvssa              16    promoter-TSS
# 7 NC_029649.1-10882184-10883095 12.85762    0.060731364 0.04393368  1.38234183 0.1668668 0.7709106      sil1              95    promoter-TSS

# load prepped gene sets
load('../../scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways_GeneSets/ENSEMBL/2025-01-17_Killified_ENSEMBL_ZebraGO_GeneSetCollections_for_Phenotest_GSEA.RData')
load('../../scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways_GeneSets/MsigDB/2024-02-23_Killified_GeneSetCollections_for_Phenotest_GSEA.RData')


################################################################################
# 2. Import data for mitch

# restrict to cell types seen in both 'omes'
cell.types.mitch <- intersect(names(rna.GRZ.deseq2), names(atac.deseq2))

rna.GRZ.deseq2  <- rna.GRZ.deseq2[cell.types.mitch]
rna.ZMZ.deseq2  <- rna.ZMZ.deseq2[cell.types.mitch]
atac.deseq2     <- atac.deseq2   [cell.types.mitch]

# create structure for the mitch data import
mitch.data.grz <- vector(mode = "list", length = length(cell.types.mitch) )
names(mitch.data.grz) <- cell.types.mitch

# create structure for the mitch results
mitch.res.grz.goall <- vector(mode = "list", length = length(cell.types.mitch) )
mitch.res.grz.react <- vector(mode = "list", length = length(cell.types.mitch) )
names(mitch.res.grz.goall) <- cell.types.mitch
names(mitch.res.grz.react) <- cell.types.mitch

for (i in 1:length(cell.types.mitch) ) {
  
  # reformat ATAC to help matching
  atac.tmp <- atac.deseq2[cell.types.mitch[i]][[1]]
  rownames(atac.tmp) <- atac.tmp$Gene.Name
  atac.tmp <- atac.tmp[,c("baseMean", "log2FoldChange","lfcSE","stat","pvalue","padj")]
  
  # create required list
  tmp.mitch.list <- list("RNA_GRZ"  = rna.GRZ.deseq2[cell.types.mitch[i]][[1]],
                         "RNA_ZMZ"  = rna.ZMZ.deseq2[cell.types.mitch[i]][[1]],
                         "ATAC_GRZ" = atac.tmp)
  
  # generate mitch input
  mitch.data.grz[[i]] <- mitch_import(tmp.mitch.list, DEtype="DESeq2")
  
  # prioritisation by significance
  mitch.res.grz.goall[[i]] <- mitch_calc(mitch.data.grz[[i]], Killi.ENS.GO_ALL , priority = "significance", minsetsize = 10, resrows = 20, cores = 1)
  mitch.res.grz.react[[i]] <- mitch_calc(mitch.data.grz[[i]], Killi.m2.reactome, priority = "significance", minsetsize = 10, resrows = 20, cores = 1)

  # head(mitch.res.grz.react[[i]]$enrichment_result)
  
  mitch_report(mitch.res.grz.goall[[i]] , paste(Sys.Date(),cell.types.mitch[i],"GO_ALL_mitch_report.html", sep = "_"))
  mitch_plots(mitch.res.grz.goall[[i]]  , outfile = paste(Sys.Date(),cell.types.mitch[i],"GO_ALL_mitch_charts.pdf", sep = "_"))

  mitch_report(mitch.res.grz.react[[i]] , paste(Sys.Date(),cell.types.mitch[i],"REACTOME_mitch_report.html", sep = "_"))
  mitch_plots(mitch.res.grz.react[[i]]  , outfile = paste(Sys.Date(),cell.types.mitch[i],"REACTOME_mitch_charts.pdf", sep = "_"))

}

save(mitch.data.grz,
     mitch.res.grz.goall, 
     mitch.res.grz.react,
     file = paste0(Sys.Date(),"_GRZ_ZMZ_mitch_snRNA_snATAC_Brain_Aging_Results.RData"))
###############################################################################################


###############################################################################################
# 3. Export to xlsx
options(java.parameters = "-Xmx16g" )
require(openxlsx)

get_enrich_res <- function (mitch.res) {
  mitch.res$enrichment_result[mitch.res$enrichment_result$p.adjustMANOVA < 0.05,]
}


##
GOALL.out <-  paste0(Sys.Date(),"_GRZ_GOALL_mitch_snRNA_snATAC_Brain_Aging_Results_FDR5.xlsx")
write.xlsx(lapply(mitch.res.grz.goall,get_enrich_res), rowNames = F, file = GOALL.out)

##
REACT.out <-  paste0(Sys.Date(),"_GRZ_REACTOME_mitch_snRNA_snATAC_Brain_Aging_Results_FDR5.xlsx")
write.xlsx(lapply(mitch.res.grz.react,get_enrich_res), rowNames = F, file = REACT.out)
###############################################################################################

#######################
sink(file = paste(Sys.Date(),"_MuscatDEseq2_PB_DESeq2_GSEA_scRNAseq_BrainAtlas_session_Info.txt", sep =""))
sessionInfo()
sink()
