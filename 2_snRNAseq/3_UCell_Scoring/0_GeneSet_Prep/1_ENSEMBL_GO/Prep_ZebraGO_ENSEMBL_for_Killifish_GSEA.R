setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways_GeneSets/ENSEMBL')
options(stringsAsFactors = F)

# 2024-03-20
# prep ENSEMBL GO gmt for fast loading for GSEA with phenotest
# "Killify" gene sets

# 2025-01-17
# save correct output

library(phenoTest)
library(qusage)

# Read in BLAST homology file for killifish/mouse (best mouse hit to killifish to get conversion)
zebra.homol <- read.csv("/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Zebrafish_alignment/2022-03-15_Zebrafish_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-5_Minimal_HOMOLOGY_TABLE_REV.txt", sep = "\t", header = T)

ENS.GO_ALL <- read.gmt("2024-03-20_Zebrafish_Ens111_GO_ALL.gmt")
ENS.GO_BP  <- read.gmt("2024-03-20_Zebrafish_Ens111_GO_BP.gmt")
ENS.GO_CC  <- read.gmt("2024-03-20_Zebrafish_Ens111_GO_CC.gmt")
ENS.GO_MF  <- read.gmt("2024-03-20_Zebrafish_Ens111_GO_MF.gmt")


### function for conversion
killify_gsets <- function(my.gset.list, homol.table = zebra.homol) {
  
  # prepare list
  my.gset.list.killi <- vector(length = length(my.gset.list), mode = "list")
  names(my.gset.list.killi) <- names(my.gset.list) 
  
  # grab all killi homologs to mouse genes in gene set
  for (i in 1:length(my.gset.list)) {
    my.gset.list.killi[[i]] <- unique(homol.table$Nfur_Symbol[homol.table$DanRer_Symbol %in% my.gset.list[[i]]])
  }
  return(my.gset.list.killi)
}

# Killify gene sets
Killi.ENS.GO_ALL  <- killify_gsets(ENS.GO_ALL)
Killi.ENS.GO_BP   <- killify_gsets(ENS.GO_BP )
Killi.ENS.GO_CC   <- killify_gsets(ENS.GO_CC )
Killi.ENS.GO_MF   <- killify_gsets(ENS.GO_MF )



save(Killi.ENS.GO_ALL ,
     Killi.ENS.GO_BP  ,
     Killi.ENS.GO_CC  ,
     Killi.ENS.GO_MF  ,
     file = paste0(Sys.Date(),"_Killified_ENSEMBL_ZebraGO_GeneSetCollections_for_Phenotest_GSEA.RData"))
