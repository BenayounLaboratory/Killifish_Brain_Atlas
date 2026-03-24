setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways_GeneSets/ENSEMBL')
options(stringsAsFactors = F)

# 2022-12-23
# prep MSigDB mouse gmts for fast loading for clusterprofiler ORA
# "Killify" gene sets

library(clusterProfiler)

# Read in BLAST homology file for killifish/mouse (best mouse hit to killifish to get conversion)
zebra.homol <- read.csv("/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Zebrafish_alignment/2022-03-15_Zebrafish_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-5_Minimal_HOMOLOGY_TABLE_REV.txt", sep = "\t", header = T)

ENS.GO_ALL <- read.gmt("2024-03-20_Zebrafish_Ens111_GO_ALL.gmt")
ENS.GO_BP  <- read.gmt("2024-03-20_Zebrafish_Ens111_GO_BP.gmt")
ENS.GO_CC  <- read.gmt("2024-03-20_Zebrafish_Ens111_GO_CC.gmt")
ENS.GO_MF  <- read.gmt("2024-03-20_Zebrafish_Ens111_GO_MF.gmt")


### function for conversion
# my.gset.tab <- Sym.m2.biocarta
# homol.table <- mouse.homol

killify_gsets <- function(my.gset.tab, homol.table = zebra.homol) {
  
  # prepare dataframe
  my.gset.tab.killi <- data.frame(colnames = colnames(my.gset.tab))
  my.gset.tab.killi <- unique(merge(my.gset.tab, homol.table[,3:4], by.x = "gene", by.y = "DanRer_Symbol"))
  my.gset.tab.killi <- my.gset.tab.killi[, c("term",  "Nfur_Symbol")]
  colnames(my.gset.tab.killi) <- c("term",  "gene")
  return(my.gset.tab.killi)
}

# Killify gene sets
Killi.ENS.GO_ALL  <- killify_gsets(ENS.GO_ALL)
Killi.ENS.GO_BP   <- killify_gsets(ENS.GO_BP )
Killi.ENS.GO_CC   <- killify_gsets(ENS.GO_CC )
Killi.ENS.GO_MF   <- killify_gsets(ENS.GO_MF )

# order
Killi.ENS.GO_ALL  <-  unique(Killi.ENS.GO_ALL[order(Killi.ENS.GO_ALL$term),])
Killi.ENS.GO_BP   <-  unique(Killi.ENS.GO_BP [order(Killi.ENS.GO_BP $term),])
Killi.ENS.GO_CC   <-  unique(Killi.ENS.GO_CC [order(Killi.ENS.GO_CC $term),])
Killi.ENS.GO_MF   <-  unique(Killi.ENS.GO_MF [order(Killi.ENS.GO_MF $term),])


# save for export
save(Killi.ENS.GO_ALL ,
     Killi.ENS.GO_BP  ,
     Killi.ENS.GO_CC  ,
     Killi.ENS.GO_MF  ,
     file = paste0(Sys.Date(),"_Killified_GeneSetCollections_for_ClusterProfilerORA.RData"))
