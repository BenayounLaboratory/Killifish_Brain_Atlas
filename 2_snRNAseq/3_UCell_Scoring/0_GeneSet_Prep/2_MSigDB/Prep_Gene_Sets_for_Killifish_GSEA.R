setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways/')
options(stringsAsFactors = F)

# 2024-02-21
# prep MSigDB mouse gmts for fast loading for GSEA with phenotest
# "Killify" gene sets

library(phenoTest)
library(qusage)

# Read in BLAST homology file for killifish/mouse (best mouse hit to killifish to get conversion)
mouse.homol <- read.csv("../../Mouse_alignment/2022-10-11_Mouse_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-3_Minimal_HOMOLOGY_TABLE_REV.txt", sep = "\t", header = T)

#########

Sym.m2.biocarta  <- read.gmt("/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways/m2.cp.biocarta.v2023.2.Mm.symbols.gmt")
Sym.m2.wiki      <- read.gmt("/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways/m2.cp.wikipathways.v2023.2.Mm.symbols.gmt")
Sym.m2.reactome  <- read.gmt("/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways/m2.cp.reactome.v2023.2.Mm.symbols.gmt")
Sym.m5.gobp      <- read.gmt("/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways/m5.go.bp.v2023.2.Mm.symbols.gmt")
Sym.m5.goall     <- read.gmt("/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways/m5.go.v2023.2.Mm.symbols.gmt")
Sym.mh.hallmarks <- read.gmt("/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways/mh.all.v2023.2.Mm.symbols.gmt")


### function for conversion
killify_gsets <- function(my.gset.list, homol.table = mouse.homol) {
  
  # prepare list
  my.gset.list.killi <- vector(length = length(my.gset.list), mode = "list")
  names(my.gset.list.killi) <- names(my.gset.list) 
  
  # grab all killi homologs to mouse genes in gene set
  for (i in 1:length(my.gset.list)) {
    my.gset.list.killi[[i]] <- unique(homol.table$Nfur_Symbol[homol.table$Mmu_Symbol %in% my.gset.list[[i]]])
  }
  return(my.gset.list.killi)
}

# Killify gene sets
Killi.m2.biocarta  <- killify_gsets(Sym.m2.biocarta )
Killi.m2.wiki      <- killify_gsets(Sym.m2.wiki     )
Killi.m2.reactome  <- killify_gsets(Sym.m2.reactome )
Killi.m5.gobp      <- killify_gsets(Sym.m5.gobp     )
Killi.m5.goall     <- killify_gsets(Sym.m5.goall    )
Killi.mh.hallmarks <- killify_gsets(Sym.mh.hallmarks)



######### Clean up GMT pathway names
dataup <- function(x, y) {
  substr(x, 1, y) <- toupper(substr(x, 1, y))
  x
}

clean_pathname <- function(my.coll, offset) {
  names(my.coll) <- dataup( tolower(names(my.coll)), offset )
  return(my.coll)
}


Killi.m2.biocarta  <- clean_pathname(Killi.m2.biocarta , 10 )
Killi.m2.wiki      <- clean_pathname(Killi.m2.wiki     , 4)
Killi.m2.reactome  <- clean_pathname(Killi.m2.reactome , 10)
Killi.m5.gobp      <- clean_pathname(Killi.m5.gobp     , 6)
Killi.m5.goall     <- clean_pathname(Killi.m5.goall    , 6)
Killi.mh.hallmarks <- clean_pathname(Killi.mh.hallmarks, 10)




save(Killi.m2.biocarta  ,
     Killi.m2.wiki      ,
     Killi.m2.reactome  ,
     Killi.m5.gobp      ,
     Killi.m5.goall     ,
     Killi.mh.hallmarks ,
     file = paste0(Sys.Date(),"_Killified_GeneSetCollections_for_Phenotest_GSEA.RData"))
