setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways/')
options(stringsAsFactors = F)

# 2022-12-23
# prep MSigDB mouse gmts for fast loading for clusterprofiler ORA
# "Killify" gene sets

library(clusterProfiler)

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
# my.gset.tab <- Sym.m2.biocarta
# homol.table <- mouse.homol

killify_gsets <- function(my.gset.tab, homol.table = mouse.homol) {
  
  # prepare dataframe
  my.gset.tab.killi <- data.frame(colnames = colnames(my.gset.tab))
  my.gset.tab.killi <- unique(merge(my.gset.tab, homol.table[,3:4], by.x = "gene", by.y = "Mmu_Symbol"))
  my.gset.tab.killi <- my.gset.tab.killi[, c("term",  "Nfur_Symbol")]
  colnames(my.gset.tab.killi) <- c("term",  "gene")
  return(my.gset.tab.killi)
}

# Killify gene sets
Killi.m2.biocarta  <- killify_gsets(Sym.m2.biocarta )
Killi.m2.wiki      <- killify_gsets(Sym.m2.wiki     )
Killi.m2.reactome  <- killify_gsets(Sym.m2.reactome )
Killi.m5.gobp      <- killify_gsets(Sym.m5.gobp     )
Killi.m5.goall     <- killify_gsets(Sym.m5.goall    )
Killi.mh.hallmarks <- killify_gsets(Sym.mh.hallmarks)

# order
Killi.m2.biocarta  <-  Killi.m2.biocarta [order(Killi.m2.biocarta $term),]
Killi.m2.wiki      <-  Killi.m2.wiki     [order(Killi.m2.wiki     $term),]
Killi.m2.reactome  <-  Killi.m2.reactome [order(Killi.m2.reactome $term),]
Killi.m5.gobp      <-  Killi.m5.gobp     [order(Killi.m5.gobp     $term),]
Killi.m5.goall     <-  Killi.m5.goall    [order(Killi.m5.goall    $term),]
Killi.mh.hallmarks <-  Killi.mh.hallmarks[order(Killi.mh.hallmarks$term),]



######### Clean up GMT pathway names
dataup <- function(x, y) {
  substr(x, 1, y) <- toupper(substr(x, 1, y))
  x
}

clean_pathname <- function(my.coll, offset) {
  my.coll$term <- dataup( tolower(my.coll$term), offset )
  return(my.coll)
}


Killi.m2.biocarta  <- clean_pathname(Killi.m2.biocarta , 10 )
Killi.m2.wiki      <- clean_pathname(Killi.m2.wiki     , 4)
Killi.m2.reactome  <- clean_pathname(Killi.m2.reactome , 10)
Killi.m5.gobp      <- clean_pathname(Killi.m5.gobp     , 6)
Killi.m5.goall     <- clean_pathname(Killi.m5.goall    , 6)
Killi.mh.hallmarks <- clean_pathname(Killi.mh.hallmarks, 10)


# save for export
save(Killi.m2.biocarta  ,
     Killi.m2.wiki      ,
     Killi.m2.reactome  ,
     Killi.m5.gobp      ,
     Killi.m5.goall     ,
     Killi.mh.hallmarks ,
     file = paste0(Sys.Date(),"_Killified_GeneSetCollections_for_ClusterProfilerORA.RData"))
