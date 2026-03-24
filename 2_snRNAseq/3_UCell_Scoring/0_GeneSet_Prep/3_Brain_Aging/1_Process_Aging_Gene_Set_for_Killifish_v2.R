setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/Pathways_GeneSets')
options(stringsAsFactors = F)

library(readxl)
library(bitops)

# 2023-03-20
# Killify aging related gene sets


# 2024-04-12
# Clean/simplify output
####################################################################################################################
# 1.Prepare gene lists of interest

# Read in BLAST homology file for killifish/mouse (best mouse hit to killifish to get conversion)
mouse.homol <- read.csv("/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Mouse_alignment/2022-10-11_Mouse_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-3_Minimal_HOMOLOGY_TABLE_REV.txt", sep = "\t", header = T)
human.homol <- read.csv("/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Human_alignment/2024-03-13_Human_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-3_Minimal_HOMOLOGY_TABLE_REV.txt", sep = "\t", header = T)

### read gene list files
senmayo.mouse    <- read_xlsx('Saul-2022-SenMayoPanel.xlsx', sheet = "mouse")                  # mouse
hahn.cag         <- read_xlsx('Hahn-2023-TableS1.xlsx', sheet = 2, skip = 1)                   # mouse

peng.aging.up    <- read_xls('Peng-2021_GTEX_Brain_Aging_Table_2.XLS', sheet = "Aging_OL4_UP_27"  , skip = 1)   # human (at least 4 datasets)
peng.aging.dwn   <- read_xls('Peng-2021_GTEX_Brain_Aging_Table_2.XLS', sheet = "Aging_OL4_DOWN_64", skip = 1)   # human (at least 4 datasets)


# create gene lists
NAMES.mm.aging.glists <- c("SenMayo"    ,
                           "Hahn_CAG"   )

NAMES.hs.aging.glists <- c("GTEx_UP"    ,
                           "GTEx_DWN"   )


aging.glists.mouse        <- vector(mode = "list", length = length(NAMES.mm.aging.glists))
names(aging.glists.mouse) <- NAMES.mm.aging.glists

aging.glists.mouse$SenMayo             <- unique(senmayo.mouse$`Gene(murine)`)
aging.glists.mouse$Hahn_CAG            <- hahn.cag$`Gene Symbol`[hahn.cag$`Signature class` %in% "CAS gene"]

aging.glists.human        <- vector(mode = "list", length = length(NAMES.hs.aging.glists))
names(aging.glists.human) <- NAMES.hs.aging.glists

aging.glists.human$GTEx_UP            <- unique(peng.aging.up$`Gene Symbol`)
aging.glists.human$GTEx_DWN           <- unique(peng.aging.dwn$`Gene Symbol`)



### functions for conversion
killify_mm_gsets <- function(my.gset.list, homol.table = mouse.homol) {
  
  # prepare list
  my.gset.list.killi <- vector(length = length(my.gset.list), mode = "list")
  names(my.gset.list.killi) <- names(my.gset.list) 
  
  # grab all killi homologs to mouse genes in gene set
  for (i in 1:length(my.gset.list)) {
    my.gset.list.killi[[i]] <- unique(homol.table$Nfur_Symbol[homol.table$Mmu_Symbol %in% my.gset.list[[i]]])
  }
  return(my.gset.list.killi)
}

killify_hs_gsets <- function(my.gset.list, homol.table = human.homol) {
  
  # prepare list
  my.gset.list.killi <- vector(length = length(my.gset.list), mode = "list")
  names(my.gset.list.killi) <- names(my.gset.list) 
  
  # grab all killi homologs to mouse genes in gene set
  for (i in 1:length(my.gset.list)) {
    my.gset.list.killi[[i]] <- unique(homol.table$Nfur_Symbol[homol.table$Human_Symbol %in% my.gset.list[[i]]])
  }
  return(my.gset.list.killi)
}

mm.aging.glists.killi <- killify_mm_gsets(aging.glists.mouse)
hs.aging.glists.killi <- killify_hs_gsets(aging.glists.human)

aging.glists.killi <- c(mm.aging.glists.killi,hs.aging.glists.killi)

save(aging.glists.killi,
     file = paste(Sys.Date(),"Aging_Gene_lists_killified.RData",sep = "_"))
####################################################################################################################
