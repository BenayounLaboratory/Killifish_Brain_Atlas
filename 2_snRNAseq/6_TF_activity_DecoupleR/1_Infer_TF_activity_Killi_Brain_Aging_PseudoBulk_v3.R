setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Downstream_Analyses/TF_activity_analysis')
options(stringsAsFactors = F)

library("DESeq2")        #
library("decoupleR")     #
library("OmnipathR")     #

library(dplyr)
library(tibble)
library(ggplot2)
library(beeswarm)
library(ggrepel)

library(bitops)
library(ComplexHeatmap)

theme_set(theme_bw())   

# 2024-03-08
# Try using decoupleR bioconductor package (and mouse TF target DB)
# to infer differential TF activity in brain aging data
# https://www.bioconductor.org/packages/release/bioc/vignettes/decoupleR/inst/doc/tf_sc.html

# 2024-03-21
# Try pseudobulk analysis since:
#    - single cell seems to notbe able to run in reasonable time
#    - PB is likely more robust and taking care of true biological replication anyway
# Keep already computed killified mouse omnipath database
# https://www.bioconductor.org/packages/release/bioc/vignettes/decoupleR/inst/doc/tf_bk.html

# 2024-04-11
# do additional summarizing steps for data reporting

# 2024-11-20
# also shade the data for FDR10

# 2025-07-07
# output all for jitter

#####################################################################################################################
#### 0. Load up annotated DEseq2/VST objects

# Load necessary objects
load('../Differential_Expression_Age/DGE_Analysis/2024-02-28_pseudobulk_killi_cell_types_AGING_GRZ_DEseq2_objects.RData'   )
load('../Differential_Expression_Age/DGE_Analysis/2024-02-28_pseudobulk_killi_cell_types_AGING_GRZ_VST_data_objects.RData' )
load('../Differential_Expression_Age/DGE_Analysis/2024-02-28_pseudobulk_killi_cell_types_AGING_ZMZ_DEseq2_objects.RData'   )
load('../Differential_Expression_Age/DGE_Analysis/2024-02-28_pseudobulk_killi_cell_types_AGING_ZMZ_VST_data_objects.RData' )

ls()
# [1] "deseq.res.list.genes.grz" "deseq.res.list.genes.zmz" "deseq.res.list.grz"       "deseq.res.list.zmz"       "vst.cts.grz"             
# [6] "vst.cts.zmz"    
#####################################################################################################################


#####################################################################################################################
#### 1. Process TF network data from mouse model

######################################################
####### Translate mouse network to killifish   #######
######################################################
# Read in BLAST homology file for killifish/mouse (best mouse hit to killifish to get conversion)
mouse.homol <- read.csv("../../Mouse_alignment/2022-10-11_Mouse_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-3_Minimal_HOMOLOGY_TABLE_REV.txt", sep = "\t", header = T)

# Get mouse CollecTRI network
mouse.net <- get_collectri(organism='mouse', split_complexes=FALSE)
mouse.net

# Translate mouse genes to killi genes
killi.net <- data.frame(matrix(0,0,3))
colnames(killi.net) <- colnames(mouse.net)

for (i in 1:nrow(mouse.net)) {
  
  mm.source <- mouse.net[i,]$source
  mm.target <- mouse.net[i,]$target
  mm.mor    <- mouse.net[i,]$mor
  
  nf.source <- mouse.homol$Nfur_Symbol[toupper(mouse.homol$Mmu_Symbol) %in% mm.source]
  nf.target <- mouse.homol$Nfur_Symbol[toupper(mouse.homol$Mmu_Symbol) %in% mm.target]
  
  if ((length(nf.target) > 0) &&  (length(nf.source) > 0)) { # if at least one target and one source ortholog exist
    
    if (length(nf.source) ==  1) { # if only one TF source
      nf.tmp <- data.frame("source" = rep(nf.source, length(nf.target)),
                           "target" = nf.target,
                           "mor"    = rep(mm.mor, length(mm.mor)))
      
      killi.net <- rbind(killi.net,nf.tmp)
      
    } else if (length(nf.source) >  1) { # if more than one TF source
      
      # get the first one
      nf.tmp <- data.frame("source" = rep(nf.source[1], length(nf.target)),
                           "target" = nf.target,
                           "mor"    = rep(mm.mor, length(nf.target)))
      
      # loop over the other TF orthologs
      for (j in 1:(length(nf.source)-1) ) {
        nf.tmp <- rbind(nf.tmp,
                        data.frame("source" = rep(nf.source[j+1], length(nf.target)),
                                   "target" = nf.target,
                                   "mor"    = rep(mm.mor, length(nf.target))) )
      }
      killi.net <- rbind(killi.net,nf.tmp)
      
    } 
  }
}

# make tibble and remove redundancy
killi.net <- as_tibble(killi.net) 
killi.net <- unique(killi.net)

# some edges are weird now, remove edges with contradicting mor
killi.net$chain <- paste0(killi.net$source, "_",killi.net$target)
killi.net.dup   <- duplicated(killi.net[,-3])

killi.net.cl <- killi.net[!(killi.net$chain %in% killi.net[killi.net.dup,]$chain),]


# save file for future runs
save(killi.net.cl, file = paste0(Sys.Date(),"_killifish_omnipath_object.RData"))
#####################################################################################################################

#####################################################################################################################
#### 2. Calculate TF network data

load('2024-03-08_killifish_omnipath_object.RData')

killi.net.cl
# # A tibble: 43,175 × 4
# source       target         mor chain                    
# <chr>        <chr>        <dbl> <chr>                    
#   1 myc          tert           1 myc_tert                 
# 2 LOC107381469 LOC107382206     1 LOC107381469_LOC107382206
# 3 LOC107388645 LOC107382206     1 LOC107388645_LOC107382206
# 4 smad3        LOC107396464     1 smad3_LOC107396464       
# 5 LOC107389156 LOC107396464     1 LOC107389156_LOC107396464

n_tfs <- 10

#################################################
####### Activity inference   ++++   GRZ   #######
#################################################

# infer pathway activities from the t-values of the DEGs with aging
# for each cell type

decoupler.list.grz <- vector(mode = "list", length = length(deseq.res.list.genes.grz))
names(decoupler.list.grz) <- names(deseq.res.list.genes.grz)

Top_TFs_res.grz <- data.frame(matrix(0,0,7))
colnames(Top_TFs_res.grz) <- c("statistic", "source","condition","score","p_value","Reg_Rank","Cell_Type")

# Loop over DEseq2 results
for (i in 1:length(deseq.res.list.genes.grz)) {
  
  # Run fgsea scoring
  contrast_acts <- run_fgsea(mat     = deseq.res.list.genes.grz[[i]][, 'stat', drop=FALSE],
                             net     = killi.net.cl                                  , 
                             minsize = 5                                             )
  # contrast_acts
  # # A tibble: 1,214 × 5
  # statistic  source       condition  score p_value
  # <chr>      <chr>        <chr>      <dbl>   <dbl>
  # 1 fgsea      LOC107372367 stat       0.263  0.430 
  # 2 norm_fgsea LOC107372367 stat       0.965  0.430 
  # 3 fgsea      LOC107372434 stat       0.329  0.353 
  # 4 norm_fgsea LOC107372434 stat       1.08   0.353 
  # 5 fgsea      LOC107372679 stat       0.427  0.539 
  
  # We select the norm_fgsea activities and then we show changes in activity with aging
  # Filter norm_fgsea
  tf.scores         <- contrast_acts[contrast_acts$statistic == 'norm_fgsea',]
  tf.sig            <- tf.scores[tf.scores$p_value < 0.10,]
  tf.sig$Reg_Rank   <- NA
  tf.sig$Cell_Type  <- names(deseq.res.list.genes.grz)[i]
  
  decoupler.list.grz[[i]] <- tf.scores
  
  # Only keep TF regulons who TF is expressed in cell type
  tf.sig <- tf.sig[tf.sig$source %in% row.names(deseq.res.list.genes.grz[[i]]),]
  
  # Add rank information
  msk <- tf.sig$score > 0
  tf.sig$Reg_Rank[msk ]  <- round(rank(-tf.sig[msk, 'score']))
  tf.sig$Reg_Rank[!msk]  <- round(rank(-abs(tf.sig[!msk, 'score'])))
  
  Top_TFs_res.grz <- rbind(Top_TFs_res.grz,tf.sig)
  
  # Filter top significant TFs in both signs
  tfs <- tf.sig %>%
    arrange(Reg_Rank) %>%
    head(n_tfs) %>%
    pull(source)
  
  f_contrast_acts <- tf.sig %>% filter(source %in% tfs)
  
  # Plot top each direction
  tf.reg.plot <- ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + geom_bar(aes(fill = score), stat = "identity") 
  tf.reg.plot <- tf.reg.plot + scale_fill_gradient2(low = "#333399", mid = "whitesmoke", high = "#CC3333", midpoint = 0, limits = c(-2.5,2.5)) 
  tf.reg.plot <- tf.reg.plot  + theme(axis.title = element_text(face = "bold", size = 12),
                                      axis.text.x = element_text(angle = 45, hjust = 1, size =10),
                                      axis.text.y = element_text(size =10, face= "bold"),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank() ) 
  tf.reg.plot <- tf.reg.plot + xlab("Enriched TF regulons") + ylab("decoupleR score") + theme(plot.title = element_text(hjust = 0.5))
  tf.reg.plot <- tf.reg.plot + ylim(c(-2.5,2.5)) + ggtitle(paste0(names(deseq.res.list.genes.grz)[i]," (GRZ)")) + coord_flip()
  
  pdf(paste0(Sys.Date(),"_BarPLot_decoupleR_fgsea_Top_", n_tfs , "_TF_regulons_", names(deseq.res.list.genes.grz)[i], "Aging_GRZ_FDR10.pdf"), width = 4, height = 5)
  print(tf.reg.plot)
  dev.off()
  
}

write.table(Top_TFs_res.grz, file = paste0(Sys.Date(),"_decoupleR_fgsea_TFregulons_by_CellType_GRZ_FDR10.txt"), sep = "\t", quote = F, row.names = F)

# Reccurrent TF in more than half the cell types (FDR5)
concat <- function (x) {
  paste0(x, collapse = ",")
}

Top_TFs_res.grz.5 <- Top_TFs_res.grz[Top_TFs_res.grz$p_value < 0.05,]

Top_TFs_res.grz.5$sign <- sign(Top_TFs_res.grz.5$score)
my.num.grz   <- aggregate(Top_TFs_res.grz.5$Cell_Type, by = list("TF_regulon" = Top_TFs_res.grz.5$source), FUN = length)
colnames(my.num.grz)[2] <- "N_CellType"
my.signs.grz <- aggregate(Top_TFs_res.grz.5$sign, by = list("TF_regulon" =Top_TFs_res.grz.5$source), FUN = concat)
colnames(my.signs.grz)[2] <- "Signs"
my.sum.grz <- merge(my.num.grz,my.signs.grz)

# recurrent perturbed regulons
my.sum.grz[my.sum.grz$N_CellType > 3,]
#       TF_regulon N_CellType           Signs
# 16         ddit3          5       1,1,1,1,1
# 65  LOC107375066          8 1,1,1,1,1,1,1,1 nr3c1 nuclear receptor subfamily 3, group C, member 1 (glucocorticoid receptor) [ Nothobranchius furzeri (turquoise killifish) ]
# 73  LOC107376749          5       1,1,1,1,1 rela v-rel avian reticuloendotheliosis viral oncogene homolog A [ Nothobranchius furzeri (turquoise killifish) ]
# 200        stat1          5       1,1,1,1,1
# 91  LOC107384049          4         1,1,1,1 LOC107384049 class E basic helix-loop-helix protein 40 [ Nothobranchius furzeri (turquoise killifish) ]
# 139 LOC107396464          4        1,1,1,-1 jun Jun proto-oncogene, AP-1 transcription factor subunit [ Nothobranchius furzeri (turquoise killifish) ]
# 182         rest          4        -1,1,1,1

pdf(paste0(Sys.Date(),"_decoupleR_fgsea_TFregulons_CellTypeSharing_Boxplot_GRZ_FDR5.pdf"), height = 5, width = 3.5)
boxplot(list("GRZ"= my.sum.grz$N_CellType), col = "goldenrod1", main = "DecoupleR Regulon Sharing", ylab = "Number of cell types", outline = F, ylim = c(0,10))
beeswarm(my.sum.grz$N_CellType[my.sum.grz$N_CellType > 3], method = "compactswarm",corral = "gutter", add = T, pch = 16)
text(1.1, 8, "nr3c1"          , pos = 4, cex = 0.5)
text(1.1, 5,"ddit3,rela,stat1", pos = 4, cex = 0.5)
text(1.1, 4,"bhlhe40,jun,rest", pos = 4, cex = 0.5)
dev.off()


#################################################
####### Activity inference   ++++   ZMZ   #######
#################################################


decoupler.list.zmz <- vector(mode = "list", length = length(deseq.res.list.genes.zmz))
names(decoupler.list.zmz) <- names(deseq.res.list.genes.zmz)


Top_TFs_res.zmz <- data.frame(matrix(0,0,7))
colnames(Top_TFs_res.zmz) <- c("statistic", "source","condition","score","p_value","Reg_Rank","Cell_Type")

# Loop over DEseq2 results
for (i in 1:length(deseq.res.list.genes.zmz)) {
  
  # Run fgsea scoring
  contrast_acts <- run_fgsea(mat     = deseq.res.list.genes.zmz[[i]][, 'stat', drop=FALSE],
                             net     = killi.net.cl                                  , 
                             minsize = 5                                             )
  # contrast_acts
  # # A tibble: 1,214 × 5
  # statistic  source       condition  score p_value
  # <chr>      <chr>        <chr>      <dbl>   <dbl>
  # 1 fgsea      LOC107372367 stat       0.263  0.430 
  # 2 norm_fgsea LOC107372367 stat       0.965  0.430 
  # 3 fgsea      LOC107372434 stat       0.329  0.353 
  # 4 norm_fgsea LOC107372434 stat       1.08   0.353 
  # 5 fgsea      LOC107372679 stat       0.427  0.539 
  
  # We select the norm_fgsea activities and then we show changes in activity with aging
  # Filter norm_fgsea
  tf.scores         <- contrast_acts[contrast_acts$statistic == 'norm_fgsea',]
  tf.sig            <- tf.scores[tf.scores$p_value < 0.10,]
  tf.sig$Reg_Rank   <- NA
  tf.sig$Cell_Type  <- names(deseq.res.list.genes.zmz)[i]
  
  decoupler.list.zmz[[i]] <- tf.scores
  
  # Only keep TF regulons who TF is expressed in cell type
  tf.sig <- tf.sig[tf.sig$source %in% row.names(deseq.res.list.genes.zmz[[i]]),]
  
  # Add rank information
  msk <- tf.sig$score > 0
  tf.sig$Reg_Rank[msk ]  <- round(rank(-tf.sig[msk, 'score']))
  tf.sig$Reg_Rank[!msk]  <- round(rank(-abs(tf.sig[!msk, 'score'])))
  
  Top_TFs_res.zmz <- rbind(Top_TFs_res.zmz,tf.sig)
  
  # Filter top significant TFs in both signs
  tfs <- tf.sig %>%
    arrange(Reg_Rank) %>%
    head(n_tfs) %>%
    pull(source)
  
  f_contrast_acts <- tf.sig %>% filter(source %in% tfs)
  
  # Plot top each direction
  tf.reg.plot <- ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + geom_bar(aes(fill = score), stat = "identity") 
  tf.reg.plot <- tf.reg.plot + scale_fill_gradient2(low = "#333399", mid = "whitesmoke", high = "#CC3333", midpoint = 0, limits = c(-3,3)) 
  tf.reg.plot <- tf.reg.plot  + theme(axis.title = element_text(face = "bold", size = 12),
                                      axis.text.x = element_text(angle = 45, hjust = 1, size =10),
                                      axis.text.y = element_text(size =10, face= "bold"),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank() ) 
  tf.reg.plot <- tf.reg.plot + xlab("Enriched TF regulons") + ylab("decoupleR score") + theme(plot.title = element_text(hjust = 0.5))
  tf.reg.plot <- tf.reg.plot + ylim(c(-3,3)) + ggtitle(paste0(names(deseq.res.list.genes.zmz)[i]," (ZMZ)")) + coord_flip()
  
  pdf(paste0(Sys.Date(),"_BarPLot_decoupleR_fgsea_Top_", n_tfs , "_TF_regulons_", names(deseq.res.list.genes.zmz)[i], "Aging_ZMZ_FDR10.pdf"), width = 4, height = 5)
  print(tf.reg.plot)
  dev.off()
  
}

write.table(Top_TFs_res.zmz, file = paste0(Sys.Date(),"_decoupleR_fgsea_TFregulons_by_CellType_ZMZ_FDR10.txt"), sep = "\t", quote = F, row.names = F)

# Reccurrent TF in more than half the cell types
concat <- function (x) {
  paste0(x, collapse = ",")
}

Top_TFs_res.zmz.5 <- Top_TFs_res.zmz[Top_TFs_res.zmz$p_value < 0.05,]

Top_TFs_res.zmz.5$sign <- sign(Top_TFs_res.zmz.5$score)
my.num.zmz   <- aggregate(Top_TFs_res.zmz.5$Cell_Type, by = list("TF_regulon" = Top_TFs_res.zmz.5$source), FUN = length)
colnames(my.num.zmz)[2] <- "N_CellType"
my.signs.zmz <- aggregate(Top_TFs_res.zmz.5$sign, by = list("TF_regulon" =Top_TFs_res.zmz.5$source), FUN = concat)
colnames(my.signs.zmz)[2] <- "Signs"
my.sum.zmz <- merge(my.num.zmz,my.signs.zmz)

# recurrent perturbed regulons
my.sum.zmz[my.sum.zmz$N_CellType > 3,]
#       TF_regulon N_CellType          Signs
# 25          egr1          5 -1,-1,-1,-1,-1
# 62          irf2          5     1,-1,1,1,1
# 63          irf7          4        1,1,1,1
# 68         kmt2b          4     -1,-1,1,-1
# 112 LOC107385678          5      1,1,1,1,1 stat2 signal transducer and activator of transcription 2 [ Nothobranchius furzeri (turquoise killifish) ]
# 113 LOC107385700          5  -1,1,-1,-1,-1 vhl von Hippel-Lindau tumor suppressor [ Nothobranchius furzeri (turquoise killifish) 
# 120 LOC107386825          5     1,-1,1,1,1 irf1b interferon regulatory factor 1b [ Nothobranchius furzeri (turquoise killifish) ]
# 123 LOC107387991          5      1,1,1,1,1 irf3 interferon regulatory factor 3 [ Nothobranchius furzeri (turquoise killifish) ]
# 180        nfkb1          4        1,1,1,1
# 182         nfya          4    -1,-1,-1,-1
# 195         pax6          5 -1,-1,-1,-1,-1
# 230         sox2          4    -1,-1,-1,-1
# 237          srf          4    -1,-1,-1,-1
# 238        stat1          4        1,1,1,1

pdf(paste0(Sys.Date(),"_decoupleR_fgsea_TFregulons_CellTypeSharing_Boxplot_ZMZ_FDR5.pdf"), height = 5, width = 3.5)
boxplot(list("ZMZ"= my.sum.zmz$N_CellType), col = "forestgreen", main = "DecoupleR Regulon Sharing", ylab = "Number of cell types", outline = F, ylim = c(0,10))
beeswarm(my.sum.zmz$N_CellType[my.sum.zmz$N_CellType > 3], method = "compactswarm",corral = "gutter", add = T, pch = 16)
points(1,3,pch = 16, col = "darkgrey")
text(1.1, 5,paste(sort(my.sum.zmz[my.sum.zmz$N_CellType == 5,]$TF_regulon), collapse = ","), pos = 4, cex = 0.5)
text(1.1, 4,paste(sort(my.sum.zmz[my.sum.zmz$N_CellType == 4,]$TF_regulon), collapse = ","), pos = 4, cex = 0.5)
text(1.1, 3,"nr3c1", pos = 4, cex = 0.5, col = "darkgrey") # in astrocytes, microglia, GABAergic_neurons
dev.off()
#####################################################################################################################

save(decoupler.list.grz, decoupler.list.zmz, file = paste0(Sys.Date(), "_DecoupleR_results.RData"))


#####################################################################################################################
Top_TFs_res.grz <- read.table('2024-11-20_decoupleR_fgsea_TFregulons_by_CellType_GRZ_FDR10.txt', header = T)
Top_TFs_res.zmz <- read.table('2024-11-20_decoupleR_fgsea_TFregulons_by_CellType_ZMZ_FDR10.txt', header = T)

# Reccurrent TF in more than half the cell types
concat <- function (x) {
  paste0(x, collapse = ",")
}


# get recurrent data
Top_TFs_res.grz.5 <- Top_TFs_res.grz[Top_TFs_res.grz$p_value < 0.05,]
Top_TFs_res.grz.5$sign <- sign(Top_TFs_res.grz.5$score)
my.num.grz   <- aggregate(Top_TFs_res.grz.5$Cell_Type, by = list("TF_regulon" = Top_TFs_res.grz.5$source), FUN = length)
colnames(my.num.grz)[2] <- "N_CellType"
my.signs.grz <- aggregate(Top_TFs_res.grz.5$sign, by = list("TF_regulon" =Top_TFs_res.grz.5$source), FUN = concat)
colnames(my.signs.grz)[2] <- "Signs"
my.cells.grz <- aggregate(Top_TFs_res.grz.5$Cell_Type, by = list("TF_regulon" =Top_TFs_res.grz.5$source), FUN = concat)
colnames(my.cells.grz)[2] <- "Cell_Types"
my.sum.grz <- merge(merge(my.cells.grz,my.num.grz),my.signs.grz)


Top_TFs_res.zmz.5 <- Top_TFs_res.zmz[Top_TFs_res.zmz$p_value < 0.05,]
Top_TFs_res.zmz.5$sign <- sign(Top_TFs_res.zmz.5$score)
my.num.zmz   <- aggregate(Top_TFs_res.zmz.5$Cell_Type, by = list("TF_regulon" = Top_TFs_res.zmz.5$source), FUN = length)
colnames(my.num.zmz)[2] <- "N_CellType"
my.signs.zmz <- aggregate(Top_TFs_res.zmz.5$sign, by = list("TF_regulon" =Top_TFs_res.zmz.5$source), FUN = concat)
colnames(my.signs.zmz)[2] <- "Signs"
my.cells.zmz <- aggregate(Top_TFs_res.zmz.5$Cell_Type, by = list("TF_regulon" =Top_TFs_res.zmz.5$source), FUN = concat)
colnames(my.cells.zmz)[2] <- "Cell_Types"
my.sum.zmz <- merge(merge(my.cells.zmz,my.num.zmz),my.signs.zmz)


# recurrent perturbed regulons in 4 or more cell types in GRZ or ZMZ data
my.rec.grz <- my.sum.grz[my.sum.grz$N_CellType > 3,]
my.rec.zmz <- my.sum.zmz[my.sum.zmz$N_CellType > 3,]

# grab in the other strain
my.grz.rec_zmz <- my.sum.grz[my.sum.grz$TF_regulon %in% my.rec.zmz$TF_regulon,]
my.zmz.rec_grz <- my.sum.zmz[my.sum.zmz$TF_regulon %in% my.rec.grz$TF_regulon,]

# consolidate
my.rec.grz.v2 <- unique(rbind(my.rec.grz, my.grz.rec_zmz))
my.rec.zmz.v2 <- unique(rbind(my.rec.zmz, my.zmz.rec_grz))

for (i in 1:nrow(my.rec.grz.v2)) {
  my.rec.grz.v2$Direction[i] <- ifelse(length(unique(Top_TFs_res.grz.5[Top_TFs_res.grz.5$source %in% my.rec.grz.v2$TF_regulon[i],]$sign)) == 1,"consistent", "inconsistent")
}

my.rec.grz.v2[,-2]
# TF_regulon N_CellType           Signs    Direction
# 16         ddit3          5       1,1,1,1,1   consistent
# 65  LOC107375066          8 1,1,1,1,1,1,1,1   consistent # nr3c1 nuclear receptor subfamily 3, group C, member 1 (glucocorticoid receptor) [ Nothobranchius furzeri (turquoise killifish) ]
# 73  LOC107376749          5       1,1,1,1,1   consistent # rela v-rel avian reticuloendotheliosis viral oncogene homolog A [ Nothobranchius furzeri (turquoise killifish) ]
# 91  LOC107384049          4         1,1,1,1   consistent # bhlhe40 class E basic helix-loop-helix protein 40 [ Nothobranchius furzeri (turquoise killifish) ]
# 139 LOC107396464          4        1,1,1,-1 inconsistent # jun Jun proto-oncogene, AP-1 transcription factor subunit [ Nothobranchius furzeri (turquoise killifish) ]
# 182         rest          4        -1,1,1,1 inconsistent
# 200        stat1          5       1,1,1,1,1   consistent
# 27          egr1          2           -1,-1   consistent
# 50          irf2          3           1,1,1   consistent
# 52          irf7          1               1   consistent
# 95  LOC107385678          3           1,1,1   consistent # stat2 signal transducer and activator of transcription 2 [ Nothobranchius furzeri (turquoise killifi
# 96  LOC107385700          1               1   consistent # vhl von Hippel-Lindau tumor suppressor [ Nothobranchius furzeri (turquoise killifish) 
# 101 LOC107386825          1               1   consistent # irf1b interferon regulatory factor 1b [ Nothobranchius furzeri (turquoise killifish) ]
# 105 LOC107387991          3           1,1,1   consistent # irf3 interferon regulatory factor 3 [ Nothobranchius furzeri (turquoise killifish) ]
# 153        nfkb1          3           1,1,1   consistent
# 156         nfya          3        -1,-1,-1   consistent
# 166         pax6          2           -1,-1   consistent
# 191         sox2          2           -1,-1   consistent


for (i in 1:nrow(my.rec.zmz.v2)) {
  my.rec.zmz.v2$Direction[i] <- ifelse(length(unique(Top_TFs_res.zmz.5[Top_TFs_res.zmz.5$source %in% my.rec.zmz.v2$TF_regulon[i],]$sign)) == 1,"consistent", "inconsistent")
}

my.rec.zmz.v2[,-2]
# TF_regulon N_CellType          Signs    Direction
# 25          egr1          5 -1,-1,-1,-1,-1   consistent
# 62          irf2          5     1,-1,1,1,1 inconsistent
# 63          irf7          4        1,1,1,1   consistent
# 68         kmt2b          4     -1,-1,1,-1 inconsistent
# 112 LOC107385678          5      1,1,1,1,1   consistent # stat2 signal transducer and activator of transcription 2 [ Nothobranchius furzeri (turquoise killifi
# 113 LOC107385700          5  -1,1,-1,-1,-1 inconsistent # vhl von Hippel-Lindau tumor suppressor [ Nothobranchius furzeri (turquoise killifish) 
# 120 LOC107386825          5     1,-1,1,1,1 inconsistent # irf1b interferon regulatory factor 1b [ Nothobranchius furzeri (turquoise killifish) ]
# 123 LOC107387991          5      1,1,1,1,1   consistent # irf3 interferon regulatory factor 3 [ Nothobranchius furzeri (turquoise killifish) ]
# 180        nfkb1          4        1,1,1,1   consistent
# 182         nfya          4    -1,-1,-1,-1   consistent
# 195         pax6          5 -1,-1,-1,-1,-1   consistent
# 230         sox2          4    -1,-1,-1,-1   consistent
# 237          srf          4    -1,-1,-1,-1   consistent
# 238        stat1          4        1,1,1,1   consistent
# 80  LOC107375066          3          1,1,1   consistent # nr3c1 nuclear receptor subfamily 3, group C, member 1 (glucocorticoid receptor) [ Nothobranchius furzeri (turquoise killifish) ]
# 85  LOC107376749          3          1,1,1   consistent # rela v-rel avian reticuloendotheliosis viral oncogene homolog A [ Nothobranchius furzeri (turquoise killifish) ]
# 107 LOC107384049          2            1,1   consistent # bhlhe40 class E basic helix-loop-helix protein 40 [ Nothobranchius furzeri (turquoise killifish) ]
# 161 LOC107396464          3        -1,1,-1 inconsistent # jun Jun proto-oncogene, AP-1 transcription factor subunit [ Nothobranchius furzeri (turquoise killifish) ]
# 214         rest          2          -1,-1   consistent

# merge (TFs that aren't significant in at least one cell type of one strain are eliminated at this stage)
my.rec.merged <- merge(my.rec.grz.v2, my.rec.zmz.v2, by = "TF_regulon", suffixes = c(".grz",".zmz"))
my.rec.merged[,-c(2,6)]
#    TF_regulon N_CellType.grz       Signs.grz Direction.grz N_CellType.zmz      Signs.zmz Direction.zmz
#         egr1              2           -1,-1    consistent              5 -1,-1,-1,-1,-1    consistent
#         irf2              3           1,1,1    consistent              5     1,-1,1,1,1  inconsistent
#         irf7              1               1    consistent              4        1,1,1,1    consistent
# LOC107375066              8 1,1,1,1,1,1,1,1    consistent              3          1,1,1    consistent # nr3c1 nuclear receptor subfamily 3, group C, member 1 (glucocorticoid receptor) [ Nothobranchius furzeri (turquoise killifish) ]
# LOC107376749              5       1,1,1,1,1    consistent              3          1,1,1    consistent # rela v-rel avian reticuloendotheliosis viral oncogene homolog A [ Nothobranchius furzeri (turquoise killifish) ]
# LOC107384049              4         1,1,1,1    consistent              2            1,1    consistent # bhlhe40 class E basic helix-loop-helix protein 40 [ Nothobranchius furzeri (turquoise killifish) ]
# LOC107385678              3           1,1,1    consistent              5      1,1,1,1,1    consistent # stat2 signal transducer and activator of transcription 2 [ Nothobranchius furzeri (turquoise killifi
# LOC107385700              1               1    consistent              5  -1,1,-1,-1,-1  inconsistent # vhl von Hippel-Lindau tumor suppressor [ Nothobranchius furzeri (turquoise killifish) 
# LOC107386825              1               1    consistent              5     1,-1,1,1,1  inconsistent # irf1b interferon regulatory factor 1b [ Nothobranchius furzeri (turquoise killifish) ]
# LOC107387991              3           1,1,1    consistent              5      1,1,1,1,1    consistent # irf3 interferon regulatory factor 3 [ Nothobranchius furzeri (turquoise killifish) ]
# LOC107396464              4        1,1,1,-1  inconsistent              3        -1,1,-1  inconsistent # jun Jun proto-oncogene, AP-1 transcription factor subunit [ Nothobranchius furzeri (turquoise killifish) ]
#        nfkb1              3           1,1,1    consistent              4        1,1,1,1    consistent
#         nfya              3        -1,-1,-1    consistent              4    -1,-1,-1,-1    consistent
#         pax6              2           -1,-1    consistent              5 -1,-1,-1,-1,-1    consistent
#         rest              4        -1,1,1,1  inconsistent              2          -1,-1    consistent
#         sox2              2           -1,-1    consistent              4    -1,-1,-1,-1    consistent
#        stat1              5       1,1,1,1,1    consistent              4        1,1,1,1    consistent

# Remove any TFs with inconsistent directionality
my.rec.merged.clean <- my.rec.merged[bitAnd(my.rec.merged$Direction.grz %in% "consistent", my.rec.merged$Direction.zmz  %in% "consistent")>0,]

# my.rec.merged.clean$TF_regulon
# [1] "egr1"         "irf7"         "LOC107375066" "LOC107376749" "LOC107384049" "LOC107385678" "LOC107387991" "nfkb1"        "nfya"         "pax6"        
# [11] "sox2"         "stat1" 

my.rec.merged.clean$TF_regulon_NAME <- c("egr1","irf7","nr3c1", "rela","bhlhe40", "stat2", "irf3", "nfkb1","nfya","pax6","sox2","stat1") 

# order
my.rec.merged.clean <- my.rec.merged.clean[order(my.rec.merged.clean$N_CellType.grz, my.rec.merged.clean$N_CellType.zmz, decreasing = T),]

write.table(my.rec.merged.clean, file = paste0(Sys.Date(),"_decoupleR_fgsea_Recurent_GRZ_ZMZ_TFregulons_summary_FDR5.txt"), sep = "\t", quote = F, row.names = F)



# make a table for plotting
my.cell.types <- sort(unique(Top_TFs_res.grz$Cell_Type))

# grab FDR10 for plotting
grz.top.tf.data <- Top_TFs_res.grz[Top_TFs_res.grz$source %in% my.rec.merged.clean$TF_regulon,]
zmz.top.tf.data <- Top_TFs_res.zmz[Top_TFs_res.zmz$source %in% my.rec.merged.clean$TF_regulon,]

my.plot.table.grz <- data.frame(matrix(0,nrow(my.rec.merged.clean),length(my.cell.types)))
my.plot.table.zmz <- data.frame(matrix(0,nrow(my.rec.merged.clean),length(my.cell.types)))
colnames(my.plot.table.grz) <- my.cell.types
colnames(my.plot.table.zmz) <- my.cell.types
rownames(my.plot.table.grz) <- my.rec.merged.clean$TF_regulon_NAME
rownames(my.plot.table.zmz) <- my.rec.merged.clean$TF_regulon_NAME

for (i in 1:nrow(my.rec.merged.clean)) {

  # regulon
  my.reg  <- my.rec.merged.clean$TF_regulon[i]
  
  # regulon data in GRZ and ZMZ
  grz.reg <- grz.top.tf.data[grz.top.tf.data$source %in% my.reg, ]
  zmz.reg <- zmz.top.tf.data[zmz.top.tf.data$source %in% my.reg, ]
  
  for (j in 1:length(my.cell.types)) {
    
    # subset for cell type
    grz.reg.ct <- grz.reg[grz.reg$Cell_Type %in% my.cell.types[j],]
    zmz.reg.ct <- zmz.reg[zmz.reg$Cell_Type %in% my.cell.types[j],]
    
    # go over GRZ data
    if (nrow(grz.reg.ct) > 0) {
      if ((grz.reg.ct$score > 0) && (grz.reg.ct$p_value < 0.05))        {
        my.plot.table.grz[i,j]  <- 1
      } else if ((grz.reg.ct$score > 0) && (grz.reg.ct$p_value < 0.1))  {
        my.plot.table.grz[i,j]   <- 0.5
      } else if ((grz.reg.ct$score < 0) && (grz.reg.ct$p_value < 0.05)) {
        my.plot.table.grz[i,j]  <- -1
      } else if ((grz.reg.ct$score < 0) && (grz.reg.ct$p_value < 0.1))  {
        my.plot.table.grz[i,j]   <- -0.5
      }
    }
    
    # go over data
    if (nrow(zmz.reg.ct) > 0) {
      if ((zmz.reg.ct$score > 0) && (zmz.reg.ct$p_value < 0.05))        {
        my.plot.table.zmz[i,j]  <- 1
      } else if ((zmz.reg.ct$score > 0) && (zmz.reg.ct$p_value < 0.1))  {
        my.plot.table.zmz[i,j]   <- 0.5
      } else if ((zmz.reg.ct$score < 0) && (zmz.reg.ct$p_value < 0.05)) {
        my.plot.table.zmz[i,j]  <- -1
      } else if ((zmz.reg.ct$score < 0) && (zmz.reg.ct$p_value < 0.1))  {
        my.plot.table.zmz[i,j]   <- -0.5
      }
    }
    
  }
  
}

grz.heat <- Heatmap(my.plot.table.grz, border = T, rect_gp = gpar(col = "grey", lwd = 0.5),
                    cluster_rows = F, cluster_columns = F,
                    column_title = "recurrent decoupleR (GRZ)")

zmz.heat <- Heatmap(my.plot.table.zmz, border = T, rect_gp = gpar(col = "grey", lwd = 0.5),
                    cluster_rows = F, cluster_columns = F,
                    column_title = "recurrent decoupleR (ZMZ)")

pdf(paste0(Sys.Date(),"_Recurrent_decoupleR_Regulons_byCellTypes_FDR5_FDR10_overlap_ZMZ_GRZ_PANELED.pdf"), width = 5.5, height = 5)
grz.heat + zmz.heat
dev.off()


#####################################################################################################################

#######################
sink(file = paste(Sys.Date(),"_R_session_Info_decoupleR_PseudoBulk_KilliBrain_Aging_Per_Cell_Type.txt", sep =""))
sessionInfo()
sink()

