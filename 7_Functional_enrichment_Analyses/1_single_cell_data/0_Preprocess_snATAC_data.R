setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/snATAC_Brain_Aging_Meta/Downstream_Analyses/GSEA')
options(stringsAsFactors = F)

#### Packages
library('muscat')         # 
library('DESeq2')         # 
library('sva')            # 
library('limma')          # 

library(ggplot2)          # 
library(scales)           # 
library("bitops")         # 
library(Vennerable)       # 
library(data.table)       #

library(ComplexHeatmap)   #
library(circlize)         #

library(phenoTest)        #
library(qusage)           #

theme_set(theme_bw())   

library(ChIPpeakAnno)
library(rtracklayer)

# 2024-06-24
# Process scRNAseq brain aging cohorts for differential gene analysis

# 2025-06-25
# rerun with correct cleaned up GO definitions

###############################################################################################
# 0. Annotate snATAC results to genes

## Load granges annotation file used for Signac run
load('/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/snATAC_Brain_Aging_Meta/Signac/2023-08-02_Killi2015_withFishTEDB_GRangesObject_Annotation.RData')
gtf3

load('../Differential_Accessiblity/2024-04-15_pseudobulk_killi_cell_types_snATAC_AGING_GRZ_DEseq2_objects.RData')
deseq.res.list.grz

# create list of granges for the coordinates since they're required by ChIPpeakAnno
deseq.granges.grz <- vector(mode = "list", length = length(deseq.res.list.grz))
names(deseq.granges.grz) <- names(deseq.res.list.grz)

for (i in 1:length(deseq.res.list.grz)) {
  
  # remove the TE sequences (since they can't be GSEAed anyway)
  deseq.res.list.grz[[i]] <- deseq.res.list.grz[[i]][-grep('NotFur',rownames(deseq.res.list.grz[[i]])),]
  
  rownames(deseq.res.list.grz[[i]]) <- gsub("NC-","NC_",rownames(deseq.res.list.grz[[i]]) )
  rownames(deseq.res.list.grz[[i]]) <- gsub("NW-","NW_",rownames(deseq.res.list.grz[[i]]) )
  
  deseq.granges.grz[[i]] <- toGRanges(data.frame("chr"   = unlist(lapply(strsplit(rownames(deseq.res.list.grz[[i]]),"-"),'[[',1)),
                                                 "start" = unlist(lapply(strsplit(rownames(deseq.res.list.grz[[i]]),"-"),'[[',2)),
                                                 "end"   = unlist(lapply(strsplit(rownames(deseq.res.list.grz[[i]]),"-"),'[[',3)) ),
                                      genome = "Killi2015")
  
  export.bed(deseq.granges.grz[[i]], con = paste0(Sys.Date(),"_", names(deseq.res.list.grz)[i],'_peaks.bed'))
}

# ChIPpeakAnno does not cooperate with custom genome
# formatted peaks will be annotated outside of R using HOMER
my.annot.files <- c("HOMER_2024-06-24_Astrocytes_Radial_Glia_peaks.xls",
                    "HOMER_2024-06-24_OPCs_peaks.xls",
                    "HOMER_2024-06-24_Oligodendrocytes_peaks.xls",
                    "HOMER_2024-06-24_NSPCs_peaks.xls",
                    "HOMER_2024-06-24_Microglia_peaks.xls",
                    "HOMER_2024-06-24_Granule_Excitatory_Neurons_peaks.xls",
                    "HOMER_2024-06-24_GABAergic_neurons_peaks.xls",
                    "HOMER_2024-06-24_Ependymal_cells_peaks.xls",
                    "HOMER_2024-06-24_PV_interneurons_peaks.xls")

#create annotation dataframe
peaks.annots.grz <- vector(mode = "list", length = length(deseq.res.list.grz))
names(peaks.annots.grz) <- names(deseq.res.list.grz)

for (i in 1:length(my.annot.files)) {
  # get HOMER annotations 
  peaks.annots.grz[[i]] <- read.csv(my.annot.files[i], sep = "\t", header = T)
  colnames(peaks.annots.grz[[i]] )[1] <- "PeakID"
  peaks.annots.grz[[i]] $PeakName <- paste(peaks.annots.grz[[i]] $Chr,peaks.annots.grz[[i]] $Start,peaks.annots.grz[[i]] $End,sep = "-")
  
  # clean gene names and genomic annotations
  peaks.annots.grz[[i]] $Gene.Name <- gsub("gene-","",peaks.annots.grz[[i]] $Nearest.PromoterID)
  peaks.annots.grz[[i]] $Genomic_Context <- unlist(lapply(strsplit(peaks.annots.grz[[i]] $Annotation, " "), '[[',1))
  # unique(peaks.annots.grz[[i]] $Genomic_Context)
}

save(peaks.annots.grz, file = paste0(Sys.Date(),"_parsed_peak_annotations_HOMER.RData"))

# use annotation to annotate DESeq2 results
deseq.res.list.grz.annot        <- vector(mode = "list", length = length(deseq.res.list.grz))
names(deseq.res.list.grz.annot) <- names(deseq.res.list.grz)

for (i in 1:length(deseq.res.list.grz)) {
  # annotate
  deseq.res.list.grz.annot[[i]] <- merge(deseq.res.list.grz[[i]], peaks.annots.grz[[i]][,c("PeakName", "Gene.Name", "Distance.to.TSS", "Genomic_Context")], by.x = 0, by.y = "PeakName")
  
}

# parse to keep only the peak closest to gene
parsed.deseq.res.list.grz.annot        <- vector(mode = "list", length = length(deseq.res.list.grz))
names(parsed.deseq.res.list.grz.annot) <- names(deseq.res.list.grz)

for (i in 1:length(parsed.deseq.res.list.grz.annot)) {

  # grab genes
  my.genes <- unique(deseq.res.list.grz.annot[[i]]$Gene.Name)
  my.keep <- c()
  
  # keep only closest peak for each gene
  for (j in 1:length(my.genes)) {
    my.rows    <- which(deseq.res.list.grz.annot[[i]]$Gene.Name %in% my.genes[j])
    my.closest <- which.min(abs(deseq.res.list.grz.annot[[i]]$Distance.to.TSS[my.rows]))
    my.keep <- c(my.keep,my.rows[my.closest])
  }
  
  # select closest
  parsed.deseq.res.list.grz.annot[[i]] <- deseq.res.list.grz.annot[[i]][my.keep,]
  
  # filter anything further than 10kb to the TSS
  parsed.deseq.res.list.grz.annot[[i]] <- parsed.deseq.res.list.grz.annot[[i]][abs(parsed.deseq.res.list.grz.annot[[i]]$Distance.to.TSS) < 10000,]

}
save(parsed.deseq.res.list.grz.annot, file = paste0(Sys.Date(),"_Annotated_Parsed_snATAC_DESeq2_aging_results.RData"))

###############################################################################################

#######################
sink(file = paste(Sys.Date(),"_MuscatDEseq2_PB_DESeq2_GSEA_scRNAseq_BrainAtlas_session_Info.txt", sep =""))
sessionInfo()
sink()