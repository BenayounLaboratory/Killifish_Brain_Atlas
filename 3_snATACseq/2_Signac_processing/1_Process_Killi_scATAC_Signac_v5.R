setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/snATAC_Brain_Aging_Meta/Signac')
options(stringsAsFactors = F)

#### Packages
library(Seurat)        # single cell general package
library(Signac)        # scATAC processing

library(GenomeInfoDb)  # for genome info
library(GenomicRanges) # for genome info
library(seqinr)        # for genome info

library(scDblFinder)   # doublet finding for scATAC
library(ggplot2)       # plotting

library(Polychrome)

# 2023-08-02
# Process scATAC cohorts 1 and 2 together for the first time
# use signac to analyze ATAC-seq data
# https://stuartlab.org/signac/0.2/articles/pbmc_vignette
# https://stuartlab.org/signac/0.2/articles/mouse_brain_vignette

# 2023-08-23
# transfer annotatiosn from scRNAseq object
# also update gft to include tx_id column to allow plotting of coverage plots

# 2024-06-05
# add a co-embedding step for visual representation

# 2025-06-27
# plot covariate UMAPs to match snRNAseq

###########################################################################################
# 0. Generate necessary annotation files for Seurat and Signac

###### a. create killi seqinfo object
killi.genome.fa <- seqinr::read.fasta(file = "/Volumes/BB_Home_HQ2/SIngle_Cell_analysis/2021-08-19_Killifish_tissue_scRNAseq/Reference/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa",
                                      seqtype = "DNA", as.string = T,
                                      set.attributes = TRUE, seqonly = FALSE, strip.desc = TRUE,
                                      whole.header = FALSE)
kili.seqnames  <- names(killi.genome.fa)
kili.seqlength <- unlist(lapply(killi.genome.fa,nchar))
killi.seqinfo  <- Seqinfo(seqnames = kili.seqnames, seqlengths = kili.seqlength, isCircular= rep(FALSE,length(killi.genome.fa)), genome="Killi2015")
rm(killi.genome.fa)
save(killi.seqinfo, file = paste0(Sys.Date(), "_Killi2015_withFishTEDB_seqinfoobject.RData") )

killi.seqinfo
# Seqinfo object with 8350 sequences from Killi2015 genome:


###### b. create granges annotation object
# We can also add gene annotations to the signac object for the killifish genome.
# This will allow downstream functions to pull the gene annotation information directly from the object.
# since killifish genomic data is not in annotation hub, need to import from gtf annotation file
gtf  <- rtracklayer::import('/Volumes/BB_Home_HQ2/SIngle_Cell_analysis/2021-08-19_Killifish_tissue_scRNAseq/Reference/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.gtf')
gtf2 <- as.data.frame(gtf)[,c("seqnames","start","end","width","strand","transcript_id","gene_name","gene_id","gene_biotype","type")]

# need to propagate biotypes
for (i in 1:nrow(gtf2)) {
  my.gid <- gtf2$gene_id[i]
  if (is.na(gtf2$gene_biotype[i])) {
    my.replacement <- unique(na.omit(gtf2$gene_biotype[gtf2$gene_id == my.gid]))
    gtf2$gene_biotype[i] <- ifelse(length(my.replacement) > 0,my.replacement,"unknown")
  }
}

gtf3 <- makeGRangesFromDataFrame(gtf2, keep.extra.columns = TRUE, seqinfo = killi.seqinfo)

# save genome annotation file
save(gtf3, file = paste0(Sys.Date(), "_Killi2015_withFishTEDB_GRangesObject_Annotation.RData") )
###########################################################################################


########################################################################################################################################################
#### 0. Assumed doublet information/to calculate %age for prediction
# Multiplet Rate (%)  # of Nuclei Loaded  # of Nuclei Recovered
# 0.4%                 ~775                ~500
# 0.8%                 ~1,550              ~1,000
# 1.6%                 ~3,075              ~2,000
# 2.3%                 ~4,625              ~3,000
# 3.1%                 ~6,150              ~4,000
# 3.9%                 ~7,700              ~5,000
# 4.6%                 ~9,250              ~6,000
# 5.4%                 ~10,750             ~7,000
# 6.2%                 ~12,300             ~8,000
# 6.9%                 ~13,850             ~9,000
# 7.7%                 ~15,400             ~10,000

pred.10x.dblt <- data.frame( "nuc_number" = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000),
                             "dblt_rate"   = c(0.4, 0.8, 1.6, 2.3, 3.1, 3.9, 4.6, 5.4, 6.2, 6.9, 7.7))

pred_dblt_lm <- lm(dblt_rate ~ nuc_number, data = pred.10x.dblt)

pdf(paste0(Sys.Date(), "_10xATAC_nuclei_number_vs_doublet_rate.pdf"))
plot(dblt_rate ~ nuc_number, data = pred.10x.dblt)
abline(pred_dblt_lm, col = "red", lty = "dashed")
dev.off()
########################################################################################################################################################


###########################################################################################
# 1. Load cellranger-atac data for each sample separately 

# Merge individual objects into one (cell ranger aggregate doesn't work for killi genome)
# https://satijalab.org/signac/articles/merging.html


# When pre-processing chromatin data, Signac uses information from two related input files, 
# both of which can be created using CellRanger output: 
#   - Peak/Cell matrix.
#   - Fragment file.

#### Read and combine peaksets
# read in CellRanger peak sets
# https://stackoverflow.com/questions/28433328/skip-comment-line-in-csv-file-using-r
peaks.YF1 <- read.table(file = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1/CellRanger/Killi_Brain_YF1_FishTEDB_NR/outs/peaks.bed", col.names = c("chr", "start", "end"), comment.char = '#')
peaks.OF1 <- read.table(file = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1/CellRanger/Killi_Brain_OF1_FishTEDB_NR/outs/peaks.bed", col.names = c("chr", "start", "end"), comment.char = '#')
peaks.YM1 <- read.table(file = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1/CellRanger/Killi_Brain_YM1_FishTEDB_NR/outs/peaks.bed", col.names = c("chr", "start", "end"), comment.char = '#')
peaks.OM1 <- read.table(file = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1/CellRanger/Killi_Brain_OM1_FishTEDB_NR/outs/peaks.bed", col.names = c("chr", "start", "end"), comment.char = '#')
peaks.YF2 <- read.table(file = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2/CellRanger/Killi_Brain_YF2_FishTEDB_NR/outs/peaks.bed", col.names = c("chr", "start", "end"), comment.char = '#')
peaks.OF2 <- read.table(file = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2/CellRanger/Killi_Brain_OF2_FishTEDB_NR/outs/peaks.bed", col.names = c("chr", "start", "end"), comment.char = '#')
peaks.YM2 <- read.table(file = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2/CellRanger/Killi_Brain_YM2_FishTEDB_NR/outs/peaks.bed", col.names = c("chr", "start", "end"), comment.char = '#')
peaks.OM2 <- read.table(file = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2/CellRanger/Killi_Brain_OM2_FishTEDB_NR/outs/peaks.bed", col.names = c("chr", "start", "end"), comment.char = '#')

# convert peaks to genomic ranges
gr.YF1 <- makeGRangesFromDataFrame(peaks.YF1, seqinfo = killi.seqinfo)
gr.OF1 <- makeGRangesFromDataFrame(peaks.OF1, seqinfo = killi.seqinfo)
gr.YM1 <- makeGRangesFromDataFrame(peaks.YM1, seqinfo = killi.seqinfo)
gr.OM1 <- makeGRangesFromDataFrame(peaks.OM1, seqinfo = killi.seqinfo)
gr.YF2 <- makeGRangesFromDataFrame(peaks.YF2, seqinfo = killi.seqinfo)
gr.OF2 <- makeGRangesFromDataFrame(peaks.OF2, seqinfo = killi.seqinfo)
gr.YM2 <- makeGRangesFromDataFrame(peaks.YM2, seqinfo = killi.seqinfo)
gr.OM2 <- makeGRangesFromDataFrame(peaks.OM2, seqinfo = killi.seqinfo)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.YF1, gr.OF1, gr.YM1, gr.OM1,
                               gr.YF2, gr.OF2, gr.YM2, gr.OM2))

# Filter out bad peaks based on length
peakwidths       <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
# GRanges object with 126692 ranges and 0 metadata columns:

####  Create Fragment objects
# load metadata
md.YF1 <- read.csv(file = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1/CellRanger/Killi_Brain_YF1_FishTEDB_NR/outs/singlecell.csv", header = TRUE, row.names = 1)[-1, ] # remove the first row
md.OF1 <- read.csv(file = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1/CellRanger/Killi_Brain_OF1_FishTEDB_NR/outs/singlecell.csv", header = TRUE, row.names = 1)[-1, ] # remove the first row
md.YM1 <- read.csv(file = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1/CellRanger/Killi_Brain_YM1_FishTEDB_NR/outs/singlecell.csv", header = TRUE, row.names = 1)[-1, ] # remove the first row
md.OM1 <- read.csv(file = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1/CellRanger/Killi_Brain_OM1_FishTEDB_NR/outs/singlecell.csv", header = TRUE, row.names = 1)[-1, ] # remove the first row
md.YF2 <- read.csv(file = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2/CellRanger/Killi_Brain_YF2_FishTEDB_NR/outs/singlecell.csv", header = TRUE, row.names = 1)[-1, ] # remove the first row
md.OF2 <- read.csv(file = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2/CellRanger/Killi_Brain_OF2_FishTEDB_NR/outs/singlecell.csv", header = TRUE, row.names = 1)[-1, ] # remove the first row
md.YM2 <- read.csv(file = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2/CellRanger/Killi_Brain_YM2_FishTEDB_NR/outs/singlecell.csv", header = TRUE, row.names = 1)[-1, ] # remove the first row
md.OM2 <- read.csv(file = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2/CellRanger/Killi_Brain_OM2_FishTEDB_NR/outs/singlecell.csv", header = TRUE, row.names = 1)[-1, ] # remove the first row


# perform an initial filtering of low count cells
md.YF1 <- md.YF1[md.YF1$passed_filters > 250, ]
md.OF1 <- md.OF1[md.OF1$passed_filters > 250, ]
md.YM1 <- md.YM1[md.YM1$passed_filters > 250, ]
md.OM1 <- md.OM1[md.OM1$passed_filters > 250, ]
md.YF2 <- md.YF2[md.YF2$passed_filters > 250, ]
md.OF2 <- md.OF2[md.OF2$passed_filters > 250, ]
md.YM2 <- md.YM2[md.YM2$passed_filters > 250, ]
md.OM2 <- md.OM2[md.OM2$passed_filters > 250, ]

# create fragment objects
frags.YF1 <- CreateFragmentObject(path = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1/CellRanger/Killi_Brain_YF1_FishTEDB_NR/outs/fragments.tsv.gz", cells = rownames(md.YF1))
frags.OF1 <- CreateFragmentObject(path = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1/CellRanger/Killi_Brain_OF1_FishTEDB_NR/outs/fragments.tsv.gz", cells = rownames(md.OF1))
frags.YM1 <- CreateFragmentObject(path = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1/CellRanger/Killi_Brain_YM1_FishTEDB_NR/outs/fragments.tsv.gz", cells = rownames(md.YM1))
frags.OM1 <- CreateFragmentObject(path = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1/CellRanger/Killi_Brain_OM1_FishTEDB_NR/outs/fragments.tsv.gz", cells = rownames(md.OM1))
frags.YF2 <- CreateFragmentObject(path = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2/CellRanger/Killi_Brain_YF2_FishTEDB_NR/outs/fragments.tsv.gz", cells = rownames(md.YF2))
frags.OF2 <- CreateFragmentObject(path = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2/CellRanger/Killi_Brain_OF2_FishTEDB_NR/outs/fragments.tsv.gz", cells = rownames(md.OF2))
frags.YM2 <- CreateFragmentObject(path = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2/CellRanger/Killi_Brain_YM2_FishTEDB_NR/outs/fragments.tsv.gz", cells = rownames(md.YM2))
frags.OM2 <- CreateFragmentObject(path = "/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2/CellRanger/Killi_Brain_OM2_FishTEDB_NR/outs/fragments.tsv.gz", cells = rownames(md.OM2))

#### Quantify peaks in each dataset
# We can now create a matrix of peaks x cell for each sample using the FeatureMatrix function. 
# This function is parallelized using the future package. 
# See the parallelization vignette for more information about using future.

YF1.counts <- FeatureMatrix( fragments = frags.YF1, features = combined.peaks, cells = rownames(md.YF1))
OF1.counts <- FeatureMatrix( fragments = frags.OF1, features = combined.peaks, cells = rownames(md.OF1))
YM1.counts <- FeatureMatrix( fragments = frags.YM1, features = combined.peaks, cells = rownames(md.YM1))
OM1.counts <- FeatureMatrix( fragments = frags.OM1, features = combined.peaks, cells = rownames(md.OM1))
YF2.counts <- FeatureMatrix( fragments = frags.YF2, features = combined.peaks, cells = rownames(md.YF2))
OF2.counts <- FeatureMatrix( fragments = frags.OF2, features = combined.peaks, cells = rownames(md.OF2))
YM2.counts <- FeatureMatrix( fragments = frags.YM2, features = combined.peaks, cells = rownames(md.YM2))
OM2.counts <- FeatureMatrix( fragments = frags.OM2, features = combined.peaks, cells = rownames(md.OM2))

#### Create the objects
# We will now use the quantified matrices to create a Seurat object for each dataset, 
# storing the Fragment object for each dataset in the assay.
YF1.assay  <- CreateChromatinAssay(YF1.counts, fragments = frags.YF1)
OF1.assay  <- CreateChromatinAssay(OF1.counts, fragments = frags.OF1)
YM1.assay  <- CreateChromatinAssay(YM1.counts, fragments = frags.YM1)
OM1.assay  <- CreateChromatinAssay(OM1.counts, fragments = frags.OM1)
YF2.assay  <- CreateChromatinAssay(YF2.counts, fragments = frags.YF2)
OF2.assay  <- CreateChromatinAssay(OF2.counts, fragments = frags.OF2)
YM2.assay  <- CreateChromatinAssay(YM2.counts, fragments = frags.YM2)
OM2.assay  <- CreateChromatinAssay(OM2.counts, fragments = frags.OM2)


YF1.seurat <- CreateSeuratObject(YF1.assay, assay = "ATAC", meta.data=md.YF1)
OF1.seurat <- CreateSeuratObject(OF1.assay, assay = "ATAC", meta.data=md.OF1)
YM1.seurat <- CreateSeuratObject(YM1.assay, assay = "ATAC", meta.data=md.YM1)
OM1.seurat <- CreateSeuratObject(OM1.assay, assay = "ATAC", meta.data=md.OM1)
YF2.seurat <- CreateSeuratObject(YF2.assay, assay = "ATAC", meta.data=md.YF2)
OF2.seurat <- CreateSeuratObject(OF2.assay, assay = "ATAC", meta.data=md.OF2)
YM2.seurat <- CreateSeuratObject(YM2.assay, assay = "ATAC", meta.data=md.YM2)
OM2.seurat <- CreateSeuratObject(OM2.assay, assay = "ATAC", meta.data=md.OM2)

# Merge objects
# Now that the objects each contain an assay with the same set of features, 
# we can use the standard merge function to merge the objects. 
# This will also merge all the fragment objects so that we retain the fragment information
# for each cell in the final merged object.

# add information to identify dataset of origin
YF1.seurat$sample  <- 'YF1'
OF1.seurat$sample  <- 'OF1'
YM1.seurat$sample  <- 'YM1'
OM1.seurat$sample  <- 'OM1'
YF2.seurat$sample  <- 'YF2'
OF2.seurat$sample  <- 'OF2'
YM2.seurat$sample  <- 'YM2'
OM2.seurat$sample  <- 'OM2'

YF1.seurat$sex  <- 'F'
OF1.seurat$sex  <- 'F'
YM1.seurat$sex  <- 'M'
OM1.seurat$sex  <- 'M'
YF2.seurat$sex  <- 'F'
OF2.seurat$sex  <- 'F'
YM2.seurat$sex  <- 'M'
OM2.seurat$sex  <- 'M'

YF1.seurat$Batch  <- 'Cohort1'
OF1.seurat$Batch  <- 'Cohort1'
YM1.seurat$Batch  <- 'Cohort1'
OM1.seurat$Batch  <- 'Cohort1'
YF2.seurat$Batch  <- 'Cohort2'
OF2.seurat$Batch  <- 'Cohort2'
YM2.seurat$Batch  <- 'Cohort2'
OM2.seurat$Batch  <- 'Cohort2'


# get average ages per group/Cohort from Ari's file (cohorts 1 and 3) [each sample derived from pool of 3 animals]
YF1.seurat$age_weeks  <- 6.29
OF1.seurat$age_weeks  <- 15.14
YM1.seurat$age_weeks  <- 6.29
OM1.seurat$age_weeks  <- 15.14
YF2.seurat$age_weeks  <- 6.14
OF2.seurat$age_weeks  <- 16.23
YM2.seurat$age_weeks  <- 6.14
OM2.seurat$age_weeks  <- 16.57

# merge all datasets, adding a cell ID to make sure cell names are unique
brain.atac.combined <- merge(x = YF1.seurat,
                             y = list(OF1.seurat, YM1.seurat, OM1.seurat,
                                      YF2.seurat,OF2.seurat, YM2.seurat, OM2.seurat),
                             add.cell.ids = c("YF1", "OF1", "YM1", "OM1",
                                              "YF2", "OF2", "YM2", "OM2"))

# add the gene information to the object
Annotation(brain.atac.combined) <- gtf3

brain.atac.combined[["ATAC"]]
# ChromatinAssay data with 126692 features for 51212 cells
# Variable features: 0 
# Genome: 
#   Annotation present: TRUE 
# Motifs present: FALSE 
# Fragment files: 8 

save(brain.atac.combined, file = paste0(Sys.Date(), "_Raw_combined_killi_scATAC_Signac_object.RData"))
###########################################################################################



###########################################################################################
# 2. Cell QC and filtering

# compute nucleosome signal score per cell
brain.atac.combined <- NucleosomeSignal(object = brain.atac.combined)
# Found 13644 cell barcodes
# Done Processing 75 million lines
# Found 5174 cell barcodes
# Done Processing 83 million lines
# Found 4729 cell barcodes
# Done Processing 37 million lines
# Found 2089 cell barcodes
# Done Processing 46 million lines
# Found 6838 cell barcodes
# Done Processing 75 million lines
# Found 5837 cell barcodes
# Done Processing 66 million lines
# Found 6491 cell barcodes
# Done Processing 88 million lines
# Found 6410 cell barcodes
# Done Processing 70 million lines

# compute TSS enrichment score per cell
brain.atac.combined <- TSSEnrichment(object = brain.atac.combined, fast = FALSE)

# add fraction of reads in peaks
brain.atac.combined$pct_reads_in_peaks <- brain.atac.combined$peak_region_fragments / brain.atac.combined$passed_filters * 100
brain.atac.combined$high.tss           <- ifelse(brain.atac.combined$TSS.enrichment > 2, 'High', 'Low')

pdf(paste0(Sys.Date(),"_TSS_enrichment_plot_Killi_brain_ATAC_2Cohorts.pdf"))
TSSPlot(brain.atac.combined, group.by = 'high.tss') + NoLegend()
TSSPlot(brain.atac.combined, group.by = 'sample') + NoLegend()
dev.off()

pdf(paste0(Sys.Date(),"_QC_Violing_plots_Killi_brain_ATAC_2Cohorts.pdf"), height = 5, width = 13)
VlnPlot(object = brain.atac.combined,
        features = c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'nucleosome_signal', 'mitochondrial'),
        pt.size = 0.01, ncol = 5, group.by = 'sample')
dev.off()


pdf(paste0(Sys.Date(),"_QC_Violing_plots_Killi_brain_ATAC_2Cohorts_noPoints.pdf"), height = 5, width = 13)
VlnPlot(object = brain.atac.combined,
        features = c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'nucleosome_signal', 'mitochondrial'),
        pt.size = 0, ncol = 5, group.by = 'sample')
dev.off()

brain.atac.combined$nucleosome_group <- ifelse(brain.atac.combined$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
table(brain.atac.combined$nucleosome_group)
# NS < 2 NS > 2
#  50604    608 

# pdf(paste0(Sys.Date(),"_QC_FragmentHistogram_Killi_brain_ATAC_2Cohorts.pdf"), height = 5, width = 13)
# FragmentHistogram(object = brain.atac.combined, group.by = 'nucleosome_group', region = combined.peaks) ### too much memory
# dev.off()

# Finally we remove cells that are outliers for these QC metrics.
brain.atac.combined <- subset(x      = brain.atac.combined,
                              subset = pct_reads_in_peaks >  50     &
                                peak_region_fragments     >  3000   &
                                peak_region_fragments     <  40000  &
                                nucleosome_signal         <  2      &
                                TSS.enrichment            >  2      &
                                mitochondrial             <  1000   )
brain.atac.combined
# An object of class Seurat 
# 126692 features across 29829 samples within 1 assay 
# Active assay: ATAC (126692 features, 0 variable features)

table(brain.atac.combined@meta.data$sample)
#  OF1  OF2  OM1  OM2  YF1  YF2  YM1  YM2 
# 4120 4013 1526 4707 3536 4537 3390 4000 

table(brain.atac.combined@meta.data$Batch)
# Cohort1 Cohort2 
# 12572   17257 

pdf(paste0(Sys.Date(),"_QC_Violing_plots_Killi_brain_ATAC_2Cohorts_POST_Filter.pdf"), height = 5, width = 13)
VlnPlot(object = brain.atac.combined,
        features = c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'nucleosome_signal', 'mitochondrial'),
        pt.size = 0.01, ncol = 5, group.by = 'sample', raster = T)
dev.off()

pdf(paste0(Sys.Date(),"_QC_Violing_plots_Killi_brain_ATAC_2Cohorts_POST_Filter_NoPoints.pdf"), height = 5, width = 13)
VlnPlot(object = brain.atac.combined,
        features = c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'nucleosome_signal', 'mitochondrial'),
        pt.size = 0, ncol = 5, group.by = 'sample')
dev.off()

save(brain.atac.combined, file = paste0(Sys.Date(), "_QCFiltered_combined_killi_scATAC_Signac_object.RData"))
###########################################################################################


###########################################################################################
# 3. Doublet removal using scDoublet
# https://www.bioconductor.org/packages/devel/bioc/vignettes/scDblFinder/inst/doc/scATAC.html
### error on gird - need to updater delayedarray https://github.com/plger/scDblFinder/issues/81

# load QC filter object
load('2023-08-03_QCFiltered_combined_killi_scATAC_Signac_object.RData')

# add biological group
brain.atac.combined$Group <- NA
brain.atac.combined$Group[grep("YF", brain.atac.combined$sample)] <- "YF"
brain.atac.combined$Group[grep("OF", brain.atac.combined$sample)] <- "OF"
brain.atac.combined$Group[grep("YM", brain.atac.combined$sample)] <- "YM"
brain.atac.combined$Group[grep("OM", brain.atac.combined$sample)] <- "OM"

# ## Assume doublet rate based on 10x information (add a 10% fudge factor due to nuclei being more sticky)
# pred.dblt.rate <- 1.10 * predict(pred_dblt_lm, data.frame("nuc_number" = unlist(lapply(atac.list, ncol))))/100
# #     YF1        OF1        YM1        OM1        YF2        OF2        YM2        OM2
# # 0.03017409 0.03509768 0.02894320 0.01322819 0.03861332 0.03419558 0.03408598 0.04004656

# convert to SCE object
brain.atac.sce <- as.SingleCellExperiment(brain.atac.combined)

# Run scDblFinder
brain.atac.sce <- scDblFinder(brain.atac.sce,
                              artificialDoublets = 10000            ,
                              samples            = "sample"         , # split per library
                              aggregateFeatures  = TRUE             ,
                              nfeatures          = 1000             ,
                              processing         = "normFeatures"   ,
                              verbose            = T                )

# Check order
sum(rownames(brain.atac.combined@meta.data) == colnames(brain.atac.sce)) # 29829

# transfer doublet/singlet calls
brain.atac.combined@meta.data$scDblFinder <- brain.atac.sce$scDblFinder.class
# brain.atac.combined@meta.data

save(brain.atac.combined, file = paste0(Sys.Date(),"_Killi_brain_ATAC_2Cohorts_FiltNorm_doublets_annotated.RData") )
###########################################################################################


###########################################################################################
# 4. Dimensionality reduction and clustering

# Normalization and linear dimensional reduction
# The combined steps of TF-IDF followed by SVD are known as latent semantic indexing (LSI),
# and were first introduced for the analysis of scATAC-seq data by Cusanovich et al. 2015.
brain.atac.combined <- RunTFIDF(brain.atac.combined)
brain.atac.combined <- FindTopFeatures(brain.atac.combined, min.cutoff = 'q70') # top 30% features for UMAP
brain.atac.combined <- RunSVD(brain.atac.combined)

# The first LSI component often captures sequencing depth (technical variation) rather than biological variation.
# If this is the case, the component should be removed from downstream analysis. We can assess the correlation between each LSI component and sequencing depth using the DepthCor() function:
pdf(paste0(Sys.Date(),"_Depth_SVD_corplot_Killi_brain_ATAC_2Cohorts.pdf"), height = 5, width = 13)
DepthCor(brain.atac.combined)
dev.off()

# run non linear dim reduction / use number of dimensions from RNAseq (19)
brain.atac.combined <- RunUMAP(brain.atac.combined, dims = 2:19, reduction = 'lsi')

# plot by Group and by doublet call
pdf(paste0(Sys.Date(),"_UMAP_by_sample_Killi_brain_ATAC_2Cohorts.pdf"), height = 5, width = 6)
DimPlot(brain.atac.combined, group.by = 'Group', shuffle = T, cols = c("deeppink4","deepskyblue4","deeppink","deepskyblue"), raster = T, raster.dpi = c(600,600))
dev.off()

pdf(paste0(Sys.Date(),"_UMAP_DOublet_labels_Killi_brain_ATAC_2Cohorts.pdf"), height = 5, width = 6)
DimPlot(brain.atac.combined, group.by = 'scDblFinder', shuffle = T, raster = T, raster.dpi = c(600,600))
dev.off()

pdf(paste0(Sys.Date(),"_UMAP_bySex_Killi_brain_ATAC_2Cohorts.pdf"), height = 5, width = 6)
DimPlot(brain.atac.combined, group.by = 'sex', shuffle = T, cols = c("deeppink","deepskyblue"), raster = T, raster.dpi = c(600,600))
dev.off()

# How many doublets
table(brain.atac.combined@meta.data$scDblFinder)
# singlet doublet 
# 27171    2658 

# remove doublets
brain.atac.singlets <- subset(x = brain.atac.combined, subset = scDblFinder %in% 'singlet'  )
brain.atac.singlets
# An object of class Seurat 
# 126692 features across 27171 samples within 1 assay 
# Active assay: ATAC (126692 features, 38018 variable features)
#  2 dimensional reductions calculated: lsi, umap


# rerun non linear dim reduction / use number of dimensions from RNAseq (19)
brain.atac.singlets <- RunUMAP(brain.atac.singlets, dims = 2:19, reduction = 'lsi')

pdf(paste0(Sys.Date(),"_UMAP_by_sample_Killi_brain_ATAC_2Cohorts_SINGLETS.pdf"), height = 5, width = 6)
DimPlot(brain.atac.singlets, group.by = 'Group', shuffle = T, cols = c("deeppink4","deepskyblue4","deeppink","deepskyblue"), raster = T, raster.dpi = c(600,600))
dev.off()


# Clustering of singlets
# Now that the cells are embedded in a low-dimensional space, we can use methods commonly applied
# for the analysis of scRNA-seq data to perform graph-based clustering and non-linear dimension reduction for visualization.
brain.atac.singlets <- FindNeighbors(object = brain.atac.singlets, reduction = 'lsi', dims = 2:19)
brain.atac.singlets <- FindClusters(object = brain.atac.singlets, verbose = FALSE, algorithm = 3)

pdf(paste0(Sys.Date(),"_UMAP_withSeurat_Clusters_Killi_brain_ATAC_2Cohorts_SINGLETS.pdf"), height = 5, width = 6)
DimPlot(object = brain.atac.singlets, label = TRUE, raster = T, raster.dpi = c(600,600))
dev.off()

save(brain.atac.singlets, file = paste0(Sys.Date(),"_Killi_brain_ATAC_2Cohorts_FiltNorm_Singlets.RData") )
###########################################################################################


###########################################################################################
# 5. Calculate a gene activity matrix

# split Seurat by tissue of origin
# (Full object runs but fails with a matrix merging error - https://github.com/stuart-lab/signac/issues/1399)
atac.samp.list <- SplitObject(brain.atac.singlets, split.by = "sample")

gene.act.list  <- vector(mode = "list", length = 8)

for (i in 1:8) {
  # Create a gene activity matrix
  gene.act.list[[i]] <- GeneActivity(atac.samp.list[[i]], assay = "ATAC" ) 
  
}

# some genes are not found in all
lapply(gene.act.list,nrow)
# [[1]] 22143
# [[2]] 22143
# [[3]] 22140
# [[4]] 22141
# [[5]] 22145
# [[6]] 22145
# [[7]] 22147
# [[8]] 22142

# get the genes that could be identified in all
genes.use <- names(table(unlist(lapply(gene.act.list,row.names))))[table(unlist(lapply(gene.act.list,row.names))) == 8] # 22135

# select rows for genes available in all libraries
select_genes <- function(matrix, genes = genes.use) {
  matrix[genes.use,]
}

gene.act.list.filt <- lapply(gene.act.list, select_genes)
lapply(gene.act.list.filt,nrow) # all are 22135 now

# merge to creat one output
gene.activities <- Reduce(cbind, gene.act.list.filt[-1], gene.act.list.filt[[1]])

# add the gene activity matrix to the Seurat object as a new assay and normalize it
brain.atac.singlets[['RNA']] <- CreateAssayObject(counts = gene.activities)
brain.atac.singlets          <- NormalizeData(object = brain.atac.singlets,
                                              assay = 'RNA',
                                              normalization.method = 'LogNormalize',
                                              scale.factor = median(brain.atac.singlets$nCount_RNA) )

DefaultAssay(brain.atac.singlets) <- 'RNA'

#### Plot markers
pdf(paste(Sys.Date(),"Killi_brain_ATAC_2Cohorts_SINGLETS_Brain_Markers_UMAP_GeneActivity.pdf", sep = "_"), height = 5, width = 6.5)
FeaturePlot(brain.atac.singlets, features = c("olig1","olig2","mpz")                               , max.cutoff = 'q95', pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Oligodendrocyte/OPC
FeaturePlot(brain.atac.singlets, features = c("marco","csf1r", "ptprc")                            , max.cutoff = 'q95', pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # mph/microglia
FeaturePlot(brain.atac.singlets, features = c("s100b", "slc1a2")                                   , max.cutoff = 'q95', pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # astrocyte/radial glia
FeaturePlot(brain.atac.singlets, features = c("rbfox3", "map2", "eno2", "ncam1")                   , max.cutoff = 'q95', pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # mature neuron
FeaturePlot(brain.atac.singlets, features = c("dcx")                                               , max.cutoff = 'q95', pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # immature neuron
FeaturePlot(brain.atac.singlets, features = c("fat2", "neurod1", "eomes", "pax6")                  , max.cutoff = 'q95', pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # Granule excitatory neuron
FeaturePlot(brain.atac.singlets, features = c("pvalb")                                             , max.cutoff = 'q95', pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # pvalb interneurons
FeaturePlot(brain.atac.singlets, features = c("gad1", "gad2",  "LOC107384443", "LOC107391088")     , max.cutoff = 'q95', pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # GABAergic neurons
FeaturePlot(brain.atac.singlets, features = c("LOC107386767")                                      , max.cutoff = 'q95', pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # dopaminergic neurons
FeaturePlot(brain.atac.singlets, features = c("slc17a6","LOC107381463")                            , max.cutoff = 'q95', pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # glutamatergic neurons
FeaturePlot(brain.atac.singlets, features = c("LOC107386535",  "clu")                              , max.cutoff = 'q95', pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # ependymal cells
FeaturePlot(brain.atac.singlets, features = c("sox5",  "sox2")                                     , max.cutoff = 'q95', pt.size	= 2 , raster = T, raster.dpi = c(1024, 1024)) # NSPCs
dev.off()

save(brain.atac.singlets, file = paste0(Sys.Date(),"_Killi_brain_ATAC_2Cohorts_FiltNorm_Singlets_withGeneActivity.RData") )
###########################################################################################


###########################################################################################
# 6. Integrating with scRNA-seq data

# To help interpret the scATAC-seq data, we can classify cells based on an scRNA-seq experiment
# from the same biological system (killifish brain snRNAseq).
# Identify shared correlation patterns in the gene activity matrix and scRNA-seq dataset

# Load the ATAC object for killibrain
load('2023-08-03_Killi_brain_ATAC_2Cohorts_FiltNorm_Singlets_withGeneActivity.RData')
brain.atac.singlets
# An object of class Seurat 
# 148827 features across 27171 samples within 2 assays 
# Active assay: RNA (22135 features, 0 variable features)
# 1 other assay present: ATAC
# 2 dimensional reductions calculated: lsi, umap

# Load the annotated RNA Seurat object for killibrain
load('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Atlas_Analysis/2023-08-23_Seurat_object_with_Manual_annotation_FINAL.RData')
killi.brain.clean
# An object of class Seurat 
# 21160 features across 209939 samples within 1 assay 
# Active assay: RNA (21160 features, 5000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony

kill.brain.rna <- killi.brain.clean; rm(killi.brain.clean)
DefaultAssay(kill.brain.rna)      <- "RNA"
DefaultAssay(brain.atac.singlets) <- "RNA"

# runs out of memory (likely due to file sizes)
# for annotation transfer, will just use a downsampled RNA object
# First subset GRZ (since ATAC is just GRZ)
killi.brain.rna.grz   <- subset(killi.brain.rna, subset = Strain %in% "GRZ")
killi.brain.rna.grz
# An object of class Seurat 
# 21160 features across 104665 samples within 1 assay 
# Active assay: RNA (21160 features, 5000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony

# Then DS per condition, 2500 cells per conditions
killi.brain.rna.grz       <- SetIdent(object = killi.brain.rna.grz, value = 'Group')
killi.brain.rna.grz.small <- subset(killi.brain.rna.grz, downsample = 2500)
killi.brain.rna.grz.small
# An object of class Seurat 
# 21160 features across 15000 samples within 1 assay 
# Active assay: RNA (21160 features, 5000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony

# Check representation
set.seed(123123) # stabilize
P16 = createPalette(16+3,  c("#ff0000", "#00ff00", "#0000ff"))
swatch(P16)
killi.brain.rna.grz.small@meta.data$Cell_Identity <- factor(x = killi.brain.rna.grz.small@meta.data$Cell_Identity, levels = sort(levels(killi.brain.rna.grz.small@meta.data$Cell_Identity)))

pdf(paste(Sys.Date(),"Killifish_Brain_RNA_UMAP_color_by_deNovo_annotation_RASTER_GRZ_15K_subset.pdf", sep = "_"), height = 4, width = 7)
DimPlot(killi.brain.rna.grz.small, reduction = "umap", group.by = "Cell_Identity", raster = T,
        raster.dpi = c(350,350), cols = as.vector(P16[-c(1:3)]))
dev.off()

# find transfer anchors
transfer.anchors <- FindTransferAnchors(reference            = killi.brain.rna.grz.small,
                                        query                = brain.atac.singlets,
                                        reduction            = 'cca',
                                        normalization.method = 'LogNormalize')
# Running CCA
# Merging objects
# Finding neighborhoods
# Finding anchors
# Found 28500 anchors
# Filtering anchors
# Retained 9044 anchors

predicted.labels <- TransferData(anchorset        = transfer.anchors,
                                 refdata          = kill.brain.rna.grz.small$Cell_Identity,
                                 weight.reduction = brain.atac.singlets[['lsi']],
                                 dims             = 2:19)

brain.atac.singlets <- AddMetaData(object = brain.atac.singlets, metadata = predicted.labels)

brain.atac.singlets@meta.data$predicted.id <- factor(x = brain.atac.singlets@meta.data$predicted.id, levels = sort(levels(brain.atac.singlets@meta.data$predicted.id)))
brain.atac.singlets@meta.data$Cell_Identity <- as.factor(brain.atac.singlets@meta.data$predicted.id)

# plot annotated UMAPs
plot1 <- DimPlot( object = kill.brain.rna.grz.small,
                  group.by = 'Cell_Identity',
                  label = TRUE,
                  repel = TRUE,
                  raster = T,
                  raster.dpi = c(500,500), 
                  cols = as.vector(P16[-c(1:3)]) ) + NoLegend() + ggtitle('scRNA-seq (GRZ, Downsampled)')

plot2 <- DimPlot( object = brain.atac.singlets,
                  group.by = 'Cell_Identity',
                  label = TRUE,
                  repel = TRUE,
                  raster = T,
                  raster.dpi = c(500,500), 
                  cols = as.vector(P16[-c(1:3)])) + NoLegend() + ggtitle('scATAC-seq')

pdf(paste0(Sys.Date(),"_Side_by_side_UMAPs_Killi_brain_aging_snRNAseq_Cp1234_snATACseq_C12.pdf"), height = 6, width = 13)
plot1 + plot2
dev.off()


plot1 <- DimPlot( object = kill.brain.rna.grz.small,
                  group.by = 'Cell_Identity',
                  label = F,
                  raster = T,
                  raster.dpi = c(500,500), 
                  cols = as.vector(P16[-c(1:3)]) ) + NoLegend() + ggtitle('scRNA-seq (GRZ, Downsampled)')

plot2 <- DimPlot( object = brain.atac.singlets,
                  group.by = 'Cell_Identity',
                  label = F,
                  raster = T,
                  raster.dpi = c(500,500), 
                  cols = as.vector(P16[-c(1:3)])) + NoLegend() + ggtitle('scATAC-seq')

pdf(paste0(Sys.Date(),"_Side_by_side_UMAPs_Killi_brain_aging_snRNAseq_Cp1234_snATACseq_C12_NO_LABEL.pdf"), height = 6, width = 13)
plot1 + plot2
dev.off()

table(brain.atac.singlets$predicted.id)
#       Astrocytes_Radial_Glia              Ependymal_cells                 Erythrocytes            GABAergic_neurons   Granule_Excitatory_Neurons 
#                         1084                          229                           85                          815                        10667 
#                    Microglia               Neurons_misc_1               Neurons_misc_2               Neurons_misc_3               Neurons_misc_4 
#                          566                          124                          758                          510                         9358 
#                        NSPCs             Oligodendrocytes                         OPCs               Purkinje_cells              PV_interneurons 
#                          496                         1491                          283                           57                          449 
# Vascular_smooth_muscle_cells 
#                          199 

save(brain.atac.singlets, file = paste0(Sys.Date(),"_Killi_brain_ATAC_2Cohorts_FiltNorm_Singlets_ANNOTATEDfromRNA.RData") )
###########################################################################################


###########################################################################################
# 7. Co-embedding

# https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette
# similar to the AD brain paper Morabito et al, 2021

# load ATAC
load('2023-08-23_Killi_brain_ATAC_2Cohorts_FiltNorm_Singlets_ANNOTATEDfromRNA.RData')
brain.atac.singlets
# An object of class Seurat 
# 148827 features across 27171 samples within 2 assays 
# Active assay: RNA (22135 features, 0 variable features)
# 1 other assay present: ATAC
# 2 dimensional reductions calculated: lsi, umap

# Load the annotated RNA Seurat object for killibrain
load('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Atlas_Analysis/2023-08-23_Seurat_object_with_Manual_annotation_FINAL.RData')
killi.brain.clean
# An object of class Seurat 
# 21160 features across 209939 samples within 1 assay 
# Active assay: RNA (21160 features, 5000 variable features)
# 3 dimensional reductions calculated: pca, umap, harmony

killi.brain.rna <- killi.brain.clean; rm(killi.brain.clean)
DefaultAssay(killi.brain.rna)      <- "RNA"
DefaultAssay(brain.atac.singlets) <- "RNA"

# find transfer anchors
transfer.anchors <- FindTransferAnchors(reference            = killi.brain.rna.grz,
                                        query                = brain.atac.singlets,
                                        reduction            = 'cca',
                                        normalization.method = 'LogNormalize')

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the full transcriptome if we wanted to
genes.use         <- VariableFeatures(killi.brain.rna.grz)
killi.brain.ref   <- GetAssayData(killi.brain.rna.grz, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation                   <- TransferData(anchorset = transfer.anchors, refdata = killi.brain.ref, weight.reduction = brain.atac.singlets[["lsi"]], dims = 2:19)
brain.atac.singlets[["RNA"]] <- imputation

# move annotation
brain.atac.singlets@meta.data$Cell_Identity <- brain.atac.singlets@meta.data$predicted.id
brain.atac.singlets@meta.data$assay <- "ATAC"
killi.brain.rna.grz@meta.data$assay <- "RNA"

# Merge for joint representation
coembed <- merge(x = killi.brain.rna.grz, y = brain.atac.singlets)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)


# Check representation
set.seed(123123) # stabilize
P16 = createPalette(16+3,  c("#ff0000", "#00ff00", "#0000ff"))
swatch(P16)
# coembed@meta.data$Cell_Identity <- factor(x = coembed@meta.data$Cell_Identity, levels = sort(levels(coembed@meta.data$Cell_Identity)))

pdf(paste(Sys.Date(),"Killifish_Brain_RNA_ATAC_coEMBED_UMAP_color_by_ASSAY_RASTER_GRZOnly.pdf", sep = "_"), height = 4, width = 5)
DimPlot(coembed, group.by = "assay", raster = T,raster.dpi = c(350,350), cols = c("darkorange","cyan4"), shuffle = T)
dev.off()

pdf(paste(Sys.Date(),"Killifish_Brain_RNA_ATAC_coEMBED_UMAP_color_by_deNovo_annotation_RASTER_GRZOnly.pdf", sep = "_"), height = 4, width = 7)
DimPlot(coembed , group.by = "Cell_Identity", raster = T,raster.dpi = c(350,350), cols = as.vector(P16[-c(1:3)]))
dev.off()

pdf(paste(Sys.Date(),"Killifish_Brain_RNA_ATAC_coEMBED_UMAP_color_by_ASSAY_RASTER_GRZOnly_SPLIT.pdf", sep = "_"), height = 4, width = 7.5)
DimPlot(coembed, group.by = "assay", split.by = "assay", raster = T,raster.dpi = c(350,350), cols = c("darkorange","cyan4"), shuffle = T)
dev.off()

save(coembed, file = paste0(Sys.Date(),"_coembedded_object.RData"))
###########################################################################################

###########################################################################################
# 8. Cell type markers plotting chromatin plotting

# change back to working with peaks instead of gene activities
DefaultAssay(brain.atac.singlets) <- 'ATAC'
brain.atac.singlets <- SetIdent(brain.atac.singlets, value = 'Group')


###############################################################
# to fix track plotting
# https://github.com/stuart-lab/signac/issues/1159
load('2023-08-02_Killi2015_withFishTEDB_GRangesObject_Annotation.RData')
gtf3$tx_id <- gtf3$transcript_id

# update annotation
Annotation(brain.atac.singlets) <- gtf3
###############################################################


pdf(paste0(Sys.Date(),"_s100b_Astrocytes_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "s100b",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()

pdf(paste0(Sys.Date(),"_mpz_oligo_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "mpz",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()

pdf(paste0(Sys.Date(),"_ncam1_neuron_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "ncam1",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()

pdf(paste0(Sys.Date(),"_olig2_OPC_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "olig2",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()


pdf(paste0(Sys.Date(),"_csf1r_microglia_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "csf1r",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()

pdf(paste0(Sys.Date(),"_LOC107387973_itgam_microglia_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "LOC107387973",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()

pdf(paste0(Sys.Date(),"_LOC107379395_apoeb_microglia_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "LOC107379395",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()

pdf(paste0(Sys.Date(),"_LOC107386535_ependymin-2-like_ependymal_cells_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "LOC107386535",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()

pdf(paste0(Sys.Date(),"_LOC107378372_hbba_Erythrocytes_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "LOC107378372",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()

pdf(paste0(Sys.Date(),"_gli3_NSPCs_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "gli3",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()

pdf(paste0(Sys.Date(),"_sema5a_OPCs_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "sema5a",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()

pdf(paste0(Sys.Date(),"_tpm1_VascularSmoothMuscle_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "tpm1",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()

pdf(paste0(Sys.Date(),"_LOC107384443_gad1l_GABAergic_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "LOC107384443",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()

pdf(paste0(Sys.Date(),"_fat2_GranuleExcitatoryNeuron_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "fat2",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()

pdf(paste0(Sys.Date(),"_itpr1_PurkinjeCells_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "itpr1",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()


pdf(paste0(Sys.Date(),"_pvalb_PVInterneuron_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "pvalb",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()

pdf(paste0(Sys.Date(),"_map2_Neuron_by_cell_type_scATAC_GRZ.pdf"), height = 3, width = 5)
CoveragePlot(object            = brain.atac.singlets,
             region            = "map2",
             extend.upstream   = 10000,
             extend.downstream = 10000,
             group.by = "Cell_Identity")
dev.off()
###########################################################################################


###########################################################################################
# 8. Plot co-variates in UMAP form to match snRNAseq

load("2023-08-23_Killi_brain_ATAC_2Cohorts_FiltNorm_Singlets_ANNOTATEDfromRNA.RData")
brain.atac.singlets
# An object of class Seurat 
# 148827 features across 27171 samples within 2 assays 
# Active assay: RNA (22135 features, 0 variable features)
# 2 layers present: counts, data
# 1 other assay present: ATAC
# 2 dimensional reductions calculated: lsi, umap

# Check representation
set.seed(123123) # stabilize
P16 = createPalette(16+3,  c("#ff0000", "#00ff00", "#0000ff"))
swatch(P16)

# add age group
brain.atac.singlets$Age_Group <- ifelse(brain.atac.singlets$Group %in% c("YF","YM"), "Young", "Old")

# Plot annotation output
pdf(paste(Sys.Date(),"Killifish_Brain_UMAP_snATAC_GRZ_color_by_transfer_annotation_RASTER.pdf", sep = "_"), height = 4, width = 7)
DimPlot(brain.atac.singlets, reduction = "umap", group.by = "predicted.id", raster = T,
        raster.dpi = c(350,350), cols = as.vector(P16[-c(1:3)]))
dev.off()



pdf(paste(Sys.Date(),"Killifish_Brain_UMAP_snATAC_GRZ_Singlets_UMAP_color_by_Sex.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(brain.atac.singlets, reduction = "umap", group.by = "sex", pt.size	= 2, shuffle = T,
        cols = c(alpha("deeppink", alpha = 0.5 ), alpha("deepskyblue", alpha = 0.5 ) ), raster = T, raster.dpi = c(1024, 1024) )
dev.off()


pdf(paste(Sys.Date(),"Killifish_Brain_UMAP_snATAC_GRZ_Singlets_UMAP_color_by_Age.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(brain.atac.singlets, reduction = "umap", group.by = "Age_Group", pt.size	= 2, shuffle = T,
        cols = c(alpha("darkorange", alpha = 0.5 ), alpha("darkolivegreen2", alpha = 0.5 ) ),
        raster = T, raster.dpi = c(1024, 1024))
dev.off()

pdf(paste(Sys.Date(),"Killifish_Brain_UMAP_snATAC_GRZ_Singlets_UMAP_color_by_Batch.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(brain.atac.singlets, reduction = "umap", group.by = "Batch", pt.size	= 2, shuffle = T, raster = T, raster.dpi = c(1024, 1024))
dev.off()




# reorder cell types for ease of plotting (non neuronal/neuronal)
brain.atac.singlets@meta.data$predicted.id <- factor(x = brain.atac.singlets@meta.data$predicted.id, 
                                                     levels = c("Astrocytes_Radial_Glia"      , 
                                                                "Ependymal_cells"             ,
                                                                "Erythrocytes"                ,
                                                                "Microglia"                   ,
                                                                "NSPCs"                       ,
                                                                "Oligodendrocytes"            ,
                                                                "OPCs"                        ,
                                                                "Vascular_smooth_muscle_cells",
                                                                "GABAergic_neurons"           ,
                                                                "Granule_Excitatory_Neurons"  ,       
                                                                "Neurons_misc_1"              ,
                                                                "Neurons_misc_2"              ,
                                                                "Neurons_misc_3"              ,
                                                                "Neurons_misc_4"              ,       
                                                                "Purkinje_cells"              ,
                                                                "PV_interneurons"              ))

###### Plot marker genes - mix of previously known and from the marker calculation above
pdf(paste(Sys.Date(),"Killifish_Brain_snATAC_GRZ_Dotplot_KnownAndPredicted_Cell_type_markers_by_TransferAnnotation.pdf", sep = "_"), height = 8, width = 16)
DotPlot(brain.atac.singlets,
        features = c("s100b", "slc1a2","kcnj10",                              # astrocyte/radial glia
                     "LOC107386535", "LOC107383970","clu",                    # ependymal cells               ### LOC107386535 (ependymin-2-like), LOC107383970	(serotransferrin-like)
                     "LOC107378372", "LOC107378374","LOC107378381",           # Erythrocytes (hemoglobin)
                     "LOC107387973", "LOC107379395", "csf1r",                 # mph/microglia                 ### LOC107387973 (itgam/cd11b), LOC107379395	(apoeb)
                     "gli3","efna2", "LOC107385497",                          # NSPCs                         ### gli3, efna2, LOC107385497 (msi1)
                     "mpz","LOC107394899", "LOC107386530",                    # Oligodendrocyte               ### LOC107394899 (mbpa), LOC107386530 (plp)
                     "olig1","olig2", "sema5a","sox5",                        # OPCs
                     "tpm1", "LOC107375895", "ahnak",                         # Vascular smooth muscle        ### LOC107375895 (kcnq5)
                     
                     "LOC107384443","gad2", "LOC107391088",                   # GABAergic neurons             ### LOC107384443 (gad1l), LOC107391088 (scl6a1; GABA transporter)
                     "fat2", "LOC107392205", "LOC107376653",                  # Granule excitatory neuron     ### LOC107392205 (qka), LOC107376653 (Slc17a7/VGLUT1)
                     "LOC107392873","LOC107397057","LOC107387921",            # Neurons_misc_1                ### LOC107392873 (fcho1),LOC107397057 (gjc1),LOC107387921 (ngfr)
                     "LOC107385736","LOC107385201", "cdh24",                  # Neurons_misc_2                ### LOC107385736 (grm2), LOC107385201 (camkv)
                     "elavl4", "tmem163","LOC107380370",                      # Neurons_misc_3                ### LOC107380370 (elnl)
                     "nrp1","meis2","LOC107390816",                           # Neurons_misc_4                ### LOC107375999 (ebf3) / LOC107390816 cpne5a
                     "itpr1","mpped1","LOC107379976",                         # Purkinje cells                ### LOC107379976 (aldoca)
                     "pvalb", "pde2a",                                        # PV interneurons
                     "rbfox3", "map2", "ncam1"                                # mature neuron markers
        ),
        group.by = "predicted.id") + RotatedAxis()
dev.off()


###########################################################################################

#######################
sink(file = paste(Sys.Date(),"_Seurat_Signac_scATAC_session_Info.txt", sep =""))
sessionInfo()
sink()
