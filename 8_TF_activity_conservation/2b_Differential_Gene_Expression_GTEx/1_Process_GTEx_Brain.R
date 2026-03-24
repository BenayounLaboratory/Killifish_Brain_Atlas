setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/Misc/GTEX_Human_Brain/DEseq2/')
options(stringsAsFactors = F)

library("DESeq2")        #
library("sva")           #
library("limma")         #
library("pheatmap")      #
library("bitops")      #

library(ggplot2) 
library(scales) 
theme_set(theme_bw())


# 2025-03-04
# Analyze GTEX Brain data

##################################################################
## 0a. Read in count files

# Read files in from directory that have relevant extension
files <- paste0("../Data/",list.files(path = "../Data/", pattern = "\\.gct$"))

# extract brain region name
brain.reg <- gsub("-","",gsub(".gct","",gsub("../Data/gene_reads_v10_", "", files)))

### prepare output object
my.gtex.cts <- vector(mode = "list", length(files))
names(my.gtex.cts) <- brain.reg

## read files
for (i in 1:length(files)) {
  
  # read counts
  tmp <- read.csv(files[i], header = T, sep = "\t", skip = 2)
  
  # change column name to subject
  colnames(tmp) <- unlist(lapply(strsplit(gsub("GTEX.", "GTEX_",colnames(tmp)), ".", fixed = T), '[[', 1))
  
  # aggregate by gene name
  my.gtex.cts [[i]] <- aggregate(tmp[,-c(1:2)], by = list("Gene" = tmp$Description), FUN = "sum")
  
}

save(my.gtex.cts, file = paste0(Sys.Date(), "_GTEX_v10_Brain_Counts.RData"))
##################################################################


##################################################################
## 0b. Read meta data
gtex.mdata <- read.csv("../Data/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt", header = T, sep = "\t")
# 1=Male	2=Female

gtex.mdata$SampleID <- gsub("-","_",gtex.mdata$SUBJID)
gtex.mdata$Sex      <- factor(ifelse(gtex.mdata$SEX == 1, "Male", "Female"))

# infer age to be at the middle point of the interval
gtex.mdata$Age <- NA

unique(gtex.mdata$AGE) 
# "60-69" "50-59" "40-49" "20-29" "30-39" "70-79"
gtex.mdata$Age[gtex.mdata$AGE %in% "20-29"] <- 25
gtex.mdata$Age[gtex.mdata$AGE %in% "30-39"] <- 35
gtex.mdata$Age[gtex.mdata$AGE %in% "40-49"] <- 45
gtex.mdata$Age[gtex.mdata$AGE %in% "50-59"] <- 55
gtex.mdata$Age[gtex.mdata$AGE %in% "60-69"] <- 65
gtex.mdata$Age[gtex.mdata$AGE %in% "70-79"] <- 75

# get rpw
rownames(gtex.mdata) <- gtex.mdata$SampleID
##################################################################

##############################################################################################
# 1. Use SVA to clean up batch effects on expression and DEseq2 for DE analysis

# will run SVA to clean up noise

###############################################
#######   DEG analysis   ++++   GTEX    #######
###############################################

# Create list object to receive clean SVA counts
sva.cts.gtex        <- vector(mode = "list", length = length(my.gtex.cts))
names(sva.cts.gtex) <- names(my.gtex.cts)

# Create list object to receive VST normalized counts
vst.cts.gtex        <- vector(mode = "list", length = length(my.gtex.cts))
names(vst.cts.gtex) <- names(my.gtex.cts)

# Create list object to receive DESeq2 results
deseq.res.list.gtex        <- vector(mode = "list", length = length(my.gtex.cts))
names(deseq.res.list.gtex) <- names(my.gtex.cts)

# loop over pseudobulk data
for  (i in 1:length(my.gtex.cts)) {
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_GTEX_",names(my.gtex.cts)[[i]])
  
  # format count matriX
  rownames(my.gtex.cts[[i]]) <- my.gtex.cts[[i]]$Gene
  
  # subset meta data
  gtex.mdata.TISSUE <- gtex.mdata[colnames(my.gtex.cts[[i]])[-1],]
  
  ###################################
  #######       Run SVA      #######
  
  # build design matrix
  sva.dataDesign = data.frame( row.names = gtex.mdata.TISSUE$SampleID , 
                               sex       = gtex.mdata.TISSUE$Sex      ,
                               age       = gtex.mdata.TISSUE$Age      )
  
  # Set null and alternative models (ignore batch)
  mod1    = model.matrix(~ sex + age , data = sva.dataDesign)
  n.sv.be = num.sv(my.gtex.cts[[i]][,-1], mod1, method="be") # brain_amygdala is 1
  
  # apply SVAseq algortihm
  my.svseq = svaseq(as.matrix(my.gtex.cts[[i]][,-1]), mod1, n.sv=n.sv.be, constant = 0.1)
  
  # remove RIN and SV, preserve age and sex
  my.clean <- removeBatchEffect(log2(my.gtex.cts[[i]][,-1] + 0.1), 
                                batch      = NULL, 
                                covariates = cbind(my.svseq$sv),
                                design     = mod1)
  
  # delog and round data for DEseq2 processing
  my.filtered.sva <- round(2^my.clean-0.1)
  
  # keep only robustly expressed genes
  sva.cts.gtex[[i]] <- my.filtered.sva
  
  # get matrix using age as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = sva.cts.gtex[[i]],
                                colData   = gtex.mdata.TISSUE,
                                design    = ~ Age + Sex)
  
  # run DESeq normalizations and export results
  dds.deseq <- DESeq(dds)
  
  # plot dispersion
  my.disp.out <- paste(my.outprefix,"_dispersion_plot.pdf")
  
  pdf(my.disp.out)
  plotDispEsts(dds.deseq)
  dev.off()
  
  # get DESeq2 normalized expression value
  vst.cts.gtex[[i]] <- getVarianceStabilizedData(dds.deseq)
  
  # extract gene significance by DEseq2
  res.age <- results(dds.deseq, name = "Age") # FC per week
  
  # exclude genes with NA FDR value
  res.age <- res.age[!is.na(res.age$padj),]

  # store results
  deseq.res.list.gtex[[i]]       <- data.frame(res.age)
  
  # output result tables of combined analysis to text files
  my.out.ct.mat <- paste0(my.outprefix,"_AGING_VST_log2_counts_matrix.txt")
  write.table(vst.cts.gtex[[i]], file = my.out.ct.mat , sep = "\t" , row.names = T, quote = F)
  
  my.out.stats.age <- paste0(my.outprefix,"_AGING_all_genes_statistics.txt")
  write.table(deseq.res.list.gtex[[i]], file = my.out.stats.age , sep = "\t" , row.names = T, quote = F)
  
  my.out.fdr5.age <- paste0(my.outprefix,"_AGING_FDR5_genes_statistics.txt")
  write.table(deseq.res.list.gtex[[i]][res.age$padj < 0.05,], file = my.out.fdr5.age, sep = "\t" , row.names = T, quote = F)
  

}

# save R object with all DEseq2 results
my.rdata.age <- paste0(Sys.Date(),"GTEX_Brain_AGING_DEseq2_objects.RData")
save(deseq.res.list.gtex, file = my.rdata.age)

my.vst.age <- paste0(Sys.Date(),"_GTEX_Brain_AGING_VST_data_objects.RData")
save(vst.cts.gtex, file = my.vst.age)


######### Make jitter plot of DE genes ########
## Order by pvalue:
age.results <- lapply(deseq.res.list.gtex,function(x) {x[order(x$padj),]})
n        <- sapply(age.results, nrow)
names(n) <- names(age.results)

cols <- list()
xlab <- character(length = length(age.results))
for(i in seq(along = age.results)){
  cols[[i]]      <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), n[i]) # grey60
  ind.sig.i      <- age.results[[i]]$padj < 0.05
  ind.sig.i.up   <- bitAnd(age.results[[i]]$padj < 0.05, age.results[[i]]$log2FoldChange >0)>0
  ind.sig.i.down <- bitAnd(age.results[[i]]$padj < 0.05, age.results[[i]]$log2FoldChange <0)>0
  
  cols[[i]][ind.sig.i.up]   <- "#CC3333"
  cols[[i]][ind.sig.i.down] <- "#333399"
  xlab[i] <- paste(names(age.results)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}
names(cols) <- names(age.results)

pdf(paste0(Sys.Date(),"_GTEX_Brain_Aging_perRegion_stripplot_DESeq2_with_reg_colors_FDR5.pdf"), width = 6, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 13.5),
     ylim = c(-0.2, 0.2),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change per year of life"
)
abline(h = 0)
abline(h = seq(-0.2, 0.2, by = 0.1)[-3],
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
         y = rev(age.results[[i]]$log2FoldChange),
         pch = 16,
         col = rev(cols[[i]]),
         bg = rev(cols[[i]]),
         cex = 0.75)
}
axis(1,
     at = 1:13,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-0.2, 0.2, by = 0.1))
box()
dev.off()

png(paste0(Sys.Date(),"_GTEX_Brain_Aging_perRegion_stripplot_DESeq2_with_reg_colors_FDR5.png"), width = 2400, height = 2200, res = 300, units = "px")
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 13.5),
     ylim = c(-0.2, 0.2),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change per year of life"
)
abline(h = 0)
abline(h = seq(-0.2, 0.2, by = 0.1)[-3],
       lty = "dotted",
       col = "grey")
for(i in 1:length(age.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
         y = rev(age.results[[i]]$log2FoldChange),
         pch = 16,
         col = rev(cols[[i]]),
         bg = rev(cols[[i]]),
         cex = 0.75)
}
axis(1,
     at = 1:13,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-0.2, 0.2, by = 0.1))
box()
dev.off()


pdf(paste0(Sys.Date(),"_GTEX_Brain_Aging_perRegion_stripplot_DESeq2_with_reg_colors_FDR5_NO_DOTS.pdf"), width = 6, height = 5.5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 13.5),
     ylim = c(-0.2, 0.2),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change per year of life"
)
abline(h = 0)
abline(h = seq(-0.2, 0.2, by = 0.1)[-3],
       lty = "dotted",
       col = "grey")
# for(i in 1:length(age.results)){
#   set.seed(1234)
#   points(x = jitter(rep(i, nrow(age.results[[i]])), amount = 0.2),
#          y = rev(age.results[[i]]$log2FoldChange),
#          pch = 16,
#          col = rev(cols[[i]]),
#          bg = rev(cols[[i]]),
#          cex = 0.75)
# }
axis(1,
     at = 1:13,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = gsub("brain_","",xlab),
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-0.2, 0.2, by = 0.1))
box()
dev.off()
###############################################



#######################
sink(file = paste(Sys.Date(),"_GTEX_BrainAging_DESeq2_session_Info.txt", sep =""))
sessionInfo()
sink()