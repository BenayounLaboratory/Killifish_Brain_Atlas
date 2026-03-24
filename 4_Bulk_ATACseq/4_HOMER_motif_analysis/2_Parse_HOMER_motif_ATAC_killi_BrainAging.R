setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/ATAC/Bulk_ATAC_Analysis_Folder/HOMER_Motif_Analysis/')
options(stringsAsFactors = F)

# 2024-04-12
# parse/plot ATAC motif enrichment data

# 2024-05-01
# replot only top 5 each direction
# correct code error where some motifs were enriched in up and down peaks

############################################
# 1. Read in HOMER motif enrichment in sex-specifc changes with age ATAC peaks
ATAC.motif.GRZ.up  <- read.csv("Killi_Brain_Aging_GRZ_UP/knownResults.txt" , sep = "\t", header = T)
ATAC.motif.GRZ.dwn <- read.csv("Killi_Brain_Aging_GRZ_DOWN/knownResults.txt"  , sep = "\t", header = T)
ATAC.motif.ZMZ.up  <- read.csv("Killi_Brain_Aging_ZMZ_UP/knownResults.txt" , sep = "\t", header = T)
ATAC.motif.ZMZ.dwn <- read.csv("Killi_Brain_Aging_ZMZ_DOWN/knownResults.txt"  , sep = "\t", header = T)

# get clean colnames
colnames(ATAC.motif.GRZ.up ) <- c("Motif.Name"            ,
                                "Consensus"             ,
                                "pvalue"                ,
                                "Log_pvalue"            ,
                                "q.value"               ,
                                "Num_FG_with_motif"     ,
                                "Percent_FG_with_motif" ,
                                "Num_BG_with_motif"     ,
                                "Percent_BG_with_motif" )
colnames(ATAC.motif.GRZ.dwn) <- colnames(ATAC.motif.GRZ.up)
colnames(ATAC.motif.ZMZ.up ) <- colnames(ATAC.motif.GRZ.up)
colnames(ATAC.motif.ZMZ.dwn) <- colnames(ATAC.motif.GRZ.up)

# extract percentages as numbers
ATAC.motif.GRZ.up$Percent_FG_with_motif  <- as.numeric(sub("%", "", ATAC.motif.GRZ.up$Percent_FG_with_motif )) # get percent number
ATAC.motif.GRZ.up$Percent_BG_with_motif  <- as.numeric(sub("%", "", ATAC.motif.GRZ.up$Percent_BG_with_motif )) # get percent number
ATAC.motif.GRZ.dwn$Percent_FG_with_motif <- as.numeric(sub("%", "", ATAC.motif.GRZ.dwn$Percent_FG_with_motif)) # get percent number
ATAC.motif.GRZ.dwn$Percent_BG_with_motif <- as.numeric(sub("%", "", ATAC.motif.GRZ.dwn$Percent_BG_with_motif)) # get percent number
ATAC.motif.ZMZ.up$Percent_FG_with_motif  <- as.numeric(sub("%", "", ATAC.motif.ZMZ.up$Percent_FG_with_motif )) # get percent number
ATAC.motif.ZMZ.up$Percent_BG_with_motif  <- as.numeric(sub("%", "", ATAC.motif.ZMZ.up$Percent_BG_with_motif )) # get percent number
ATAC.motif.ZMZ.dwn$Percent_FG_with_motif <- as.numeric(sub("%", "", ATAC.motif.ZMZ.dwn$Percent_FG_with_motif)) # get percent number
ATAC.motif.ZMZ.dwn$Percent_BG_with_motif <- as.numeric(sub("%", "", ATAC.motif.ZMZ.dwn$Percent_BG_with_motif)) # get percent number

ATAC.motif.GRZ.up$EnrichFold  <- ATAC.motif.GRZ.up$Percent_FG_with_motif/ATAC.motif.GRZ.up$Percent_BG_with_motif
ATAC.motif.GRZ.dwn$EnrichFold <- ATAC.motif.GRZ.dwn$Percent_FG_with_motif/ATAC.motif.GRZ.dwn$Percent_BG_with_motif
ATAC.motif.ZMZ.up$EnrichFold  <- ATAC.motif.ZMZ.up$Percent_FG_with_motif/ATAC.motif.ZMZ.up$Percent_BG_with_motif
ATAC.motif.ZMZ.dwn$EnrichFold <- ATAC.motif.ZMZ.dwn$Percent_FG_with_motif/ATAC.motif.ZMZ.dwn$Percent_BG_with_motif

# Filter to retain only significantly enriched
sig.ATAC.motif.GRZ.up  <- ATAC.motif.GRZ.up  [ATAC.motif.GRZ.up  $q.value < 0.05,]  # 130
sig.ATAC.motif.GRZ.dwn <- ATAC.motif.GRZ.dwn [ATAC.motif.GRZ.dwn $q.value < 0.05,]  # 128
sig.ATAC.motif.ZMZ.up  <- ATAC.motif.ZMZ.up  [ATAC.motif.ZMZ.up  $q.value < 0.05,]  # 115
sig.ATAC.motif.ZMZ.dwn <- ATAC.motif.ZMZ.dwn [ATAC.motif.ZMZ.dwn $q.value < 0.05,]  # 111


# Put up and down together
sig.ATAC.motif.GRZ                <- rbind(sig.ATAC.motif.GRZ.up,sig.ATAC.motif.GRZ.dwn)
sig.ATAC.motif.GRZ$Direction      <- c(rep("UP",nrow(sig.ATAC.motif.GRZ.up)),rep("DOWN",nrow(sig.ATAC.motif.GRZ.dwn)))
sig.ATAC.motif.GRZ$minusLog10pval <- -sig.ATAC.motif.GRZ$Log_pvalue
sig.ATAC.motif.GRZ$EnrichFold[-c(1:nrow(sig.ATAC.motif.GRZ.up))] <-  -sig.ATAC.motif.GRZ$EnrichFold[-c(1:nrow(sig.ATAC.motif.GRZ.up))]   # make negative values for enrichment in regions of decreased accessibility

sig.ATAC.motif.ZMZ                <- rbind(sig.ATAC.motif.ZMZ.up,sig.ATAC.motif.ZMZ.dwn)
sig.ATAC.motif.ZMZ$Direction      <- c(rep("UP",nrow(sig.ATAC.motif.ZMZ.up)),rep("DOWN",nrow(sig.ATAC.motif.ZMZ.dwn)))
sig.ATAC.motif.ZMZ$minusLog10pval <- -sig.ATAC.motif.ZMZ$Log_pvalue
sig.ATAC.motif.ZMZ$EnrichFold[-c(1:nrow(sig.ATAC.motif.ZMZ.up))] <-  -sig.ATAC.motif.ZMZ$EnrichFold[-c(1:nrow(sig.ATAC.motif.ZMZ.up))]   # make negative values for enrichment in regions of decreased accessibility

write.table(sig.ATAC.motif.GRZ, file = paste0(Sys.Date(),"_HOMER_Motifs_Enriched_in_GRZ_ATAC_DA_peaks_FDR5.txt"), quote = F, sep = "\t", row.names = F)
write.table(sig.ATAC.motif.ZMZ, file = paste0(Sys.Date(),"_HOMER_Motifs_Enriched_in_ZMZ_ATAC_DA_peaks_FDR5.txt"), quote = F, sep = "\t", row.names = F)

# Grab TFs enriched in both strains in the same direction
up.both  <- intersect(sig.ATAC.motif.GRZ.up$Motif.Name ,sig.ATAC.motif.ZMZ.up$Motif.Name ) # 79
dwn.both <- intersect(sig.ATAC.motif.GRZ.dwn$Motif.Name,sig.ATAC.motif.ZMZ.dwn$Motif.Name) # 96

up.merged <- merge(sig.ATAC.motif.GRZ[sig.ATAC.motif.GRZ$Motif.Name %in% up.both,-c(4,12)], sig.ATAC.motif.ZMZ[sig.ATAC.motif.ZMZ$Motif.Name %in% up.both,-c(2,4,12)], by = "Motif.Name", suffixes = c('.GRZ',".ZMZ"))
dwn.merged <- merge(sig.ATAC.motif.GRZ[sig.ATAC.motif.GRZ$Motif.Name %in% dwn.both,-c(4,12)], sig.ATAC.motif.ZMZ[sig.ATAC.motif.ZMZ$Motif.Name %in% dwn.both,-c(2,4,12)], by = "Motif.Name", suffixes = c('.GRZ',".ZMZ"))

# some TFs are in up and down lists (KLF14)
# we can filter by looking for repeated names
sum(duplicated(up.merged$Motif.Name) ) # 0
sum(duplicated(dwn.merged$Motif.Name) ) #14
dwn.motif.conflict <- dwn.merged$Motif.Name[duplicated(dwn.merged$Motif.Name)]

dwn.merged <- dwn.merged[!(dwn.merged$Motif.Name %in% dwn.motif.conflict),]

# get rank
up.merged$GRZ_Rank <- rank(-up.merged$EnrichFold.GRZ)  # so it's decreasing order
up.merged$ZMZ_Rank <- rank(-up.merged$EnrichFold.ZMZ)  # so it's decreasing order
dwn.merged$GRZ_Rank <- rank(dwn.merged$EnrichFold.GRZ) # so it's decreasing order (already neg values)
dwn.merged$ZMZ_Rank <- rank(dwn.merged$EnrichFold.ZMZ) # so it's decreasing order (already neg values)

up.merged$Rank_Product  <- up.merged$GRZ_Rank  * up.merged$ZMZ_Rank
dwn.merged$Rank_Product <- dwn.merged$GRZ_Rank * dwn.merged$ZMZ_Rank

sig.ATAC.motif.both <- rbind(up.merged, dwn.merged)
write.table(sig.ATAC.motif.both, file = paste0(Sys.Date(),"_HOMER_Motifs_Enriched_in_GRZ_and_ZMZ_ATAC_DA_peaks_FDR5.txt"), quote = F, sep = "\t", row.names = F)

############################################################################################
# 2. Make bubble chart summary
library(ggplot2) 
library(scales) 
theme_set(theme_bw())

#### Grab top 10 motifs each way
up.top10  <- sort(up.merged$Rank_Product , index.return = T)$ix[1:10]
dwn.top10 <- sort(dwn.merged$Rank_Product, index.return = T)$ix[rev(1:10)]

up.merged.top10  <- up.merged[up.top10,]
dwn.merged.top10 <- dwn.merged[dwn.top10,]

# get filtered/merged datafame for ggplot
my.top.data<- rbind(up.merged.top10, dwn.merged.top10)

## clean up motif names
my.top.data$Motif.Name <- gsub("/Homer","",my.top.data$Motif.Name)
my.top.data$Motif.Name <- gsub("/Promoter","",my.top.data$Motif.Name)
my.top.data$Motif.Name <- unlist(lapply(strsplit(my.top.data$Motif.Name,"-ChIP-Seq"),'[[',1))
my.top.data$Motif.Name <- unlist(lapply(strsplit(my.top.data$Motif.Name,"-ChIP-seq"),'[[',1))
my.top.data$Motif.Name <- gsub("/","_",my.top.data$Motif.Name)

my.top.data.grz <- my.top.data[,c("Motif.Name","q.value.GRZ","EnrichFold.GRZ")]
my.top.data.zmz <- my.top.data[,c("Motif.Name","q.value.ZMZ","EnrichFold.ZMZ")]
colnames(my.top.data.grz) <- c("Motif.Name","qvalue","Enrich_Score")
colnames(my.top.data.zmz) <- c("Motif.Name","qvalue","Enrich_Score")

my.top.data.grz$Strain <- "GRZ"
my.top.data.zmz$Strain <- "ZMZ1001"

my.top.data.v2  <- rbind(my.top.data.grz,my.top.data.zmz)
my.top.data.v2$minusLog10qval <- -log10(my.top.data.v2$qvalue+1e-5)

#### 
my.max <- max(my.top.data.v2$Enrich_Score)
my.min <- min(my.top.data.v2$Enrich_Score)
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
my.top.data.v2$Motif.Name  <- factor(my.top.data.v2$Motif.Name, levels = rev(my.top.data$Motif.Name))

pdf(paste0(Sys.Date(),"HOMER_Motifs_Enriched_ATAC_DA_KilliBrainAging_bothStrains_FDR5_Top10.pdf"), height = 5, width=5.5)
my.plot <- ggplot(my.top.data.v2,aes(x=Strain,y=Motif.Name,colour=Enrich_Score,size=minusLog10qval))+ theme_bw() + geom_point(shape = 16)
my.plot <- my.plot + ggtitle("HOMER motifs") + labs(x = "Killifish Brain Aging", y = "Aging DA ATAC motifs (FDR < 5%)")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled)
my.plot <- my.plot + scale_x_discrete(guide = guide_axis(angle = 45)) + scale_size(range = c(3, 10)) + scale_size_area(limits=c(2,5))
print(my.plot)
dev.off()


#### Grab top 5 motifs each way
up.top5  <- sort(up.merged$Rank_Product , index.return = T)$ix[1:5]
dwn.top5 <- sort(dwn.merged$Rank_Product, index.return = T)$ix[rev(1:5)]

up.merged.top5  <- up.merged[up.top5,]
dwn.merged.top5 <- dwn.merged[dwn.top5,]

# get filtered/merged datafame for ggplot
my.top.data<- rbind(up.merged.top5, dwn.merged.top5)

## clean up motif names
my.top.data$Motif.Name <- gsub("/Homer","",my.top.data$Motif.Name)
my.top.data$Motif.Name <- gsub("/Promoter","",my.top.data$Motif.Name)
my.top.data$Motif.Name <- unlist(lapply(strsplit(my.top.data$Motif.Name,"-ChIP-Seq"),'[[',1))
my.top.data$Motif.Name <- unlist(lapply(strsplit(my.top.data$Motif.Name,"-ChIP-seq"),'[[',1))
my.top.data$Motif.Name <- gsub("/","_",my.top.data$Motif.Name)

my.top.data.grz <- my.top.data[,c("Motif.Name","q.value.GRZ","EnrichFold.GRZ")]
my.top.data.zmz <- my.top.data[,c("Motif.Name","q.value.ZMZ","EnrichFold.ZMZ")]
colnames(my.top.data.grz) <- c("Motif.Name","qvalue","Enrich_Score")
colnames(my.top.data.zmz) <- c("Motif.Name","qvalue","Enrich_Score")

my.top.data.grz$Strain <- "GRZ"
my.top.data.zmz$Strain <- "ZMZ1001"

my.top.data.v2  <- rbind(my.top.data.grz,my.top.data.zmz)
my.top.data.v2$minusLog10qval <- -log10(my.top.data.v2$qvalue+1e-5)

#### 
my.max <- max(my.top.data.v2$Enrich_Score)
my.min <- min(my.top.data.v2$Enrich_Score)
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
my.top.data.v2$Motif.Name  <- factor(my.top.data.v2$Motif.Name, levels = rev(my.top.data$Motif.Name))

pdf(paste0(Sys.Date(),"HOMER_Motifs_Enriched_ATAC_DA_KilliBrainAging_bothStrains_FDR5_Top5.pdf"), height = 5, width=5.5)
my.plot <- ggplot(my.top.data.v2,aes(x=Strain,y=Motif.Name,colour=Enrich_Score,size=minusLog10qval))+ theme_bw() + geom_point(shape = 16)
my.plot <- my.plot + ggtitle("HOMER motifs") + labs(x = "Killifish Brain Aging", y = "Aging DA ATAC motifs (FDR < 5%)")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled)
my.plot <- my.plot + scale_x_discrete(guide = guide_axis(angle = 45)) + scale_size(range = c(3, 10)) + scale_size_area(limits=c(2,5))
print(my.plot)
dev.off()

