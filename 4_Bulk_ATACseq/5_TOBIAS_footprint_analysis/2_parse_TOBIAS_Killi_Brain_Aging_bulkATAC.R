setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/ATAC/Bulk_ATAC_Analysis_Folder/TOBIAS')
options(stringsAsFactors = F)

# Loading necessary libraries
library(pheatmap)
library(bitops)
library(grDevices)
library(ggplot2) 
library(scales) 
theme_set(theme_bw())

# 2024-05-01
# parse TOBIAS BINDetect results

###################################################################################################################
# 0. process TOBIAS results

# read bindtect analysis results
grz.tob.Fo <- read.csv("GRZ_YF_OF_BINDetect_JASPAR/bindetect_results.txt", header = T, sep = "\t")
grz.tob.Mo <- read.csv("GRZ_YM_OM_BINDetect_JASPAR/bindetect_results.txt", header = T, sep = "\t")
zmz.tob.Fo <- read.csv("ZMZ_YF_OF_BINDetect_JASPAR/bindetect_results.txt", header = T, sep = "\t")
zmz.tob.Mo <- read.csv("ZMZ_YM_OM_BINDetect_JASPAR/bindetect_results.txt", header = T, sep = "\t")
zmz.tob.Fg <- read.csv("ZMZ_YF_GF_BINDetect_JASPAR/bindetect_results.txt", header = T, sep = "\t")
zmz.tob.Mg <- read.csv("ZMZ_YM_OG_BINDetect_JASPAR/bindetect_results.txt", header = T, sep = "\t")

# run FDR adjustment
grz.tob.Fo$FDR <- p.adjust(grz.tob.Fo$YF_OF_pvalue, method = "fdr")
grz.tob.Mo$FDR <- p.adjust(grz.tob.Mo$YM_OM_pvalue, method = "fdr")
zmz.tob.Fo$FDR <- p.adjust(zmz.tob.Fo$YF_OF_pvalue, method = "fdr")
zmz.tob.Mo$FDR <- p.adjust(zmz.tob.Mo$YM_OM_pvalue, method = "fdr")
zmz.tob.Fg$FDR <- p.adjust(zmz.tob.Fg$YF_GF_pvalue, method = "fdr")
zmz.tob.Mg$FDR <- p.adjust(zmz.tob.Mg$YM_GO_pvalue, method = "fdr")

# use neg of the score so positive means up with aging
grz.tob.Fo$Aging_Change <- - grz.tob.Fo$YF_OF_change
grz.tob.Mo$Aging_Change <- - grz.tob.Mo$YM_OM_change
zmz.tob.Fo$Aging_Change <- - zmz.tob.Fo$YF_OF_change
zmz.tob.Mo$Aging_Change <- - zmz.tob.Mo$YM_OM_change
zmz.tob.Fg$Aging_Change <- - zmz.tob.Fg$YF_GF_change
zmz.tob.Mg$Aging_Change <- - zmz.tob.Mg$YM_GO_change

# merge for comparison
grz.o.merged <- merge(grz.tob.Fo,grz.tob.Mo[,-c(2:4)], by = "output_prefix", suffixes = c(".F",".M"))
zmz.o.merged <- merge(zmz.tob.Fo,zmz.tob.Mo[,-c(2:4)], by = "output_prefix", suffixes = c(".F",".M"))
zmz.g.merged <- merge(zmz.tob.Fg,zmz.tob.Mg[,-c(2:4)], by = "output_prefix", suffixes = c(".F",".M"))
###################################################################################################################

###################################################################################################################
# 1. plot scatters

############# GRZ ############# 
# perform test
grz.test <- cor.test(grz.o.merged$Aging_Change.F, grz.o.merged$Aging_Change.M,method = "spearman")

nr3c1.grz <- which(grz.o.merged$name == "NR3C1")
grz.o.merged[nr3c1.grz,]
#       output_prefix  name motif_id cluster total_tfbs.F YF_mean_score YF_bound OF_mean_score OF_bound YF_OF_change YF_OF_pvalue        FDR.F Aging_Change.F
# 1219 NR3C1_MA0113.3 NR3C1 MA0113.3 C_NR3C2          960       1.54883      233       1.70988      268     -0.27104  1.61814e-93 2.068681e-92        0.27104
#      total_tfbs.M YM_mean_score YM_bound OM_mean_score OM_bound YM_OM_change YM_OM_pvalue         FDR.M Aging_Change.M
# 1219          960       1.57386      247       1.79983      279     -0.36849 1.31587e-106 2.236821e-105        0.36849


pdf(paste0(Sys.Date(),"_Bindetect_TOBIAS_GRZ_bulkATAC_BRAIN_FM_Scatter_YvsO.pdf"))
plot(grz.o.merged$Aging_Change.F, grz.o.merged$Aging_Change.M, 
     pch = 16, 
     col = adjustcolor( "black", alpha.f = 0.3),
     xlim = c(-0.5,0.5),
     ylim = c(-0.5,0.5),
     xlab = "Differential binding score in Young vs. Old (Females)",
     ylab = "Differential binding score in Young vs. Old (Males)",
     main = "TOBIAS footprint analysis (GRZ)")
abline(h = 0, col = "grey", lty = "dashed")
abline(v = 0, col = "grey", lty = "dashed")
abline(0 , 1, col = "red", lty = "dashed")
text(-0.5,0.5, paste0("Rho = ",signif(grz.test$estimate,3)), pos = 4)
text(-0.5,0.45, paste0("p = ",signif(grz.test$p.value,3)), pos = 4)
points(grz.o.merged$Aging_Change.F[nr3c1.grz], grz.o.merged$Aging_Change.M[nr3c1.grz],
       pch = 16, col = "#CC3333")
text(grz.o.merged$Aging_Change.F[nr3c1.grz], grz.o.merged$Aging_Change.M[nr3c1.grz],"Nr3c1",col = "#CC3333", pos = 4)
dev.off()



############# ZMZ ############# 
# perform test
zmz.test <- cor.test(zmz.o.merged$Aging_Change.F, zmz.o.merged$Aging_Change.M,method = "spearman")

nr3c1.zmz <- which(zmz.o.merged$name == "NR3C1")
zmz.o.merged[nr3c1.zmz,]
#       output_prefix  name motif_id cluster total_tfbs.F YF_mean_score YF_bound OF_mean_score OF_bound YF_OF_change YF_OF_pvalue         FDR.F Aging_Change.F
# 1219 NR3C1_MA0113.3 NR3C1 MA0113.3 C_NR3C2          960       1.52188      236        1.7799      275     -0.35197 1.19812e-104 8.679714e-103        0.35197
#      total_tfbs.M YM_mean_score YM_bound OM_mean_score OM_bound YM_OM_change YM_OM_pvalue        FDR.M Aging_Change.M
# 1219          960       1.78526      248       2.07339      286     -0.29921  1.00862e-97 1.078066e-96        0.29921


pdf(paste0(Sys.Date(),"_Bindetect_TOBIAS_ZMZ_bulkATAC_BRAIN_FM_Scatter_YvsO.pdf"))
plot(zmz.o.merged$Aging_Change.F, zmz.o.merged$Aging_Change.M, 
     pch = 16, 
     col = adjustcolor( "black", alpha.f = 0.3),
     xlim = c(-0.5,0.5),
     ylim = c(-0.5,0.5),
     xlab = "Differential binding score in Young vs. Old (Females)",
     ylab = "Differential binding score in Young vs. Old (Males)",
     main = "TOBIAS footprint analysis (ZMZ)")
abline(h = 0, col = "grey", lty = "dashed")
abline(v = 0, col = "grey", lty = "dashed")
abline(0 , 1, col = "red", lty = "dashed")
text(-0.5,0.5, paste0("Rho = ",signif(zmz.test$estimate,3)), pos = 4)
text(-0.5,0.45, paste0("p = ",signif(zmz.test$p.value,3)), pos = 4)
points(zmz.o.merged$Aging_Change.F[nr3c1.zmz], zmz.o.merged$Aging_Change.M[nr3c1.zmz],
       pch = 16, col = "#CC3333")
text(zmz.o.merged$Aging_Change.F[nr3c1.zmz], zmz.o.merged$Aging_Change.M[nr3c1.zmz],"Nr3c1",col = "#CC3333", pos = 4)
dev.off()


####
zmz.test.2 <- cor.test(zmz.g.merged$Aging_Change.F, zmz.g.merged$Aging_Change.M,method = "spearman")

nr3c1.zmz.2 <- which(zmz.g.merged$name == "NR3C1")
zmz.g.merged[nr3c1.zmz.2,]
#       output_prefix  name motif_id cluster total_tfbs.F YF_mean_score YF_bound GF_mean_score GF_bound YF_GF_change YF_GF_pvalue        FDR.F Aging_Change.F
# 1219 NR3C1_MA0113.3 NR3C1 MA0113.3 C_NR3C2          960       1.49523      232       1.73894      270     -0.31091  2.34376e-99 4.825678e-98        0.31091
#      total_tfbs.M YM_mean_score YM_bound GO_mean_score GO_bound YM_GO_change YM_GO_pvalue        FDR.M Aging_Change.M
# 1219          960       1.76169      250       1.94587      282      -0.2184  2.24844e-84 8.174626e-84         0.2184

pdf(paste0(Sys.Date(),"_Bindetect_TOBIAS_ZMZ_bulkATAC_BRAIN_FM_Scatter_YvsG.pdf"))
plot(zmz.g.merged$Aging_Change.F, zmz.g.merged$Aging_Change.M, 
     pch = 16, 
     col = adjustcolor( "black", alpha.f = 0.3),
     xlim = c(-0.5,0.5),
     ylim = c(-0.5,0.5),
     xlab = "Differential binding score in Young vs. Old (Females)",
     ylab = "Differential binding score in Young vs. Old (Males)",
     main = "TOBIAS footprint analysis (ZMZ)")
abline(h = 0, col = "grey", lty = "dashed")
abline(v = 0, col = "grey", lty = "dashed")
abline(0 , 1, col = "red", lty = "dashed")
text(-0.5,0.5, paste0("Rho = ",signif(zmz.test.2$estimate,3)), pos = 4)
text(-0.5,0.45, paste0("p = ",signif(zmz.test.2$p.value,3)), pos = 4)
points(zmz.g.merged$Aging_Change.F[nr3c1.zmz.2], zmz.g.merged$Aging_Change.M[nr3c1.zmz.2],
       pch = 16, col = "#CC3333")
text(zmz.g.merged$Aging_Change.F[nr3c1.zmz.2], zmz.g.merged$Aging_Change.M[nr3c1.zmz.2],"Nr3c1",col = "#CC3333", pos = 4)
dev.off()
###################################################################################################################

###################################################################################################################
# 3. Clean table output
options(java.parameters = "-Xmx16g" )
require(openxlsx)

grz.o.merged.v2 <- grz.o.merged[,c("name","motif_id","total_tfbs.F","Aging_Change.F","YF_OF_pvalue", "FDR.F","Aging_Change.M","YM_OM_pvalue", "FDR.M")]
zmz.o.merged.v2 <- zmz.o.merged[,c("name","motif_id","total_tfbs.F","Aging_Change.F","YF_OF_pvalue", "FDR.F","Aging_Change.M","YM_OM_pvalue", "FDR.M")]
zmz.g.merged.v2 <- zmz.g.merged[,c("name","motif_id","total_tfbs.F","Aging_Change.F","YF_GF_pvalue", "FDR.F","Aging_Change.M","YM_GO_pvalue", "FDR.M")]

final.cols <- c("TF_name","JASPAR_motif_id","total_tfbs","Aging_Change.F","F_aging_pvalue", "FDR.F","Aging_Change.M","M_aging_pvalue", "FDR.M")

colnames(grz.o.merged.v2) <- final.cols
colnames(zmz.o.merged.v2) <- final.cols
colnames(zmz.g.merged.v2) <- final.cols

write.xlsx(list("GRZ_Aging_YvO" = grz.o.merged.v2, "ZMZ_Aging_YvO" = zmz.o.merged.v2, "ZMZ_Aging_YvG" = zmz.g.merged.v2), rowNames = F, file = paste0(Sys.Date(),"_GRZ_ZMZ_TOBIAS_ATAC_footprint_Aging_Results.xlsx"))
###################################################################################################################

###################################################################################################################
# 4. bubble plot output

# only significant
grz.tob.Fo.v3 <- grz.tob.Fo[grz.tob.Fo$FDR < 0.05,]
grz.tob.Mo.v3 <- grz.tob.Mo[grz.tob.Mo$FDR < 0.05,]
zmz.tob.Fo.v3 <- zmz.tob.Fo[zmz.tob.Fo$FDR < 0.05,]
zmz.tob.Mo.v3 <- zmz.tob.Mo[zmz.tob.Mo$FDR < 0.05,]
zmz.tob.Fg.v3 <- zmz.tob.Fg[zmz.tob.Fg$FDR < 0.05,]
zmz.tob.Mg.v3 <- zmz.tob.Mg[zmz.tob.Mg$FDR < 0.05,]

# get rankings
grz.tob.Fo.v3$TF_Rank_UP <- rank(-grz.tob.Fo.v3$Aging_Change)  # so it's decreasing order
grz.tob.Mo.v3$TF_Rank_UP <- rank(-grz.tob.Mo.v3$Aging_Change)  # so it's decreasing order
zmz.tob.Fo.v3$TF_Rank_UP <- rank(-zmz.tob.Fo.v3$Aging_Change)  # so it's decreasing order
zmz.tob.Mo.v3$TF_Rank_UP <- rank(-zmz.tob.Mo.v3$Aging_Change)  # so it's decreasing order
zmz.tob.Fg.v3$TF_Rank_UP <- rank(-zmz.tob.Fg.v3$Aging_Change)  # so it's decreasing order
zmz.tob.Mg.v3$TF_Rank_UP <- rank(-zmz.tob.Mg.v3$Aging_Change)  # so it's decreasing order

grz.tob.Fo.v3$TF_Rank_DOWN <- rank(grz.tob.Fo.v3$Aging_Change)  # so it's increasing order
grz.tob.Mo.v3$TF_Rank_DOWN <- rank(grz.tob.Mo.v3$Aging_Change)  # so it's increasing order
zmz.tob.Fo.v3$TF_Rank_DOWN <- rank(zmz.tob.Fo.v3$Aging_Change)  # so it's increasing order
zmz.tob.Mo.v3$TF_Rank_DOWN <- rank(zmz.tob.Mo.v3$Aging_Change)  # so it's increasing order
zmz.tob.Fg.v3$TF_Rank_DOWN <- rank(zmz.tob.Fg.v3$Aging_Change)  # so it's increasing order
zmz.tob.Mg.v3$TF_Rank_DOWN <- rank(zmz.tob.Mg.v3$Aging_Change)  # so it's increasing order

# Rank products
grz.o.merged.v3 <- merge(grz.tob.Fo.v3[,c(1:4,12:15)],grz.tob.Mo.v3[,c(1,12:15)], by = "output_prefix", suffixes = c(".F",".M"))
zmz.o.merged.v3 <- merge(zmz.tob.Fo.v3[,c(1:4,12:15)],zmz.tob.Mo.v3[,c(1,12:15)], by = "output_prefix", suffixes = c(".F",".M"))
zmz.g.merged.v3 <- merge(zmz.tob.Fg.v3[,c(1:4,12:15)],zmz.tob.Mg.v3[,c(1,12:15)], by = "output_prefix", suffixes = c(".F",".M"))

colnames(zmz.g.merged.v3)[-c(1:4)] <- paste0(colnames(zmz.g.merged.v3[,-c(1:4)]) ,".zmz_g")
                                                      
grz.zmz.o   <- merge(grz.o.merged.v3,zmz.o.merged.v3[,-c(2:4)], by = "output_prefix", suffixes = c(".grz_o",".zmz_o"))
grz.zmz.all <- merge(grz.zmz.o, zmz.g.merged.v3[,-c(2:4)]      , by = "output_prefix")

grz.zmz.all$RP_UP   <- grz.zmz.all$TF_Rank_UP.F.grz_o * grz.zmz.all$TF_Rank_UP.M.grz_o * grz.zmz.all$TF_Rank_UP.F.zmz_o * grz.zmz.all$TF_Rank_UP.M.zmz_o * grz.zmz.all$TF_Rank_UP.F.zmz_g * grz.zmz.all$TF_Rank_UP.M.zmz_g
grz.zmz.all$RP_DOWN <- grz.zmz.all$TF_Rank_DOWN.F.grz_o * grz.zmz.all$TF_Rank_DOWN.M.grz_o * grz.zmz.all$TF_Rank_DOWN.F.zmz_o * grz.zmz.all$TF_Rank_DOWN.M.zmz_o * grz.zmz.all$TF_Rank_DOWN.F.zmz_g * grz.zmz.all$TF_Rank_DOWN.M.zmz_g


#### Grab top 10 motifs each way
up.top10  <- sort(grz.zmz.all$RP_UP  , index.return = T)$ix[1:10]
dwn.top10 <- sort(grz.zmz.all$RP_DOWN, index.return = T)$ix[rev(1:10)]

# get filtered/merged datafame for ggplot
grz.zmz.all.top10 <- grz.zmz.all[c(up.top10,dwn.top10),]

## clean up data for ggplot object
my.top.tob.grz.Fo <- grz.zmz.all.top10[,c("output_prefix","FDR.F.grz_o","Aging_Change.F.grz_o")]
my.top.tob.grz.Mo <- grz.zmz.all.top10[,c("output_prefix","FDR.M.grz_o","Aging_Change.M.grz_o")]
my.top.tob.zmz.Fo <- grz.zmz.all.top10[,c("output_prefix","FDR.F.zmz_o","Aging_Change.F.zmz_o")]
my.top.tob.zmz.Mo <- grz.zmz.all.top10[,c("output_prefix","FDR.M.zmz_o","Aging_Change.M.zmz_o")]
my.top.tob.zmz.Fg <- grz.zmz.all.top10[,c("output_prefix","FDR.F.zmz_g","Aging_Change.F.zmz_g")]
my.top.tob.zmz.Mg <- grz.zmz.all.top10[,c("output_prefix","FDR.M.zmz_g","Aging_Change.M.zmz_g")]

colnames(my.top.tob.grz.Fo) <- c("Motif.Name","FDR","Aging_Score")
colnames(my.top.tob.grz.Mo) <- c("Motif.Name","FDR","Aging_Score")
colnames(my.top.tob.zmz.Fo) <- c("Motif.Name","FDR","Aging_Score")
colnames(my.top.tob.zmz.Mo) <- c("Motif.Name","FDR","Aging_Score")
colnames(my.top.tob.zmz.Fg) <- c("Motif.Name","FDR","Aging_Score")
colnames(my.top.tob.zmz.Mg) <- c("Motif.Name","FDR","Aging_Score")

my.top.tob.grz.Fo$Condition <- "YvO_GRZ_F"
my.top.tob.grz.Mo$Condition <- "YvO_GRZ_M"
my.top.tob.zmz.Fo$Condition <- "YvO_ZMZ_F"
my.top.tob.zmz.Mo$Condition <- "YvO_ZMZ_M"
my.top.tob.zmz.Fg$Condition <- "YvG_ZMZ_F"
my.top.tob.zmz.Mg$Condition <- "YvG_ZMZ_M"

my.top.data.merge  <- rbind(my.top.tob.grz.Fo,
                            my.top.tob.grz.Mo,
                            my.top.tob.zmz.Fo,
                            my.top.tob.zmz.Mo,
                            my.top.tob.zmz.Fg,
                            my.top.tob.zmz.Mg)
my.top.data.merge$minusLog10FDR <- -log10(my.top.data.merge$FDR)

#### 
my.max <- max(my.top.data.merge$Aging_Score)
my.min <- min(my.top.data.merge$Aging_Score)
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
my.top.data.merge$Motif.Name  <- factor(my.top.data.merge$Motif.Name, levels = rev(my.top.tob.grz.Fo$Motif.Name))
my.top.data.merge$Condition   <- factor(my.top.data.merge$Condition, levels = c("YvO_GRZ_F","YvO_GRZ_M","YvO_ZMZ_F","YvO_ZMZ_M","YvG_ZMZ_F","YvG_ZMZ_M"))

pdf(paste0(Sys.Date(),"TOBIAS_Top20_Enriched_JASPAR_DA_footprints_ATAC_KilliBrainAging_bothStrains_FDR5.pdf"), height = 5, width=5.5)
my.plot <- ggplot(my.top.data.merge,aes(x=Condition,y=Motif.Name,colour=Aging_Score,size=minusLog10FDR))+ theme_bw() + geom_point(shape = 16)
my.plot <- my.plot + ggtitle("TOBIAS/JASPAR Footprints") + labs(x = "Killifish Brain Aging", y = "Aging DA ATAC footprints (FDR < 5%)")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled)
my.plot <- my.plot + scale_x_discrete(guide = guide_axis(angle = 45)) + scale_size(range = c(3, 10)) + scale_size_area(limits=c(10,180))
print(my.plot)
dev.off()




#### Grab top 5 motifs each way
up.top5  <- sort(grz.zmz.all$RP_UP  , index.return = T)$ix[1:5]
dwn.top5 <- sort(grz.zmz.all$RP_DOWN, index.return = T)$ix[rev(1:5)]

# get filtered/merged datafame for ggplot
grz.zmz.all.top5 <- grz.zmz.all[c(up.top5,dwn.top5),]

## clean up data for ggplot object
my.top.tob.grz.Fo <- grz.zmz.all.top5[,c("output_prefix","FDR.F.grz_o","Aging_Change.F.grz_o")]
my.top.tob.grz.Mo <- grz.zmz.all.top5[,c("output_prefix","FDR.M.grz_o","Aging_Change.M.grz_o")]
my.top.tob.zmz.Fo <- grz.zmz.all.top5[,c("output_prefix","FDR.F.zmz_o","Aging_Change.F.zmz_o")]
my.top.tob.zmz.Mo <- grz.zmz.all.top5[,c("output_prefix","FDR.M.zmz_o","Aging_Change.M.zmz_o")]
my.top.tob.zmz.Fg <- grz.zmz.all.top5[,c("output_prefix","FDR.F.zmz_g","Aging_Change.F.zmz_g")]
my.top.tob.zmz.Mg <- grz.zmz.all.top5[,c("output_prefix","FDR.M.zmz_g","Aging_Change.M.zmz_g")]

colnames(my.top.tob.grz.Fo) <- c("Motif.Name","FDR","Aging_Score")
colnames(my.top.tob.grz.Mo) <- c("Motif.Name","FDR","Aging_Score")
colnames(my.top.tob.zmz.Fo) <- c("Motif.Name","FDR","Aging_Score")
colnames(my.top.tob.zmz.Mo) <- c("Motif.Name","FDR","Aging_Score")
colnames(my.top.tob.zmz.Fg) <- c("Motif.Name","FDR","Aging_Score")
colnames(my.top.tob.zmz.Mg) <- c("Motif.Name","FDR","Aging_Score")

my.top.tob.grz.Fo$Condition <- "YvO_GRZ_F"
my.top.tob.grz.Mo$Condition <- "YvO_GRZ_M"
my.top.tob.zmz.Fo$Condition <- "YvO_ZMZ_F"
my.top.tob.zmz.Mo$Condition <- "YvO_ZMZ_M"
my.top.tob.zmz.Fg$Condition <- "YvG_ZMZ_F"
my.top.tob.zmz.Mg$Condition <- "YvG_ZMZ_M"

my.top.data.merge  <- rbind(my.top.tob.grz.Fo,
                            my.top.tob.grz.Mo,
                            my.top.tob.zmz.Fo,
                            my.top.tob.zmz.Mo,
                            my.top.tob.zmz.Fg,
                            my.top.tob.zmz.Mg)
my.top.data.merge$minusLog10FDR <- -log10(my.top.data.merge$FDR)

#### 
my.max <- max(my.top.data.merge$Aging_Score)
my.min <- min(my.top.data.merge$Aging_Score)
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
my.top.data.merge$Motif.Name  <- factor(my.top.data.merge$Motif.Name, levels = rev(my.top.tob.grz.Fo$Motif.Name))
my.top.data.merge$Condition   <- factor(my.top.data.merge$Condition, levels = c("YvO_GRZ_F","YvO_GRZ_M","YvO_ZMZ_F","YvO_ZMZ_M","YvG_ZMZ_F","YvG_ZMZ_M"))

pdf(paste0(Sys.Date(),"TOBIAS_Top10_Enriched_JASPAR_DA_footprints_ATAC_KilliBrainAging_bothStrains_FDR5.pdf"), height = 4, width=5.5)
my.plot <- ggplot(my.top.data.merge,aes(x=Condition,y=Motif.Name,colour=Aging_Score,size=minusLog10FDR))+ theme_bw() + geom_point(shape = 16)
my.plot <- my.plot + ggtitle("TOBIAS/JASPAR Footprints") + labs(x = "Killifish Brain Aging", y = "Aging DA ATAC footprints (FDR < 5%)")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled)
my.plot <- my.plot + scale_x_discrete(guide = guide_axis(angle = 45)) + scale_size(range = c(3, 10)) + scale_size_area(limits=c(10,180))
print(my.plot)
dev.off()

# to preserve the wanted order
my.top.data.merge$Condition   <- factor(my.top.data.merge$Condition, levels = c("YvO_GRZ_F","YvO_ZMZ_F","YvG_ZMZ_F","YvO_GRZ_M","YvO_ZMZ_M","YvG_ZMZ_M"))

pdf(paste0(Sys.Date(),"TOBIAS_Top10_Enriched_JASPAR_DA_footprints_ATAC_KilliBrainAging_bothStrains_FDR5_reorder.pdf"), height = 4, width=5.5)
my.plot <- ggplot(my.top.data.merge,aes(x=Condition,y=Motif.Name,colour=Aging_Score,size=minusLog10FDR))+ theme_bw() + geom_point(shape = 16)
my.plot <- my.plot + ggtitle("TOBIAS/JASPAR Footprints") + labs(x = "Killifish Brain Aging", y = "Aging DA ATAC footprints (FDR < 5%)")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled)
my.plot <- my.plot + scale_x_discrete(guide = guide_axis(angle = 45)) + scale_size(range = c(3, 10)) + scale_size_area(limits=c(10,180))
print(my.plot)
dev.off()

