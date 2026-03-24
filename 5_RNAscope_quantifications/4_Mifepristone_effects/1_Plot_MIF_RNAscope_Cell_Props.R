setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/Miscroscopy/Mifepristone_Cell_types')

library(readxl)
library(beeswarm)
library(outliers)
library(ARTool)
library(scales)

# 2026-03-24
# Plot cell abundances from RNAscope in the mif supplementation experiment
# Use QPath data from Rapheal

####################################################################
## 0. read and preprocess data

# read data into tibbles
my.data.mif <- read_xlsx('2026-03-24_Mifepristone_GRZ_Cell_Proportions_Data.xlsx', sheet = "ForR")

my.data.mif$Group        <- factor(my.data.mif$Group, levels = c("YF", "OFC", "OFM", "YM","OMC","OMM"))
my.data.mif$Group_no_sex <- factor(my.data.mif$Group_no_sex, levels = c("Young", "Old", "Mifepristone"))

### for ART_ANOVA
my.data.mif$Age <- factor(my.data.mif$Age)
my.data.mif$Sex <- factor(my.data.mif$Sex)
####################################################################


####################################################################
## 1. Plot data

summary(my.data.mif)

#### Apoeb
boxplot( apoeb_av_perc ~ Group,
         data = my.data.mif,
         ylim = c(0,12),
         ylab = "% apoeb cells",
         las = 1,
         col = c("deeppink","deeppink3","deeppink4","deepskyblue" ,"deepskyblue3","deepskyblue4"),
         main = "Microglia percentage")
beeswarm(apoeb_av_perc ~ Group, data = my.data.mif, add = T, pch = 16)

### 1 potential outlier in Mif females: not significant
grubbs.test(my.data.mif$mpz_av_perc[my.data.mif$Group %in% "OFM"])
# G = 1.1061, U = 0.4562, p-value = 0.5251
# alternative hypothesis: highest value 8.32713731851156 is an outlier

apoeb.grz.aov <- summary(aov(apoeb_av_perc ~ Group_no_sex + Sex, data = my.data.mif))
apoeb.grz.aov
#              Df Sum Sq Mean Sq F value   Pr(>F)    
# Group_no_sex  2  81.97   40.98  13.659 0.000346 ***
# Sex           1   3.68    3.68   1.225 0.284730    
# Residuals    16  48.01    3.00                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

shapiro.test(aov(apoeb_av_perc ~ Group_no_sex + Sex, data = my.data.mif)$residuals)
# W = 0.93607, p-value = 0.2019 ### Normality is not violated

# post hoc test
posthoc.grz <- TukeyHSD(aov(apoeb_av_perc ~ Group_no_sex + Sex, data = my.data.mif), which = "Group_no_sex")
posthoc.grz
#   Tukey multiple comparisons of means
#     95% family-wise confidence level
# 
# Fit: aov(formula = apoeb_av_perc ~ Group_no_sex + Sex, data = my.data.mif)
# 
# $Group_no_sex
#                         diff        lwr        upr     p adj
# Old-Young           5.227033  2.6464929  7.8075740 0.0002323
# Mifepristone-Young  2.627462  0.2135877  5.0413370 0.0320163
# Mifepristone-Old   -2.599571 -5.0134458 -0.1856965 0.0339582

pdf(paste0(Sys.Date(),"_GRZ_RNAscope_CellProportions_MIFEPRISTONE.pdf"), height = 5, width = 5.5)
boxplot( apoeb_av_perc ~ Group,
         data = my.data.mif,
         ylim = c(0,12),
         ylab = "% apoeb cells",
         las = 1,
         col = c("deeppink","deeppink3","deeppink4","deepskyblue" ,"deepskyblue3","deepskyblue4"),
         main = "Microglia percentage")
beeswarm(apoeb_av_perc ~ Group, data = my.data.mif, add = T, pch = 16)
text(0.25,12  , paste("Treatment p = ", scientific(apoeb.grz.aov[[1]]$`Pr(>F)`[1],2)), pos = 4)
text(0.25,11.5, paste("Sex p = ", scientific(apoeb.grz.aov[[1]]$`Pr(>F)`[2],2)), pos = 4)

text(4,12  , paste("Y-O p = ", scientific(posthoc.grz$Group_no_sex[,4][1],2)), pos = 4)
text(4,11.5, paste("C-Mif p = ", scientific(posthoc.grz$Group_no_sex[,4][3],2)), pos = 4)
dev.off()
####################################################################



##################
sink(paste0(Sys.Date(),"_MIF_cell_abundance_plotting_sessionInfo.txt"))
sessionInfo()
sink()