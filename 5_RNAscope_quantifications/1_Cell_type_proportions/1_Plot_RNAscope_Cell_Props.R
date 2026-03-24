setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/Miscroscopy/Cell_Type_abundances')

library(readxl)
library(beeswarm)
library(outliers)
library(ARTool)
library(scales)

# 2025-07-18
# Plot cell abundances from RNAscope
# use the average values of 3 Z stacks per animal computed
# by Rapheal


####################################################################
## 0. read and preprocess data

# read data into tibbles
my.data.grz <- read_xlsx('2025-07-cell_abundances_aging_brain_summary.xlsx', sheet = "GRZ")
my.data.zmz <- read_xlsx('2025-07-cell_abundances_aging_brain_summary.xlsx', sheet = "ZMZ")

my.data.grz$Group <- factor(my.data.grz$Group, levels = c("YF", "MF", "OF", "YM","MM","OM"))
my.data.zmz$Group <- factor(my.data.zmz$Group, levels = c("YF", "OF", "GF", "YM","OM","GM"))

### for ART_ANOVA
my.data.grz$Age <- factor(my.data.grz$Age)
my.data.grz$Sex <- factor(my.data.grz$Sex)

my.data.zmz$Age <- factor(my.data.zmz$Age)
my.data.zmz$Sex <- factor(my.data.zmz$Sex)
####################################################################


####################################################################
## 1. Plot GRZ data

summary(my.data.grz)

#### Apoeb
boxplot( apoeb_av_perc ~ Group,
         data = my.data.grz,
         ylim = c(0,25),
         ylab = "% apoeb cells",
         las = 1,
         col = c("deeppink","deeppink3","deeppink4","deepskyblue" ,"deepskyblue3","deepskyblue4"),
         main = "Microglia percentage (GRZ)")
beeswarm(apoeb_av_perc ~ Group, data = my.data.grz, add = T, pch = 16)

apoeb.grz.aov <- summary(aov(apoeb_av_perc ~ Age + Sex, data = my.data.grz))
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# Age          2  400.2  200.11  20.564 6.27e-06 ***
# Sex          1   66.8   66.77   6.862    0.015 *  
# Residuals   24  233.5    9.73                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

shapiro.test(aov(apoeb_av_perc ~ Age + Sex, data = my.data.grz)$residuals)
# W = 0.9798, p-value = 0.8456  ### Normality is not violated


#### Map2
boxplot( map2_av_perc ~ Group,
         data = my.data.grz,
         ylim = c(0,70),
         ylab = "% map2 cells",
         las = 1,
         col = c("deeppink","deeppink3","deeppink4","deepskyblue" ,"deepskyblue3","deepskyblue4"),
         main = "Neuron percentage (GRZ)")
beeswarm(map2_av_perc ~ Group, data = my.data.grz, add = T, pch = 16)

map2.grz.aov <- summary(aov(map2_av_perc ~ Age + Sex, data = my.data.grz))
#             Df Sum Sq Mean Sq F value Pr(>F)
# Age          2     37    18.6   0.115  0.892
# Sex          1    322   322.1   1.990  0.171
# Residuals   24   3884   161.8               

shapiro.test(aov(map2_av_perc ~ Age + Sex, data = my.data.grz)$residuals)
# W = 0.97597, p-value = 0.7453  ### Normality is not violated


#### s100b
boxplot( s100b_av_perc ~ Group,
         data = my.data.grz,
         ylim = c(0,60),
         ylab = "% s100b cells",
         las = 1,
         col = c("deeppink","deeppink3","deeppink4","deepskyblue" ,"deepskyblue3","deepskyblue4"),
         main = "Astrocyte percentage (GRZ)")
beeswarm(s100b_av_perc ~ Group, data = my.data.grz, add = T, pch = 16)

s100b.grz.aov <- summary(aov(s100b_av_perc ~ Age + Sex, data = my.data.grz))
#             Df Sum Sq Mean Sq F value Pr(>F)  
# Age          2    926   462.8   3.014 0.0680 .
# Sex          1    469   468.9   3.053 0.0934 .
# Residuals   24   3686   153.6                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

shapiro.test(aov(s100b_av_perc ~ Age + Sex, data = my.data.grz)$residuals)
# W = 0.94979, p-value = 0.1957 ### Normality is not violated


#### mpz
boxplot( mpz_av_perc ~ Group,
         data = my.data.grz,
         ylim = c(0,60),
         ylab = "% mpz cells",
         las = 1,
         col = c("deeppink","deeppink3","deeppink4","deepskyblue" ,"deepskyblue3","deepskyblue4"),
         main = "Oligodendrocyte percentage (GRZ)")
beeswarm(mpz_av_perc ~ Group, data = my.data.grz, add = T, pch = 16)

### 1 potential outlier in Old males
grubbs.test(my.data.grz$mpz_av_perc[my.data.grz$Group %in% "OM"])
# G = 1.77420, U = 0.24453, p-value = 0.07366
# alternative hypothesis: highest value 52.377417304955 is an outlier

mpz.grz.aov <- summary(aov(mpz_av_perc ~ Age + Sex, data = my.data.grz))
#             Df Sum Sq Mean Sq F value Pr(>F)
# Age          2  136.2    68.1   0.769  0.475
# Sex          1  162.5   162.5   1.834  0.188
# Residuals   24 2126.4    88.6               

shapiro.test(aov(mpz_av_perc ~ Age + Sex, data = my.data.grz)$residuals)
# W = 0.89603, p-value = 0.009245 ### Normality IS violated

# Transform the data to ranks and compute ART ANOVA
mpz.tr.dat  <- art(mpz_av_perc ~ Age*Sex , data=my.data.grz)
mpz.aov.res <- anova(mpz.tr.dat)
mpz.aov.res
#           Df Df.res F value  Pr(>F)  
# 1 Age      2     22 0.33187 0.72111  
# 2 Sex      1     22 1.41110 0.24754  
# 3 Age:Sex  2     22 0.18910 0.82904  
# ---
# Signif. codes:   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 


pdf(paste0(Sys.Date(),"_GRZ_RNAscope_CellProportions.pdf"), height = 4.5, width = 11)
par(mfrow = c(1,4))
boxplot( apoeb_av_perc ~ Group,
         data = my.data.grz,
         ylim = c(0,100),
         ylab = "% apoeb cells",
         las = 2,
         col = c("deeppink","deeppink3","deeppink4","deepskyblue" ,"deepskyblue3","deepskyblue4"),
         main = "Microglia (GRZ)")
beeswarm(apoeb_av_perc ~ Group, data = my.data.grz, add = T, pch = 16, cex = 1.5)
text(0.5,100, paste("Age p = ", scientific(apoeb.grz.aov[[1]]$`Pr(>F)`[1],2)), pos = 4)
text(0.5,90, paste("Sex p = ", scientific(apoeb.grz.aov[[1]]$`Pr(>F)`[2],2)), pos = 4)

boxplot( map2_av_perc ~ Group,
         data = my.data.grz,
         ylim = c(0,100),
         ylab = "% map2 cells",
         las = 2,
         col = c("deeppink","deeppink3","deeppink4","deepskyblue" ,"deepskyblue3","deepskyblue4"),
         main = "Neuron (GRZ)")
beeswarm(map2_av_perc ~ Group, data = my.data.grz, add = T, pch = 16, cex = 1.5)
text(0.5,100, paste("Age p = ", scientific(map2.grz.aov[[1]]$`Pr(>F)`[1],2)), pos = 4)
text(0.5,90, paste("Sex p = ", scientific(map2.grz.aov[[1]]$`Pr(>F)`[2],2)), pos = 4)

boxplot( s100b_av_perc ~ Group,
         data = my.data.grz,
         ylim = c(0,100),
         ylab = "% s100b cells",
         las = 2,
         col = c("deeppink","deeppink3","deeppink4","deepskyblue" ,"deepskyblue3","deepskyblue4"),
         main = "Astrocyte (GRZ)")
beeswarm(s100b_av_perc ~ Group, data = my.data.grz, add = T, pch = 16, cex = 1.5)
text(0.5,100, paste("Age p = ", scientific(s100b.grz.aov[[1]]$`Pr(>F)`[1],2)), pos = 4)
text(0.5,90, paste("Sex p = ", scientific(s100b.grz.aov[[1]]$`Pr(>F)`[2],2)), pos = 4)

boxplot( mpz_av_perc ~ Group,
         data = my.data.grz,
         ylim = c(0,100),
         ylab = "% mpz cells",
         las = 2,
         col = c("deeppink","deeppink3","deeppink4","deepskyblue" ,"deepskyblue3","deepskyblue4"),
         main = "Oligodendrocyte (GRZ)")
beeswarm(mpz_av_perc ~ Group, data = my.data.grz, add = T, pch = 16, cex = 1.5)
text(0.5,100, paste("Age p = ", scientific(mpz.aov.res$`Pr(>F)`[1],2)), pos = 4)
text(0.5,90, paste("Sex p = ", scientific(mpz.aov.res$`Pr(>F)`[2],2)), pos = 4)
par(mfrow = c(1,1))
dev.off()
####################################################################



####################################################################
## 2. Plot ZMZ data

summary(my.data.zmz)

#### Apoeb
boxplot( apoeb_av_perc ~ Group,
         data = my.data.zmz,
         ylim = c(0,25),
         ylab = "% apoeb cells",
         las = 1,
         col = c("deeppink","deeppink4","magenta4","deepskyblue" ,"deepskyblue4","royalblue4"),
         main = "Microglia percentage (zmz)")
beeswarm(apoeb_av_perc ~ Group, data = my.data.zmz, add = T, pch = 16)

apoeb.zmz.aov <- summary(aov(apoeb_av_perc ~ Age + Sex, data = my.data.zmz))
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# Age          2  570.0  284.98  21.457 2.64e-06 ***
# Sex          1    4.9    4.88   0.367    0.549    
# Residuals   27  358.6   13.28                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

shapiro.test(aov(apoeb_av_perc ~ Age + Sex, data = my.data.zmz)$residuals)
# W = 0.9798, p-value = 0.8456  ### Normality IS violated

# Transform the data to ranks and compute ART ANOVA
apoeb.tr.dat  <- art(apoeb_av_perc ~ Age*Sex , data=my.data.zmz)
apoeb.aov.res <- anova(apoeb.tr.dat)
apoeb.aov.res
#           Df Df.res  F value     Pr(>F)    
# 1 Age      2     25 24.11022 1.4667e-06 ***
# 2 Sex      1     25  0.68805    0.41468    
# 3 Age:Sex  2     25  0.73781    0.48828   


#### Map2
boxplot( map2_av_perc ~ Group,
         data = my.data.zmz,
         ylim = c(0,70),
         ylab = "% map2 cells",
         las = 1,
         col = c("deeppink","deeppink4","magenta4","deepskyblue" ,"deepskyblue4","royalblue4"),
         main = "Neuron percentage (zmz)")
beeswarm(map2_av_perc ~ Group, data = my.data.zmz, add = T, pch = 16)

map2.zmz.aov <- summary(aov(map2_av_perc ~ Age + Sex, data = my.data.zmz))
#             Df Sum Sq Mean Sq F value Pr(>F)  
# Age          2  977.2   488.6   4.430 0.0217 *
# Sex          1   12.0    12.0   0.109 0.7438  
# Residuals   27 2977.9   110.3                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

shapiro.test(aov(map2_av_perc ~ Age + Sex, data = my.data.zmz)$residuals)
# W = 0.96245, p-value = 0.3382 ### Normality is not violated


#### s100b
boxplot( s100b_av_perc ~ Group,
         data = my.data.zmz,
         ylim = c(0,60),
         ylab = "% s100b cells",
         las = 1,
         col = c("deeppink","deeppink4","magenta4","deepskyblue" ,"deepskyblue4","royalblue4"),
         main = "Astrocyte percentage (zmz)")
beeswarm(s100b_av_perc ~ Group, data = my.data.zmz, add = T, pch = 16)

s100b.zmz.aov <- summary(aov(s100b_av_perc ~ Age + Sex, data = my.data.zmz))
#             Df Sum Sq Mean Sq F value Pr(>F)
# Age          2  310.8  155.40   2.016  0.153
# Sex          1   44.3   44.29   0.574  0.455
# Residuals   27 2081.6   77.10               

shapiro.test(aov(s100b_av_perc ~ Age + Sex, data = my.data.zmz)$residuals)
# W = 0.93247, p-value = 0.05122 ### Normality is not violated


#### mpz
boxplot( mpz_av_perc ~ Group,
         data = my.data.zmz,
         ylim = c(0,60),
         ylab = "% mpz cells",
         las = 1,
         col = c("deeppink","deeppink4","magenta4","deepskyblue" ,"deepskyblue4","royalblue4"),
         main = "Oligodendrocyte percentage (zmz)")
beeswarm(mpz_av_perc ~ Group, data = my.data.zmz, add = T, pch = 16)


mpz.zmz.aov <- summary(aov(mpz_av_perc ~ Age + Sex, data = my.data.zmz))
#             Df Sum Sq Mean Sq F value Pr(>F)  
# Age          2  202.5  101.24   3.970 0.0308 *
# Sex          1    1.4    1.37   0.054 0.8187  
# Residuals   27  688.5   25.50                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

shapiro.test(aov(mpz_av_perc ~ Age + Sex, data = my.data.zmz)$residuals)
#W = 0.97756, p-value = 0.742  ### Normality is not violated


pdf(paste0(Sys.Date(),"_ZMZ1001_RNAscope_CellProportions.pdf"), height = 4.5, width = 11)
par(mfrow = c(1,4))
boxplot( apoeb_av_perc ~ Group,
         data = my.data.zmz,
         ylim = c(0,100),
         ylab = "% apoeb cells",
         las = 2,
         col = c("deeppink","deeppink4","magenta4","deepskyblue" ,"deepskyblue4","royalblue4"),
         main = "Microglia (zmz)")
beeswarm(apoeb_av_perc ~ Group, data = my.data.zmz, add = T, pch = 16, cex = 1.5)
text(0.5,100, paste("Age p = ", scientific(apoeb.aov.res$`Pr(>F)`[1],2)), pos = 4)
text(0.5,90, paste("Sex p = ", scientific(apoeb.aov.res$`Pr(>F)`[2],2)), pos = 4)

boxplot( map2_av_perc ~ Group,
         data = my.data.zmz,
         ylim = c(0,100),
         ylab = "% map2 cells",
         las = 2,
         col = c("deeppink","deeppink4","magenta4","deepskyblue" ,"deepskyblue4","royalblue4"),
         main = "Neuron (zmz)")
beeswarm(map2_av_perc ~ Group, data = my.data.zmz, add = T, pch = 16, cex = 1.5)
text(0.5,100, paste("Age p = ", scientific(map2.zmz.aov[[1]]$`Pr(>F)`[1],2)), pos = 4)
text(0.5,90, paste("Sex p = ", scientific(map2.zmz.aov[[1]]$`Pr(>F)`[2],2)), pos = 4)

boxplot( s100b_av_perc ~ Group,
         data = my.data.zmz,
         ylim = c(0,100),
         ylab = "% s100b cells",
         las = 2,
         col = c("deeppink","deeppink4","magenta4","deepskyblue" ,"deepskyblue4","royalblue4"),
         main = "Astrocyte (zmz)")
beeswarm(s100b_av_perc ~ Group, data = my.data.zmz, add = T, pch = 16, cex = 1.5)
text(0.5,100, paste("Age p = ", scientific(s100b.zmz.aov[[1]]$`Pr(>F)`[1],2)), pos = 4)
text(0.5,90, paste("Sex p = ", scientific(s100b.zmz.aov[[1]]$`Pr(>F)`[2],2)), pos = 4)

boxplot( mpz_av_perc ~ Group,
         data = my.data.zmz,
         ylim = c(0,100),
         ylab = "% mpz cells",
         las = 2,
         col = c("deeppink","deeppink4","magenta4","deepskyblue" ,"deepskyblue4","royalblue4"),
         main = "Oligodendrocyte (zmz)")
beeswarm(mpz_av_perc ~ Group, data = my.data.zmz, add = T, pch = 16, cex = 1.5)
text(0.5,100, paste("Age p = ", scientific(mpz.zmz.aov[[1]]$`Pr(>F)`[1],2)), pos = 4)
text(0.5,90, paste("Sex p = ", scientific(mpz.zmz.aov[[1]]$`Pr(>F)`[2],2)), pos = 4)
par(mfrow = c(1,1))
dev.off()
##########################################################################################

##################
sink(paste0(Sys.Date(),"_cell_abundance_plotting_sessionInfo.txt"))
sessionInfo()
sink()