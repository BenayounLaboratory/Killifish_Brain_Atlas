setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/GR_signaling//Cortisol')
options(stringsAsFactors = F)

library(beeswarm)
library(readxl)

library(EnvStats)
library(ARTool)   # ART ANOVA


# 2024-10-11
# plot aging cortisol kidney data

# 2025-06-12
# add new data (seems like average values are very different)

############################################################################################
my.C.data <- data.frame(read_xlsx("2025-06-12_Combined_cortisol_data.xlsx"))
head(my.C.data)
# Sample Shorthand Strain Group Cortisol_pg_per_mg Sex Age_Group
# 1 5w, g SnRNA F1      YGF1    GRZ    YF                184   F         Y
# 2 5w, g SnRNA F2      YGF2    GRZ    YF               2893   F         Y
# 3 5w, g SnRNA F3      YGF3    GRZ    YF                364   F         Y
# 4 5w, g SnRNA F4      YGF4    GRZ    YF                557   F         Y
# 5 5w, g SnRNA F6      YGF6    GRZ    YF                406   F         Y
# 6 5w, g SnRNA M1      YGM1    GRZ    YM                143   M         Y

my.C.data$Sex <-  factor(my.C.data$Sex, levels = c("F","M"))

### subset and factorize
grz.C.data <- my.C.data[my.C.data$Strain %in% "GRZ",]
zmz.C.data <- my.C.data[my.C.data$Strain %in% "ZMZ1001",]

grz.C.data$Group <- factor(grz.C.data$Group, levels = c("YF","OF","YM","OM"))
zmz.C.data$Group <- factor(zmz.C.data$Group, levels = c("YF","OF","GF","YM","OM","GM"))

grz.C.data$Age_Group <- factor(grz.C.data$Age_Group, levels = c("Y","O"))
zmz.C.data$Age_Group <- factor(zmz.C.data$Age_Group, levels = c("Y","O", "G"))

###################################################################################################################
#### A. GRZ
boxplot(Cortisol_pg_per_mg ~ Group,
        data = grz.C.data,
        ylim = c(0,50000),
        ylab = "Cortisol (pg/mg kidney tissue)",
        outline = F,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 1, xlab = "")
beeswarm(Cortisol_pg_per_mg ~ Group,
         data = grz.C.data,
         add = T, pch = 16, cex = 1, corral = "wrap")


########## there seems to be 1 outlier in YF and 2 in YM based on boxplot
rosnerTest(grz.C.data$Cortisol_pg_per_mg[grz.C.data$Group %in% "YF"], k = 1)
# Number of Outliers Detected:     2
# i    Mean.i     SD.i Value Obs.Num    R.i+1 lambda.i+1 Outlier
# 1 0 3423.603 3914.639 16605      16 3.367206   2.619964    TRUE
#### ygf11

rosnerTest(grz.C.data$Cortisol_pg_per_mg[grz.C.data$Group %in% "YM"], k = 2)
# Number of Outliers Detected:     2
# i    Mean.i     SD.i Value Obs.Num    R.i+1 lambda.i+1 Outlier
# 1 0 3232.307 5074.554 22185.776      15 3.735002   2.651599    TRUE
# 2 1 2117.397 1894.124  7619.112      14 2.904623   2.619964    TRUE

# grz.C.data[grz.C.data$Group %in% "YM",]
### ygm8, ygm7

### remove outliers
grz.C.data.cl <- grz.C.data[!(grz.C.data$Shorthand %in% c("ygf11","ygm8", "ygm7")),]


### Test if age-related increase (1 sided), based on hypothesis from omics
grz.f.age   <- wilcox.test(grz.C.data.cl$Cortisol_pg_per_mg[grz.C.data.cl$Group == "YF"],  grz.C.data.cl$Cortisol_pg_per_mg[grz.C.data.cl$Group == "OF"], alternative = "less") 
# p-value = 0.006781
grz.m.age   <- wilcox.test(grz.C.data.cl$Cortisol_pg_per_mg[grz.C.data.cl$Group == "YM"],  grz.C.data.cl$Cortisol_pg_per_mg[grz.C.data.cl$Group == "OM"], alternative = "less") 
# p-value = 0.07549

# summary(aov(Cortisol_pg_per_mg ~ Sex + Age_Group, data = grz.C.data.cl))
# #             Df    Sum Sq   Mean Sq F value   Pr(>F)    
# # Sex          1 2.802e+08 280193088   7.783 0.007058 ** 
# # Age_Group    1 4.953e+08 495309543  13.759 0.000456 ***
# # Residuals   60 2.160e+09  35999907                     
# # ---
# # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# shapiro.test(aov(Cortisol_pg_per_mg ~ Sex + Age_Group, data = grz.C.data.cl)$residuals)
# # 
# # Shapiro-Wilk normality test
# # data:  aov(surface ~ Sex + Age, data = my.tab.data)$residuals
# # W = 0.95801, p-value =4.672e-05


### ART-ANOVA #### https://www.geeksforgeeks.org/what-is-the-non-parametric-equivalent-of-a-two-way-anova-in-r/
# Transform the data
mod.cort.grz <- art(Cortisol_pg_per_mg ~ Sex * Age_Group, data=grz.C.data.cl)

# Conduct ANOVA on the ART model
an.grz <- anova(mod.cort.grz)
# Analysis of Variance of Aligned Rank Transformed Data
# Table Type: Anova Table (Type III tests) 
# Model: No Repeated Measures (lm)
# Response: art(Cortisol_pg_per_mg)
# 
# Df Df.res F value     Pr(>F)    
# 1 Sex            1     59  12.919 0.00066467 ***
# 2 Age_Group      1     59  18.085  7.643e-05 ***
# 3 Sex:Age_Group  1     59   9.608 0.00296767  **


pdf(paste0(Sys.Date(),"_Kidney_Cortisol_GRZ_3cohorts_outliers_removed.pdf"), width = 3, height = 5)
boxplot(Cortisol_pg_per_mg ~ Group,
        data = grz.C.data.cl,
        ylim = c(0,40000),
        ylab = "Cortisol (pg/mg kidney tissue)",
        outline = F,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 1, xlab = "",
        main = "Cortisol (GRZ)")
beeswarm(Cortisol_pg_per_mg ~ Group,
         data = grz.C.data.cl,
         add = T, pch = 16, cex = 1, corral = "wrap")
text(1.5, 38000, signif(grz.f.age$p.value,3))
text(3.5, 38000, signif(grz.m.age$p.value,3))
dev.off()


pdf(paste0(Sys.Date(),"_Kidney_Cortisol_GRZ_3cohorts_outliers_removed_ART_ANOVA.pdf"), width = 3, height = 5)
boxplot(Cortisol_pg_per_mg ~ Group,
        data = grz.C.data.cl,
        ylim = c(0,40000),
        ylab = "Cortisol (pg/mg kidney tissue)",
        outline = F,
        col = c("deeppink","deeppink4","deepskyblue","deepskyblue4"),
        las = 1, xlab = "",
        main = "Cortisol (GRZ)")
beeswarm(Cortisol_pg_per_mg ~ Group,
         data = grz.C.data.cl,
         add = T, pch = 16, cex = 1, corral = "wrap")
text(2.5, 40000, paste("Age group pval ~ ",signif(an.grz$`Pr(>F)`[2],3)))
text(2.5, 37000, paste("Sex pval ~ ",signif(an.grz$`Pr(>F)`[1],3)))
dev.off()



###################################################################################################################
#### B. ZMZ

boxplot(Cortisol_pg_per_mg ~ Group,
        data = zmz.C.data,
        ylim = c(0,50000),
        ylab = "Cortisol (pg/mg kidney tissue)",
        outline = F,
        col = c("deeppink","deeppink4","magenta4","deepskyblue","deepskyblue4","royalblue4"),
        las = 1, xlab = "")
beeswarm(Cortisol_pg_per_mg ~ Group,
         data = zmz.C.data,
         add = T, pch = 16, cex = 1, corral = "wrap")


########## there seems to be 1 outlier in GF based on boxplot
rosnerTest(zmz.C.data$Cortisol_pg_per_mg[zmz.C.data$Group %in% "GF"], k = 1)
# Number of Outliers Detected:     2
# i    Mean.i     SD.i Value Obs.Num    R.i+1 lambda.i+1 Outlier
# 1 0 14338.29 14143.28 46173.61       9 2.250914   2.215004    TRUE
#### GZF1, "26w, Z RNAish F1"

### remove outlier
zmz.C.data.cl <- zmz.C.data[!(zmz.C.data$Sample %in% c("26w, Z RNAish F1")),]

### Test if age-related increase (1 sided), based on hypothesis from omics
zmz.f.YO   <- wilcox.test(zmz.C.data.cl$Cortisol_pg_per_mg[zmz.C.data.cl$Group == "YF"],  zmz.C.data.cl$Cortisol_pg_per_mg[zmz.C.data.cl$Group == "OF"], alternative = "less") 
# p-value = 0.4091
zmz.f.YG   <- wilcox.test(zmz.C.data.cl$Cortisol_pg_per_mg[zmz.C.data.cl$Group == "YF"],  zmz.C.data.cl$Cortisol_pg_per_mg[zmz.C.data.cl$Group == "GF"], alternative = "less") 
# p-value = 0.006327

zmz.m.YO   <- wilcox.test(zmz.C.data.cl$Cortisol_pg_per_mg[zmz.C.data.cl$Group == "YM"],  zmz.C.data.cl$Cortisol_pg_per_mg[zmz.C.data.cl$Group == "OM"], alternative = "less") 
# p-value = 0.5314
zmz.m.YG   <- wilcox.test(zmz.C.data.cl$Cortisol_pg_per_mg[zmz.C.data.cl$Group == "YM"],  zmz.C.data.cl$Cortisol_pg_per_mg[zmz.C.data.cl$Group == "GM"], alternative = "less") 
# p-value = 0.04654


### ART-ANOVA #### https://www.geeksforgeeks.org/what-is-the-non-parametric-equivalent-of-a-two-way-anova-in-r/
# Transform the data
mod.cort.zmz <- art(Cortisol_pg_per_mg ~ Sex * Age_Group, data=zmz.C.data.cl)

# Conduct ANOVA on the ART model
an.zmz <- anova(mod.cort.zmz)
# Analysis of Variance of Aligned Rank Transformed Data
# Table Type: Anova Table (Type III tests) 
# Model: No Repeated Measures (lm)
# Response: art(Cortisol_pg_per_mg)
# 
# Df Df.res F value     Pr(>F)    
# 1 Sex            1     32  3.0111    0.09232   .
# 2 Age_Group      2     32 13.7459 4.9102e-05 ***
# 3 Sex:Age_Group  2     32  1.7149    0.19611    

pdf(paste0(Sys.Date(),"_Kidney_Cortisol_ZMZ1001_outlier_removed.pdf"), width = 4.25, height = 5)
boxplot(Cortisol_pg_per_mg ~ Group,
        data = zmz.C.data.cl,
        ylim = c(0,25000),
        ylab = "Cortisol (pg/mg kidney tissue)",
        outline = F,
        col = c("deeppink","deeppink4","magenta4","deepskyblue","deepskyblue4","royalblue4"),
        las = 1, xlab = "",
        main = "Cortisol (ZMZ1001)")
beeswarm(Cortisol_pg_per_mg ~ Group,
         data = zmz.C.data.cl,
         add = T, pch = 16, cex = 1, corral = "wrap")
text(1.5, 18000, signif(zmz.f.YO$p.value,3))
text(2  , 24500, signif(zmz.f.YG$p.value,3))
text(4.5, 18000, signif(zmz.m.YO$p.value,3))
text(5  , 24500, signif(zmz.m.YG$p.value,3))
dev.off()

pdf(paste0(Sys.Date(),"_Kidney_Cortisol_ZMZ1001_outlier_removed_ART_ANOVA.pdf"), width = 4.25, height = 5)
boxplot(Cortisol_pg_per_mg ~ Group,
        data = zmz.C.data.cl,
        ylim = c(0,25000),
        ylab = "Cortisol (pg/mg kidney tissue)",
        outline = F,
        col = c("deeppink","deeppink4","magenta4","deepskyblue","deepskyblue4","royalblue4"),
        las = 1, xlab = "",
        main = "Cortisol (ZMZ1001)")
beeswarm(Cortisol_pg_per_mg ~ Group,
         data = zmz.C.data.cl,
         add = T, pch = 16, cex = 1, corral = "wrap")
text(2.5, 25000, paste("Age group pval ~ ",signif(an.zmz$`Pr(>F)`[2],3)))
text(2.5, 23000, paste("Sex pval ~ ",signif(an.zmz$`Pr(>F)`[1],3)))
dev.off()

