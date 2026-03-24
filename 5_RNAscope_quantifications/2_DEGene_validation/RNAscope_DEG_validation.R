setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/Miscroscopy/Recurrent_DE')

library(readxl)

library(ComplexHeatmap)   #
library(circlize)         #

library(ARTool)


# 2026-01-23
# Rapheal DAPI normalized RNAscope expression data
# plot and evaluate

####################################################################################################
## 0. read and preprocess data

# read data into tibbles
my.data.grz <- read_xlsx('DEsets1and2_QuPath_BothStrains_121925.xlsx', sheet = "GRZ_MasterSheet")
my.data.zmz <- read_xlsx('DEsets1and2_QuPath_BothStrains_121925.xlsx', sheet = "ZMZ_MasterSheet")

####################################################################################################

####################################################################################################
## 1. Plot normalized data as heatmaps

### format for complexHeatmap
grz.rnascope <- t(my.data.grz[,6:11])
zmz.rnascope <- t(my.data.zmz[,6:11])

colnames(grz.rnascope) <- my.data.grz$Animal
colnames(zmz.rnascope) <- my.data.zmz$Animal

grz.rnascope <- grz.rnascope[sort(rownames(grz.rnascope)),]
zmz.rnascope <- zmz.rnascope[sort(rownames(grz.rnascope)),]


grz.heat <- Heatmap(t(scale(t(grz.rnascope))), 
                    col = colorRamp2(c(-2,0, 2), c("darkblue", "white", "#CC3333"), transparency = 0, space = "LAB"),
                    border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
                    cluster_rows = F,
                    cluster_columns = F,
                    column_title = "GRZ RNAscope",
                    width  = ncol(grz.rnascope)*unit(5, "mm"),
                    height = nrow(grz.rnascope)*unit(5, "mm"))

zmz.heat <- Heatmap(t(scale(t(zmz.rnascope))), 
                    col = colorRamp2(c(-2,0, 2), c("darkblue", "white", "#CC3333"), transparency = 0, space = "LAB"),
                    border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
                    cluster_rows = F,
                    cluster_columns = F,
                    column_title = "ZMZ-1001 RNAscope",
                    width  = ncol(zmz.rnascope)*unit(5, "mm"),
                    height = nrow(zmz.rnascope)*unit(5, "mm"))

pdf(paste0(Sys.Date(),"_aging_GRZ_ZMZ_RNAscope_NormDAPI_heatmap.pdf"), width = 15, height = 5)
print(grz.heat + zmz.heat)
dev.off()
####################################################################################################




####################################################################################################
## 2. Perform differential analysis

################################################################################
####################      ANOVA (normality assumption)      ####################
################################################################################
#### a. GRZ

# update colnames
colnames(my.data.grz)[-c(1:5)] <- gsub(":","_",colnames(my.data.grz)[-c(1:5)])

grz.rnascope.res <- data.frame(matrix(0,6,3))
colnames(grz.rnascope.res) <- c("ANOVA_Age_p","ANOVA_Sex_p","Shapiro_Residuals_p")
rownames(grz.rnascope.res) <- colnames(my.data.grz)[-c(1:5)]

for (i in 1:6) {
  
  # run ANOVA
  grz.aov <- summary(aov(as.formula(paste(colnames(my.data.grz)[i +5], "~ Sex + Age")), data = my.data.grz))
  
  # grab sex and age p-values
  grz.rnascope.res$ANOVA_Age_p[i] <-  grz.aov[[1]]$`Pr(>F)`[2]
  grz.rnascope.res$ANOVA_Sex_p[i] <-  grz.aov[[1]]$`Pr(>F)`[1]
  
  # verify normality
  sw.test <- shapiro.test(aov(as.formula(paste(colnames(my.data.grz)[i +5], "~ Sex + Age")), data = my.data.grz)$residuals)
  
  # grab Shapiro p-values
  grz.rnascope.res$Shapiro_Residuals_p[i] <-  sw.test$p.value
  
}

grz.rnascope.res
#              ANOVA_Age_p ANOVA_Sex_p Shapiro_Residuals_p
# znfx1_DAPI  2.326232e-05  0.09842118           0.2495819
# helz2_DAPI  2.001115e-01  0.48500001           0.8312720   --- n.s on Age
# ifih1_DAPI  6.021856e-03  0.09883077           0.3725453
# optn_DAPI   2.028193e-02  0.54928793           0.7998676
# irf3_DAPI   6.609521e-01  0.25528832           0.1089483   --- n.s on Age
# stat1a_DAPI 9.156166e-03  0.46236936           0.8122536


#### b. ZMZ-1001

# update colnames
colnames(my.data.zmz)[-c(1:5)] <- gsub(":","_",colnames(my.data.zmz)[-c(1:5)])

zmz.rnascope.res <- data.frame(matrix(0,6,3))
colnames(zmz.rnascope.res) <- c("ANOVA_Age_p","ANOVA_Sex_p","Shapiro_Residuals_p")
rownames(zmz.rnascope.res) <- colnames(my.data.zmz)[-c(1:5)]

for (i in 1:6) {
  
  # run ANOVA
  zmz.aov <- summary(aov(as.formula(paste(colnames(my.data.zmz)[i +5], "~ Sex + Age")), data = my.data.zmz))
  
  # grab sex and age p-values
  zmz.rnascope.res$ANOVA_Age_p[i] <-  zmz.aov[[1]]$`Pr(>F)`[2]
  zmz.rnascope.res$ANOVA_Sex_p[i] <-  zmz.aov[[1]]$`Pr(>F)`[1]
  
  # verify normality
  sw.test <- shapiro.test(aov(as.formula(paste(colnames(my.data.zmz)[i +5], "~ Sex + Age")), data = my.data.zmz)$residuals)
  
  # grab Shapiro p-values
  zmz.rnascope.res$Shapiro_Residuals_p[i] <-  sw.test$p.value
  
}

zmz.rnascope.res
#              ANOVA_Age_p ANOVA_Sex_p Shapiro_Residuals_p
# znfx1_DAPI  0.0121730755  0.26807777         0.459299268
# helz2_DAPI  0.0016952832  0.27797461         0.164751784
# ifih1_DAPI  0.0109993245  0.42296640         0.074826730
# optn_DAPI   0.0486361486  0.09154183         0.070413987
# irf3_DAPI   0.0003713119  0.06839633         0.004766269   -- NOT normal
# stat1a_DAPI 0.0141281242  0.14268416         0.001788722   -- NOT normal
################################################################################


################################################################################
####################   ART-ANOVA (NO normality assumption)  ####################
################################################################################
#### a. GRZ

# update colnames
colnames(my.data.grz)[-c(1:5)] <- gsub(":","_",colnames(my.data.grz)[-c(1:5)])

# make factors to ART-ANOVA
my.data.grz$Age <- factor(my.data.grz$Age)
my.data.grz$Sex <- factor(my.data.grz$Sex)

# prepare output data frame
grz.rnascope.res <- data.frame(matrix(0,6,3))
colnames(grz.rnascope.res) <- c("ARTANOVA_Age_p","ARTANOVA_Sex_p","ARTANOVA_Inter_p")
rownames(grz.rnascope.res) <- colnames(my.data.grz)[-c(1:5)]

for (i in 1:6) {
  
  # Transform the data to ranks to compute ART ANOVA
  rnasc.tr.dat  <- art(as.formula(paste(colnames(my.data.grz)[i +5], "~ Sex * Age")) , data=my.data.grz)
  rnasc.aov.res <- anova(rnasc.tr.dat)
  # rnasc.aov.res
  #          Df Df.res F value     Pr(>F)    
  # 1 Sex      1     16  3.7736   0.069881   .
  # 2 Age      1     16 48.5437 3.1727e-06 ***
  # 3 Sex:Age  1     16  6.7086   0.019734   *
  
  
  # grab sex and age p-values
  grz.rnascope.res$ARTANOVA_Age_p[i]   <-  rnasc.aov.res$`Pr(>F)`[2]
  grz.rnascope.res$ARTANOVA_Sex_p[i]   <-  rnasc.aov.res$`Pr(>F)`[1]
  grz.rnascope.res$ARTANOVA_Inter_p[i] <-  rnasc.aov.res$`Pr(>F)`[3]
}

grz.rnascope.res[sort(rownames(grz.rnascope.res)),]
#             ARTANOVA_Age_p ARTANOVA_Sex_p ARTANOVA_Inter_p 
# helz2_DAPI    1.530576e-01     0.40664383       0.03963966  --- n.s. on age (but sex/age interaction)
# ifih1_DAPI    2.703274e-03     0.11576495       0.32948629 
# irf3_DAPI     3.285764e-01     0.26276517       0.04986461  --- n.s. on age (but sex/age interaction)
# optn_DAPI     1.632979e-02     0.53527597       0.83736189 
# stat1a_DAPI   1.300118e-02     0.40284605       0.62798519 
# znfx1_DAPI    3.172747e-06     0.06988075       0.01973399 


#### b. ZMZ-1001
# update colnames
colnames(my.data.zmz)[-c(1:5)] <- gsub(":","_",colnames(my.data.zmz)[-c(1:5)])

# make factors to ART-ANOVA
my.data.zmz$Age <- factor(my.data.zmz$Age)
my.data.zmz$Sex <- factor(my.data.zmz$Sex)

# prepare output data frame
zmz.rnascope.res <- data.frame(matrix(0,6,3))
colnames(zmz.rnascope.res) <- c("ARTANOVA_Age_p","ARTANOVA_Sex_p","ARTANOVA_Inter_p")
rownames(zmz.rnascope.res) <- colnames(my.data.zmz)[-c(1:5)]

for (i in 1:6) {
  
  # Transform the data to ranks to compute ART ANOVA
  rnasc.tr.dat  <- art(as.formula(paste(colnames(my.data.zmz)[i +5], "~ Sex * Age")) , data=my.data.zmz)
  rnasc.aov.res <- anova(rnasc.tr.dat)
  # rnasc.aov.res
  #           Df Df.res F value   Pr(>F)  
  # 1 Sex      1     24 1.31309 0.263134  
  # 2 Age      2     24 4.21387 0.027011 *
  # 3 Sex:Age  2     24 0.16977 0.844860  
  
  # grab sex and age p-values
  zmz.rnascope.res$ARTANOVA_Age_p[i]   <-  rnasc.aov.res$`Pr(>F)`[2]
  zmz.rnascope.res$ARTANOVA_Sex_p[i]   <-  rnasc.aov.res$`Pr(>F)`[1]
  zmz.rnascope.res$ARTANOVA_Inter_p[i] <-  rnasc.aov.res$`Pr(>F)`[3]
}

zmz.rnascope.res[sort(rownames(zmz.rnascope.res)),]
#             ARTANOVA_Age_p ARTANOVA_Sex_p ARTANOVA_Inter_p
# helz2_DAPI    0.0023985783     0.35516230       0.70216892
# ifih1_DAPI    0.0114387794     0.58482933       0.53785757
# irf3_DAPI     0.0001033274     0.04003967       0.08838181
# optn_DAPI     0.0417001948     0.08598181       0.17224002
# stat1a_DAPI   0.0251434128     0.07158793       0.77462568
# znfx1_DAPI    0.0270105734     0.26313376       0.84486029

#############################################################################################


#######################
sink(file = paste0(Sys.Date(),"_DEGs_RNAscope_analysis_session_Info.txt"))
sessionInfo()
sink()


