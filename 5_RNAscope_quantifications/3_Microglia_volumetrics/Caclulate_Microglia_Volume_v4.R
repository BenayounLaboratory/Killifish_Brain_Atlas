setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/Miscroscopy/Microglia_Volume')

library(readxl)

library(outliers)
library(scales)

# 2025-03-13
# get microglia volume estimates from Rapheal's reconstructed surfaces
# using images from cell type analyses

# 2025-03-15
# without pre filtering of sizes

# 2025-08-08
# new file from Rapheal with the 5um size floor enforced
# includes the 2 GRZ and ZMZ cohorts for cell type quantifcation
# plot GRZ 13w as middle aged as they do not behave like 16w

####################################################################
## 0. read and preprocess data

# read data into tibbles
my.data.grz <- read_xlsx('MicrogliaVolume_Celltype_072825.xlsx', range = "B1:AC3953", sheet = "GRZ_No_Filter")
my.data.zmz <- read_xlsx('MicrogliaVolume_Celltype_072825.xlsx', range = "B1:AE2191", sheet = "ZMZ_No_Filter")

# make list to remove "NA'
my.data.grz.list <- as.list(my.data.grz)
my.data.grz.list <- lapply(my.data.grz.list,na.omit)

my.data.zmz.list <- as.list(my.data.zmz)
my.data.zmz.list <- lapply(my.data.zmz.list,na.omit)
####################################################################


####################################################################
## 1. Plot per image distribution data

pdf(paste0(Sys.Date(),"_GRZ_microglia_Volumes.pdf"))
boxplot(my.data.grz.list, log = 'y', 
        ylim = c(5,10000),
        las = 2,
        col = c(rep("deepskyblue" ,6),
                rep("deepskyblue4",8),
                rep("deeppink"    ,6),
                rep("deeppink4"   ,8)),
        main = "Apoeb surface volume u^3 (GRZ)")
dev.off()


pdf(paste0(Sys.Date(),"_ZMZ_microglia_Volumes.pdf"))
boxplot(my.data.zmz.list, log = 'y', 
        ylim = c(5,10000),
        las = 2,
        col = c(rep("deepskyblue",5),
                rep("deepskyblue4",5),
                rep("royalblue4",5),
                rep("deeppink",5),
                rep("deeppink4",5),
                rep("magenta4",5)),
        main = "Apoeb surface volume u^3 (ZMZ)")
dev.off()
####################################################################


####################################################################
## 2. Get sample medians

grz.medians <- lapply(my.data.grz.list, median)
zmz.medians <- lapply(my.data.zmz.list, median)

grz.median.samp <- list("GRZ_YF" = unlist(grz.medians[15:20]) ,
                        "GRZ_MF" = unlist(grz.medians[c(21,26)]) ,
                        "GRZ_OF" = unlist(grz.medians[c(22:25,27:28)]) ,
                        
                        "GRZ_YM" = unlist(grz.medians[1:6  ]) ,
                        "GRZ_MM" = unlist(grz.medians[c(8,11)  ]) ,
                        "GRZ_OM" = unlist(grz.medians[c(7,9:10,12:14) ])                         
)


zmz.median.samp <- list("ZMZ_YF" = unlist(zmz.medians[16:20])  ,
                        "ZMZ_OF" = unlist(zmz.medians[21:25])  ,
                        "ZMZ_GF" = unlist(zmz.medians[26:30])  ,
                        
                        "ZMZ_YM" = unlist(zmz.medians[1:5]  )  ,
                        "ZMZ_OM" = unlist(zmz.medians[6:10]  )  ,
                        "ZMZ_GM" = unlist(zmz.medians[11:15]  )  
) 


pdf(paste0(Sys.Date(),"_GRZ_microglia_Volumes_GroupMedians_NoBatchCorrection.pdf"))
boxplot(grz.median.samp,
        col = c("deeppink"     ,
                "deeppink3"    ,
                "deeppink4"    ,
                "deepskyblue"  ,
                "deepskyblue3" ,
                "deepskyblue4" ),
        las = 2,
        ylab = "Median Apoeb surface u^3 (A.U.)", 
        ylim = c(0,60))
beeswarm::beeswarm(grz.median.samp, add = T, pch = 16)
dev.off()  



pdf(paste0(Sys.Date(),"_ZMZ_microglia_Volumes_GroupMedians_NoBatchCorrection.pdf"))
boxplot(zmz.median.samp,
        col = c( "deeppink"     ,
                 "deeppink4"    ,
                 "magenta4"    ,
                 "deepskyblue"  ,
                 "deepskyblue4" ,
                 "royalblue4" 
        ),
        las = 2,
        ylab = "Median Apoeb surface u^3 (A.U.)", 
        ylim = c(0,60))
beeswarm::beeswarm(zmz.median.samp, add = T, pch = 16)
dev.off()  
####################################################################


####################################################################
############################     GRZ    ############################

# correct by processing batch
grz.b1 <- c(1, 3, 5, 7, 9, 11:13, 15, 17, 19, 21, 22, 24, 27, 28)
grz.b2 <- c(2, 4, 6, 8, 10, 14, 16, 18, 20, 23, 25, 26)

grz.medians.b1 <- grz.medians[grz.b1]
grz.medians.b2 <- grz.medians[grz.b2]

grz.median.samp.b1 <- list("GRZ_YF" = unlist(grz.medians.b1[9:11]) ,
                           "GRZ_MF" = unlist(grz.medians.b1[12]) ,
                           "GRZ_OF" = unlist(grz.medians.b1[13:16]) ,
                           
                           "GRZ_YM" = unlist(grz.medians.b1[1:3  ]) ,
                           "GRZ_MM" = unlist(grz.medians.b1[c(6) ]) ,
                           "GRZ_OM" = unlist(grz.medians.b1[c(4:5,7:8) ]) 
)

grz.median.samp.b2 <- list("GRZ_YF" = unlist(grz.medians.b2[7:9]) ,
                           "GRZ_MF" = unlist(grz.medians.b2[12]) ,
                           "GRZ_OF" = unlist(grz.medians.b2[10:11]) ,
                           
                           "GRZ_YM" = unlist(grz.medians.b2[1:3  ]) ,
                           "GRZ_MM" = unlist(grz.medians.b2[4]) ,
                           "GRZ_OM" = unlist(grz.medians.b2[5:6 ]) 
)



pdf(paste0(Sys.Date(),"_GRZ_microglia_Volumes_GroupMedians_SPLIT_Batches.pdf"))
par(mfrow = c(1,2))
boxplot(grz.median.samp.b1,
        col = c("deeppink"     ,
                "deeppink3"    ,
                "deeppink4"    ,
                "deepskyblue"  ,
                "deepskyblue3" ,
                "deepskyblue4" ),
        las = 2,
        ylab = "Median Apoeb surface u^3 (A.U.)", 
        ylim = c(0,50), 
        main = "batch 1")
beeswarm::beeswarm(grz.median.samp.b1, add = T, pch = 16)
boxplot(grz.median.samp.b2,
        col = c("deeppink"     ,
                "deeppink3"    ,
                "deeppink4"    ,
                "deepskyblue"  ,
                "deepskyblue3" ,
                "deepskyblue4" ),
        las = 2,
        ylab = "Median Apoeb surface u^3 (A.U.)", 
        ylim = c(0,50), 
        main = "batch 2")
beeswarm::beeswarm(grz.median.samp.b2, add = T, pch = 16)
par(mfrow = c(1,1))
dev.off()  


# normalize to batch
norm.grz.medians.b1 <- unlist(grz.medians.b1)/median(unlist(grz.medians.b1))
norm.grz.medians.b2 <- unlist(grz.medians.b2)/median(unlist(grz.medians.b2))


norm.grz.median.samp <- list("GRZ_YF" = c(unlist(norm.grz.medians.b1[9:11       ]), unlist(norm.grz.medians.b2[7:9  ]) ) ,
                             "GRZ_MF" = c(unlist(norm.grz.medians.b1[12         ]), unlist(norm.grz.medians.b2[12   ]) ) ,
                             "GRZ_OF" = c(unlist(norm.grz.medians.b1[13:16      ]), unlist(norm.grz.medians.b2[10:11]) ) ,
                             "GRZ_YM" = c(unlist(norm.grz.medians.b1[1:3        ]), unlist(norm.grz.medians.b2[1:3  ]) ) ,
                             "GRZ_MM" = c(unlist(norm.grz.medians.b1[c(6)       ]), unlist(norm.grz.medians.b2[4    ]) ) ,
                             "GRZ_OM" = c(unlist(norm.grz.medians.b1[c(4:5,7:8) ]), unlist(norm.grz.medians.b2[5:6  ]) ) )


pdf(paste0(Sys.Date(),"_GRZ_microglia_Volumes_GroupMedians_BATCH_NORM.pdf"))
boxplot(norm.grz.median.samp,
        col = c("deeppink"     ,
                "deeppink3"    ,
                "deeppink4"    ,
                "deepskyblue"  ,
                "deepskyblue3" ,
                "deepskyblue4" ),
        las = 2,
        ylab = "Median Apoeb surface (A.U.)", 
        ylim = c(0,2), outline = T)
beeswarm::beeswarm(norm.grz.median.samp, add = T, pch = 16)
dev.off()  

# looks like there is an outlier in YF and YM
lapply(norm.grz.median.samp, grubbs.test)

# $GRZ_YF
# Grubbs test for one outlier
# data:  X[[i]]
# G.F3_5wk = 1.958654, U = 0.079282, p-value = 0.007267                  ******
# alternative hypothesis: highest value 1.68019682072498 is an outlier
# 
# 
# $GRZ_MF
# Grubbs test for one outlier
# data:  X[[i]]
# G.F1_13wk = 0.70711, U = NA, p-value = NA
# alternative hypothesis: highest value 1.27880276443473 is an outlier
# 
# 
# $GRZ_OF
# Grubbs test for one outlier
# data:  X[[i]]
# G.F2_16wk = 1.6984, U = 0.3077, p-value = 0.1198
# alternative hypothesis: highest value 1.41843196450807 is an outlier
# 
# 
# $GRZ_YM
# Grubbs test for one outlier
# data:  X[[i]]
# G.M2_5wk = 1.84762, U = 0.18072, p-value = 0.03921                  ******
# alternative hypothesis: highest value 1.42503760705294 is an outlier
# 
# 
# $GRZ_MM
# Grubbs test for one outlier
# data:  X[[i]]
# G.M2_13wk = 0.70711, U = NA, p-value = NA
# alternative hypothesis: lowest value 0.986605292332819 is an outlier
# 
# 
# $GRZ_OM
# Grubbs test for one outlier
# data:  X[[i]]
# G.M2_16wk = 1.27392, U = 0.61051, p-value = 0.5562
# alternative hypothesis: highest value 1.45000231673653 is an outlier

# remove significant outliers
norm.grz.median.samp.no_out <- norm.grz.median.samp
norm.grz.median.samp.no_out$GRZ_YF <- norm.grz.median.samp.no_out$GRZ_YF[-3]
norm.grz.median.samp.no_out$GRZ_YM <- norm.grz.median.samp.no_out$GRZ_YM[-2]


# make tabular data for statistical testing
grz.tab.data <- data.frame("surface" = c(norm.grz.median.samp.no_out$GRZ_YF,
                                         norm.grz.median.samp.no_out$GRZ_MF,
                                         norm.grz.median.samp.no_out$GRZ_OF,
                                         norm.grz.median.samp.no_out$GRZ_YM,
                                         norm.grz.median.samp.no_out$GRZ_MM,
                                         norm.grz.median.samp.no_out$GRZ_OM ),
                           "Sex" = c(rep("Female", length(unlist(norm.grz.median.samp.no_out[1:3])) ),
                                     rep("Male"  , length(unlist(norm.grz.median.samp.no_out[4:6])) ) ),
                           "Age" = c(rep("Young" , length(norm.grz.median.samp.no_out[[1]])),
                                     rep("Middle", length(norm.grz.median.samp.no_out[[2]])),
                                     rep("Old"   , length(norm.grz.median.samp.no_out[[3]])),
                                     rep("Young" , length(norm.grz.median.samp.no_out[[4]])),
                                     rep("Middle", length(norm.grz.median.samp.no_out[[5]])),
                                     rep("Old"   , length(norm.grz.median.samp.no_out[[6]])) )  )


grz.aov <- summary(aov(surface ~ Sex + Age, data = grz.tab.data))
#             Df Sum Sq Mean Sq F value Pr(>F)  
# Sex          1 0.0001 0.00007   0.002 0.9625  
# Age          2 0.2899 0.14495   4.678 0.0203 *
# Residuals   22 0.6817 0.03099                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

shapiro.test(aov(surface ~ Sex + Age, data = my.tab.data)$residuals)
# Shapiro-Wilk normality test
# data:  aov(surface ~ Sex + Age, data = my.tab.data)$residuals
# W = 0.97554, p-value = 0.7681
###### Normality is not violated

pdf(paste0(Sys.Date(),"_GRZ_microglia_Volumes_GroupMedians_BATCH_NORM_outlier_removed.pdf"), width = 4, height = 5)
boxplot(norm.grz.median.samp.no_out,
        col = c("deeppink"     ,
                "deeppink3"    ,
                "deeppink4"    ,
                "deepskyblue"  ,
                "deepskyblue3" ,
                "deepskyblue4" ),
        las = 2,
        ylab = "Normalized median Apoeb surface (A.U.)", 
        ylim = c(0,2), outline = T)
beeswarm::beeswarm(norm.grz.median.samp.no_out, add = T, pch = 16)
text(0.5, 2  , paste("Age p = ", scientific(grz.aov[[1]]$`Pr(>F)`[2],2)), pos = 4)
text(0.5, 1.9, paste("Sex p = ", scientific(grz.aov[[1]]$`Pr(>F)`[1],2)), pos = 4)
dev.off()  


####################################################################################################



####################################################################
############################     ZMZ    ############################

# correct by processing batch
zmz.b1 <- c(1:3,  6:8, 11:13, 16:18, 21:23, 26:28)
zmz.b2 <- c(4:5, 9:10, 14:15, 19:20, 24:25, 29:30)

zmz.medians.b1 <- zmz.medians[zmz.b1]
zmz.medians.b2 <- zmz.medians[zmz.b2]

zmz.median.samp.b1 <- list("ZMZ_YF" = unlist(zmz.medians.b1[1:3  ]) ,
                           "ZMZ_OF" = unlist(zmz.medians.b1[4:6  ]) ,
                           "ZMZ_GF" = unlist(zmz.medians.b1[7:9  ]) ,
                           "ZMZ_YM" = unlist(zmz.medians.b1[10:12]) ,
                           "ZMZ_OM" = unlist(zmz.medians.b1[13:15]) ,
                           "ZMZ_GM" = unlist(zmz.medians.b1[16:18]) 
)

zmz.median.samp.b2 <- list("ZMZ_YF" = unlist(zmz.medians.b2[1:2  ]) ,
                           "ZMZ_OF" = unlist(zmz.medians.b2[3:4  ]) ,
                           "ZMZ_GF" = unlist(zmz.medians.b2[5:6  ]) ,
                           "ZMZ_YM" = unlist(zmz.medians.b2[7:8  ]) ,
                           "ZMZ_OM" = unlist(zmz.medians.b2[9:10 ]) ,
                           "ZMZ_GM" = unlist(zmz.medians.b2[11:12]) 
)



pdf(paste0(Sys.Date(),"_ZMZ_microglia_Volumes_GroupMedians_SPLIT_Batches.pdf"))
par(mfrow = c(1,2))
boxplot(zmz.median.samp.b1,
        col = c( "deeppink"     ,
                 "deeppink4"    ,
                 "magenta4"    ,
                 "deepskyblue"  ,
                 "deepskyblue4" ,
                 "royalblue4" 
        ),
        las = 2,
        ylab = "Median Apoeb surface u^3 (A.U.)", 
        ylim = c(0,60), 
        main = "batch 1")
beeswarm::beeswarm(zmz.median.samp.b1, add = T, pch = 16)
boxplot(zmz.median.samp.b2,
        col = c( "deeppink"     ,
                 "deeppink4"    ,
                 "magenta4"    ,
                 "deepskyblue"  ,
                 "deepskyblue4" ,
                 "royalblue4" 
        ),
        las = 2,
        ylab = "Median Apoeb surface u^3 (A.U.)", 
        ylim = c(0,60), 
        main = "batch 2")
beeswarm::beeswarm(zmz.median.samp.b2, add = T, pch = 16)
par(mfrow = c(1,1))
dev.off()  

# normalize to batch
norm.zmz.medians.b1 <- unlist(zmz.medians.b1)/median(unlist(zmz.medians.b1))
norm.zmz.medians.b2 <- unlist(zmz.medians.b2)/median(unlist(zmz.medians.b2))


norm.zmz.median.samp <- list("ZMZ_YF" = c(unlist(norm.zmz.medians.b1[10:12 ]), unlist(norm.zmz.medians.b2[7:8  ]) ) ,
                             "ZMZ_OF" = c(unlist(norm.zmz.medians.b1[13:15 ]), unlist(norm.zmz.medians.b2[9:10 ]) ) ,
                             "ZMZ_GF" = c(unlist(norm.zmz.medians.b1[16:18 ]), unlist(norm.zmz.medians.b2[11:12]) ) ,
                             "ZMZ_YM" = c(unlist(norm.zmz.medians.b1[1:3   ]), unlist(norm.zmz.medians.b2[1:2  ]) ) ,
                             "ZMZ_OM" = c(unlist(norm.zmz.medians.b1[4:6   ]), unlist(norm.zmz.medians.b2[3:4  ]) ) ,
                             "ZMZ_GM" = c(unlist(norm.zmz.medians.b1[7:9   ]), unlist(norm.zmz.medians.b2[5:6  ]) ) )


pdf(paste0(Sys.Date(),"_ZMZ_microglia_Volumes_GroupMedians_BATCH_NORM.pdf"))
boxplot(norm.zmz.median.samp,
        col = c( "deeppink"     ,
                 "deeppink4"    ,
                 "magenta4"    ,
                 "deepskyblue"  ,
                 "deepskyblue4" ,
                 "royalblue4" 
        ),
        las = 2,
        ylab = "Median Apoeb surface (A.U.)", 
        ylim = c(0,2), outline = T)
beeswarm::beeswarm(norm.zmz.median.samp, add = T, pch = 16)
dev.off()  

# looks like there is an outlier in OM and GM
lapply(norm.zmz.median.samp, grubbs.test)

# $ZMZ_YF
# Grubbs test for one outlier
# data:  X[[i]]
# G.F4_6w = 1.35720, U = 0.42438, p-value = 0.3426
# alternative hypothesis: lowest value 0.554290371549775 is an outlier
# 
# 
# $ZMZ_OF
# Grubbs test for one outlier
# data:  X[[i]]
# G.F1_16W = 1.41991, U = 0.36995, p-value = 0.2722
# alternative hypothesis: lowest value 0.762072009430692 is an outlier
# 
# 
# $ZMZ_GF
# Grubbs test for one outlier
# data:  X[[i]]
# G.F2_26W = 1.35573, U = 0.42562, p-value = 0.3443
# alternative hypothesis: highest value 1.29066204264515 is an outlier
# 
# 
# $ZMZ_YM
# Grubbs test for one outlier
# data:  X[[i]]
# G.M6_6w = 1.28478, U = 0.48417, p-value = 0.4294
# alternative hypothesis: highest value 1.23258094084207 is an outlier
# 
# 
# $ZMZ_OM
# Grubbs test for one outlier
# 
# data:  X[[i]]
# G.M3_16W = 1.51911, U = 0.27885, p-value = 0.1717
# alternative hypothesis: highest value 1.34419851941191 is an outlier
# 
# 
# $ZMZ_GM
# Grubbs test for one outlier
# data:  X[[i]]
# G.M4_26w = 1.68498, U = 0.11276, p-value = 0.04162                 ******
# alternative hypothesis: lowest value 0.872692794593014 is an outlier


# remove significant outliers
norm.zmz.median.samp.no_out <- norm.zmz.median.samp
norm.zmz.median.samp.no_out$ZMZ_GM <- norm.zmz.median.samp.no_out$ZMZ_GM[-4]


my.tab.data.zmz <- data.frame("surface" = c(norm.zmz.median.samp.no_out$ZMZ_YF,
                                            norm.zmz.median.samp.no_out$ZMZ_OF,
                                            norm.zmz.median.samp.no_out$ZMZ_GF,
                                            norm.zmz.median.samp.no_out$ZMZ_YM,
                                            norm.zmz.median.samp.no_out$ZMZ_OM,
                                            norm.zmz.median.samp.no_out$ZMZ_GM
                                            ),
                              "Sex" = c(rep("Female", length( unlist(norm.zmz.median.samp.no_out[1:3]) )),
                                        rep("Male"  , length( unlist(norm.zmz.median.samp.no_out[4:6]) ))   ),
                              "Age" = c(rep("Y", length(norm.zmz.median.samp.no_out[[1]])),
                                        rep("O", length(norm.zmz.median.samp.no_out[[2]])),
                                        rep("G", length(norm.zmz.median.samp.no_out[[3]])),
                                        rep("Y", length(norm.zmz.median.samp.no_out[[4]])),
                                        rep("O", length(norm.zmz.median.samp.no_out[[5]])),
                                        rep("G", length(norm.zmz.median.samp.no_out[[6]]))) )


zmz.aov <- summary(aov(surface ~ Sex + Age, data = my.tab.data.zmz))
#             Df Sum Sq Mean Sq F value  Pr(>F)   
# Sex          1 0.0502 0.05018   1.849 0.18609   
# Age          2 0.3779 0.18893   6.960 0.00395 **
# Residuals   25 0.6786 0.02714                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

shapiro.test(aov(surface ~ Sex + Age, data = my.tab.data.zmz)$residuals)
# Shapiro-Wilk normality test
# data:  aov(surface ~ Sex + Age, data = my.tab.data.zmz)$residuals
# W = 0.97728, p-value = 0.7654
###### Normality is not violated

pdf(paste0(Sys.Date(),"_ZMZ_microglia_Volumes_GroupMedians_BATCH_NORM_outlier_removed.pdf"), width = 4, height = 5)
boxplot(norm.zmz.median.samp.no_out,
        col = c( "deeppink"     ,
                 "deeppink4"    ,
                 "magenta4"    ,
                 "deepskyblue"  ,
                 "deepskyblue4" ,
                 "royalblue4" 
        ),
        las = 2,
        ylab = "Normalized median Apoeb surface (A.U.)", 
        ylim = c(0,2), outline = T)
beeswarm::beeswarm(norm.zmz.median.samp.no_out, add = T, pch = 16)
text(0.5, 2  , paste("Age p = ", scientific(zmz.aov[[1]]$`Pr(>F)`[2],2)), pos = 4)
text(0.5, 1.9, paste("Sex p = ", scientific(zmz.aov[[1]]$`Pr(>F)`[1],2)), pos = 4)
dev.off()  
#############################################################################################


#######################
sink(file = paste0(Sys.Date(),"_Microglia_Volume_analysis_session_Info.txt"))
sessionInfo()
sink()


