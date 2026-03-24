setwd('/Volumes/BB_Home_HQ2/Killifish_genome_variation/Genome/Comparison/')
options(stringsAsFactors = F)

# 2021-09-28
# add in additional samples from public data
# https://speciationgenomics.github.io/pca/

# load tidyverse package
# library(tidyverse)
library(dbplyr)
library(dplyr)

# read in data
pca      <- readr::read_table2("2021-09-28_Killifish_Strains_genotypes_PRUNED_PLINK_pca.eigenvec", col_names = FALSE)
eigenval <- scan("2021-09-28_Killifish_Strains_genotypes_PRUNED_PLINK_pca.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]

# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:18, pve = eigenval/sum(eigenval)*100)

# make plot
pdf(paste0(Sys.Date(),"_Variance_explained.pdf"))
a <- ggplot2::ggplot(pve, ggplot2::aes(PC, pve)) + ggplot2::geom_bar(stat = "identity")
a + ggplot2::ylab("Percentage variance explained") + ggplot2::theme_light()
dev.off()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)
#  [1]  30.24823  42.55538  53.00442  61.51588  69.25695  76.08351  82.50285  86.80572  89.51416  91.42524  93.08454  94.60446  95.87329
#  [14]  97.01090  98.11381  99.13036 100.02554 100.00000

# remake data.frame
pca <- data.frame(pca)
pca$ind2 <- c("GRZ1",
              "GRZ2",
              "ZMZ1",
              "ZMZ2",
              "MZM0703_stanford",
              "MZM0403_stanford",
              "GRZ_P0FA",
              "GRZ_BGI170",
              "GRZ_GRZ340",
              "GRZ_GRZ300",
              "MZZW0701_M002",
              "MZZW0701_M001",
              "MZM0410_M013",
              "MZM0410_M012",
              "MZM0403_M007",
              "MZM0403_M006",
              "GRZ_M005",
              "GRZ_M004")

# plot pca
my.cols <- c("yellow",
             "yellow",
             "gold",
             "gold",
             "red",
             "red",
             "yellow",
             "yellow",
             "yellow",
             "yellow",
             "gold",
             "gold",
             "red",
             "red",
             "red",
             "red",
             "yellow",
             "yellow")
my.cols.tr <- grDevices::adjustcolor(my.cols, alpha.f = 0.5)

my.pch <- c(16,
            16,
            16,
            16,
            15,
            15,
            15,
            15,
            15,
            15,
            17,
            17,
            17,
            17,
            17,
            17,
            17,
            17)

my.pch2 <- c(1,
             1,
             1,
             1,
             1,
             0,
             0,
             0,
             0,
             0,
             2,
             2,
             2,
             2,
             2,
             2,
             2,
             2)

pdf(paste0(Sys.Date(),"Killifish_Strain_Genotypes_PCA_plot_text.pdf"))
plot(pca$PC1,pca$PC2, 
     pch = my.pch, cex = 3, col = my.cols.tr,
     xlab = paste('PC1 (', round(pve$pve[1],1),"%)", sep=""),
     ylab = paste('PC2 (', round(pve$pve[2],1),"%)", sep=""),
     cex.lab = 1.5,
     cex.axis = 1.5,
     main = "Genotype PCA", 
     xlim = c(-0.2,0.6),
     ylim = c(-0.5,0.6) ) 
# points(pca$PC1,pca$PC2, 
#      pch = my.pch2, cex = 2, col = "grey")
text(pca$PC1,pca$PC2, pca$ind2, cex = 0.7, col = "grey")
dev.off()



pdf(paste0(Sys.Date(),"Killifish_Strain_Genotypes_PCA_plot_pict.pdf"))
plot(pca$PC1,pca$PC2, 
     pch = my.pch, cex = 3, col = my.cols.tr,
     xlab = paste('PC1 (', round(pve$pve[1],1),"%)", sep=""),
     ylab = paste('PC2 (', round(pve$pve[2],1),"%)", sep=""),
     cex.lab = 1.5,
     cex.axis = 1.5,
     main = "Genotype PCA", 
     xlim = c(-0.2,0.6),
     ylim = c(-0.5,0.6) ) 
points(pca$PC1,pca$PC2, 
       pch = my.pch2, cex = 3, col = "grey")
# text(pca$PC1,pca$PC2, pca$ind2, cex = 0.5, col = "grey")
legend("bottomleft",c("GRZ","Yellow_LL (ZMZ1001,MZZW0701)","Red_LL (MZM0403,MZM0703,MZM0410)"),col = c("yellow","gold","red"), pch = 16, cex = 1, pt.cex = 2, bty = "n")
legend("topleft",c("USC","Stanford","FLI"),col = "grey", pch = c(1,0,2), cex = 1, pt.cex = 2, bty = "n")
dev.off()
