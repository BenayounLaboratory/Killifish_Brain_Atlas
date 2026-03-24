setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/ATAC/Bulk_ATAC_Analysis_Folder/Motif_Analysis')
options(stringsAsFactors = FALSE)


# 2024-04-09
# Killifish brain bulk ATAC aging
# Plot nr3c1 site density around DE peaks


################# NR3C1 #################
grz.nr3c1.dat.up   <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_nr3c1_sites.txt', sep = "\t")
grz.nr3c1.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_nr3c1_sites.txt', sep = "\t")
grz.nr3c1.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_nr3c1_sites.txt', sep = "\t")

zmz.nr3c1.dat.up   <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_nr3c1_sites.txt', sep = "\t")
zmz.nr3c1.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_nr3c1_sites.txt', sep = "\t")
zmz.nr3c1.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_nr3c1_sites.txt', sep = "\t")



pdf(paste0(Sys.Date(),"_Nr3c1_site_density_DA_Peaks_GRZ.pdf"), height = 5, width = 6)
plot(grz.nr3c1.dat.up[,1], grz.nr3c1.dat.up[,2], type = 'l', 
     ylim = c(0,5e-4), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "Nr3c1 Sites per bp per peak",
     main = "Killi Brain ATAC (GRZ)", las = 1)
points(grz.nr3c1.dat.dwn[,1], grz.nr3c1.dat.dwn[,2], type = 'l', col = "#333399")
points(grz.nr3c1.dat.bcgd[,1], grz.nr3c1.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()


pdf(paste0(Sys.Date(),"_Nr3c1_site_density_DA_Peaks_ZMZ.pdf"), height = 5, width = 6)
plot(zmz.nr3c1.dat.up[,1], zmz.nr3c1.dat.up[,2], type = 'l', 
     ylim = c(0,5e-4), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "Nr3c1 Sites per bp per peak",
     main = "Killi Brain ATAC (ZMZ)", las = 1)
points(zmz.nr3c1.dat.dwn[,1], zmz.nr3c1.dat.dwn[,2], type = 'l', col = "#333399")
points(zmz.nr3c1.dat.bcgd[,1], zmz.nr3c1.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()



################# STAT1 #################
grz.stat1.dat.up   <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_stat1_sites.txt', sep = "\t")
grz.stat1.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_stat1_sites.txt', sep = "\t")
grz.stat1.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_stat1_sites.txt', sep = "\t")

zmz.stat1.dat.up   <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_stat1_sites.txt', sep = "\t")
zmz.stat1.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_stat1_sites.txt', sep = "\t")
zmz.stat1.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_stat1_sites.txt', sep = "\t")



pdf(paste0(Sys.Date(),"_Stat1_site_density_DA_Peaks_GRZ.pdf"), height = 5, width = 6)
plot(grz.stat1.dat.up[,1], grz.stat1.dat.up[,2], type = 'l', 
     ylim = c(0,5e-4), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "Stat1 Sites per bp per peak",
     main = "Killi Brain ATAC (GRZ)", las = 1)
points(grz.stat1.dat.dwn[,1], grz.stat1.dat.dwn[,2], type = 'l', col = "#333399")
points(grz.stat1.dat.bcgd[,1], grz.stat1.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()


pdf(paste0(Sys.Date(),"_Stat1_site_density_DA_Peaks_ZMZ.pdf"), height = 5, width = 6)
plot(zmz.stat1.dat.up[,1], zmz.stat1.dat.up[,2], type = 'l', 
     ylim = c(0,5e-4), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "Stat1 Sites per bp per peak",
     main = "Killi Brain ATAC (ZMZ)", las = 1)
points(zmz.stat1.dat.dwn[,1], zmz.stat1.dat.dwn[,2], type = 'l', col = "#333399")
points(zmz.stat1.dat.bcgd[,1], zmz.stat1.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()



################# EGR1 #################
grz.egr1.dat.up   <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_egr1_sites.txt', sep = "\t")
grz.egr1.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_egr1_sites.txt', sep = "\t")
grz.egr1.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_egr1_sites.txt', sep = "\t")

zmz.egr1.dat.up   <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_egr1_sites.txt', sep = "\t")
zmz.egr1.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_egr1_sites.txt', sep = "\t")
zmz.egr1.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_egr1_sites.txt', sep = "\t")


pdf(paste0(Sys.Date(),"_egr1_site_density_DA_Peaks_GRZ.pdf"), height = 5, width = 6)
plot(grz.egr1.dat.up[,1], grz.egr1.dat.up[,2], type = 'l', 
     ylim = c(0,1.5e-3), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "egr1 Sites per bp per peak",
     main = "Killi Brain ATAC (GRZ)", las = 1)
points(grz.egr1.dat.dwn[,1], grz.egr1.dat.dwn[,2], type = 'l', col = "#333399")
points(grz.egr1.dat.bcgd[,1], grz.egr1.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()


pdf(paste0(Sys.Date(),"_egr1_site_density_DA_Peaks_ZMZ.pdf"), height = 5, width = 6)
plot(zmz.egr1.dat.up[,1], zmz.egr1.dat.up[,2], type = 'l', 
     ylim = c(0,1.5e-3), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "egr1 Sites per bp per peak",
     main = "Killi Brain ATAC (ZMZ)", las = 1)
points(zmz.egr1.dat.dwn[,1], zmz.egr1.dat.dwn[,2], type = 'l', col = "#333399")
points(zmz.egr1.dat.bcgd[,1], zmz.egr1.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()


################# IRF3 #################
grz.irf3.dat.up   <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_irf3_sites.txt', sep = "\t")
grz.irf3.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_irf3_sites.txt', sep = "\t")
grz.irf3.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_irf3_sites.txt', sep = "\t")

zmz.irf3.dat.up   <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_irf3_sites.txt', sep = "\t")
zmz.irf3.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_irf3_sites.txt', sep = "\t")
zmz.irf3.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_irf3_sites.txt', sep = "\t")


pdf(paste0(Sys.Date(),"_irf3_site_density_DA_Peaks_GRZ.pdf"), height = 5, width = 6)
plot(grz.irf3.dat.up[,1], grz.irf3.dat.up[,2], type = 'l', 
     ylim = c(0,5e-4), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "irf3 Sites per bp per peak",
     main = "Killi Brain ATAC (GRZ)", las = 1)
points(grz.irf3.dat.dwn[,1], grz.irf3.dat.dwn[,2], type = 'l', col = "#333399")
points(grz.irf3.dat.bcgd[,1], grz.irf3.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()


pdf(paste0(Sys.Date(),"_irf3_site_density_DA_Peaks_ZMZ.pdf"), height = 5, width = 6)
plot(zmz.irf3.dat.up[,1], zmz.irf3.dat.up[,2], type = 'l', 
     ylim = c(0,5e-4), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "irf3 Sites per bp per peak",
     main = "Killi Brain ATAC (ZMZ)", las = 1)
points(zmz.irf3.dat.dwn[,1], zmz.irf3.dat.dwn[,2], type = 'l', col = "#333399")
points(zmz.irf3.dat.bcgd[,1], zmz.irf3.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()


################# Pax6 #################
grz.pax6.dat.up   <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_pax6_sites.txt', sep = "\t")
grz.pax6.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_pax6_sites.txt', sep = "\t")
grz.pax6.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_pax6_sites.txt', sep = "\t")

zmz.pax6.dat.up   <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_pax6_sites.txt', sep = "\t")
zmz.pax6.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_pax6_sites.txt', sep = "\t")
zmz.pax6.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_pax6_sites.txt', sep = "\t")


pdf(paste0(Sys.Date(),"_pax6_site_density_DA_Peaks_GRZ.pdf"), height = 5, width = 6)
plot(grz.pax6.dat.up[,1], grz.pax6.dat.up[,2], type = 'l', 
     ylim = c(0,1e-4), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "pax6 Sites per bp per peak",
     main = "Killi Brain ATAC (GRZ)", las = 1)
points(grz.pax6.dat.dwn[,1], grz.pax6.dat.dwn[,2], type = 'l', col = "#333399")
points(grz.pax6.dat.bcgd[,1], grz.pax6.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()


pdf(paste0(Sys.Date(),"_pax6_site_density_DA_Peaks_ZMZ.pdf"), height = 5, width = 6)
plot(zmz.pax6.dat.up[,1], zmz.pax6.dat.up[,2], type = 'l', 
     ylim = c(0,1e-4), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "pax6 Sites per bp per peak",
     main = "Killi Brain ATAC (ZMZ)", las = 1)
points(zmz.pax6.dat.dwn[,1], zmz.pax6.dat.dwn[,2], type = 'l', col = "#333399")
points(zmz.pax6.dat.bcgd[,1], zmz.pax6.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()


################# RELA #################
grz.rela.dat.up   <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_rela_sites.txt', sep = "\t")
grz.rela.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_rela_sites.txt', sep = "\t")
grz.rela.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_rela_sites.txt', sep = "\t")

zmz.rela.dat.up   <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_rela_sites.txt', sep = "\t")
zmz.rela.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_rela_sites.txt', sep = "\t")
zmz.rela.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_rela_sites.txt', sep = "\t")


pdf(paste0(Sys.Date(),"_rela_site_density_DA_Peaks_GRZ.pdf"), height = 5, width = 6)
plot(grz.rela.dat.up[,1], grz.rela.dat.up[,2], type = 'l',
     ylim = c(0,5e-4), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "rela Sites per bp per peak",
     main = "Killi Brain ATAC (GRZ)", las = 1)
points(grz.rela.dat.dwn[,1], grz.rela.dat.dwn[,2], type = 'l', col = "#333399")
points(grz.rela.dat.bcgd[,1], grz.rela.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()


pdf(paste0(Sys.Date(),"_rela_site_density_DA_Peaks_ZMZ.pdf"), height = 5, width = 6)
plot(zmz.rela.dat.up[,1], zmz.rela.dat.up[,2], type = 'l',
     ylim = c(0,5e-4), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "rela Sites per bp per peak",
     main = "Killi Brain ATAC (ZMZ)", las = 1)
points(zmz.rela.dat.dwn[,1], zmz.rela.dat.dwn[,2], type = 'l', col = "#333399")
points(zmz.rela.dat.bcgd[,1], zmz.rela.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()


################# NFYA #################
grz.nfya.dat.up   <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_nfy_sites.txt', sep = "\t")
grz.nfya.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_nfy_sites.txt', sep = "\t")
grz.nfya.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_nfy_sites.txt', sep = "\t")

zmz.nfya.dat.up   <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_nfy_sites.txt', sep = "\t")
zmz.nfya.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_nfy_sites.txt', sep = "\t")
zmz.nfya.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_nfy_sites.txt', sep = "\t")


pdf(paste0(Sys.Date(),"_nfya_site_density_DA_Peaks_GRZ.pdf"), height = 5, width = 6)
plot(grz.nfya.dat.up[,1], grz.nfya.dat.up[,2], type = 'l',
     ylim = c(0,1.5e-3), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "nfya Sites per bp per peak",
     main = "Killi Brain ATAC (GRZ)", las = 1)
points(grz.nfya.dat.dwn[,1], grz.nfya.dat.dwn[,2], type = 'l', col = "#333399")
points(grz.nfya.dat.bcgd[,1], grz.nfya.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()


pdf(paste0(Sys.Date(),"_nfya_site_density_DA_Peaks_ZMZ.pdf"), height = 5, width = 6)
plot(zmz.nfya.dat.up[,1], zmz.nfya.dat.up[,2], type = 'l',
     ylim = c(0,1.5e-3), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "nfya Sites per bp per peak",
     main = "Killi Brain ATAC (ZMZ)", las = 1)
points(zmz.nfya.dat.dwn[,1], zmz.nfya.dat.dwn[,2], type = 'l', col = "#333399")
points(zmz.nfya.dat.bcgd[,1], zmz.nfya.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()


################# NFKB1 #################
grz.nfkb1.dat.up   <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_nfkb1_sites.txt', sep = "\t")
grz.nfkb1.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_nfkb1_sites.txt', sep = "\t")
grz.nfkb1.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_nfkb1_sites.txt', sep = "\t")

zmz.nfkb1.dat.up   <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_nfkb1_sites.txt', sep = "\t")
zmz.nfkb1.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_nfkb1_sites.txt', sep = "\t")
zmz.nfkb1.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_nfkb1_sites.txt', sep = "\t")


pdf(paste0(Sys.Date(),"_nfkb1_site_density_DA_Peaks_GRZ.pdf"), height = 5, width = 6)
plot(grz.nfkb1.dat.up[,1], grz.nfkb1.dat.up[,2], type = 'l',
     ylim = c(0,1e-4), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "nfkb1 Sites per bp per peak",
     main = "Killi Brain ATAC (GRZ)", las = 1)
points(grz.nfkb1.dat.dwn[,1], grz.nfkb1.dat.dwn[,2], type = 'l', col = "#333399")
points(grz.nfkb1.dat.bcgd[,1], grz.nfkb1.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()


pdf(paste0(Sys.Date(),"_nfkb1_site_density_DA_Peaks_ZMZ.pdf"), height = 5, width = 6)
plot(zmz.nfkb1.dat.up[,1], zmz.nfkb1.dat.up[,2], type = 'l',
     ylim = c(0,1e-4), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "nfkb1 Sites per bp per peak",
     main = "Killi Brain ATAC (ZMZ)", las = 1)
points(zmz.nfkb1.dat.dwn[,1], zmz.nfkb1.dat.dwn[,2], type = 'l', col = "#333399")
points(zmz.nfkb1.dat.bcgd[,1], zmz.nfkb1.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()


################# SOX2 #################
grz.sox2.dat.up   <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_sox2_sites.txt', sep = "\t")
grz.sox2.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_sox2_sites.txt', sep = "\t")
grz.sox2.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_sox2_sites.txt', sep = "\t")

zmz.sox2.dat.up   <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_sox2_sites.txt', sep = "\t")
zmz.sox2.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_sox2_sites.txt', sep = "\t")
zmz.sox2.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_sox2_sites.txt', sep = "\t")


pdf(paste0(Sys.Date(),"_sox2_site_density_DA_Peaks_GRZ.pdf"), height = 5, width = 6)
plot(grz.sox2.dat.up[,1], grz.sox2.dat.up[,2], type = 'l',
     ylim = c(0,1.5e-3), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "sox2 Sites per bp per peak",
     main = "Killi Brain ATAC (GRZ)", las = 1)
points(grz.sox2.dat.dwn[,1], grz.sox2.dat.dwn[,2], type = 'l', col = "#333399")
points(grz.sox2.dat.bcgd[,1], grz.sox2.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()


pdf(paste0(Sys.Date(),"_sox2_site_density_DA_Peaks_ZMZ.pdf"), height = 5, width = 6)
plot(zmz.sox2.dat.up[,1], zmz.sox2.dat.up[,2], type = 'l',
     ylim = c(0,1.5e-3), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "sox2 Sites per bp per peak",
     main = "Killi Brain ATAC (ZMZ)", las = 1)
points(zmz.sox2.dat.dwn[,1], zmz.sox2.dat.dwn[,2], type = 'l', col = "#333399")
points(zmz.sox2.dat.bcgd[,1], zmz.sox2.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()

################# BHLHE40 #################
grz.bhlhe40.dat.up   <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_bhlhe40_sites.txt', sep = "\t")
grz.bhlhe40.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_bhlhe40_sites.txt', sep = "\t")
grz.bhlhe40.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_bhlhe40_sites.txt', sep = "\t")

zmz.bhlhe40.dat.up   <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age_HOMER_hist_bhlhe40_sites.txt', sep = "\t")
zmz.bhlhe40.dat.dwn  <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age_HOMER_hist_bhlhe40_sites.txt', sep = "\t")
zmz.bhlhe40.dat.bcgd <- read.csv('./TF_site_histo/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_Background_HOMER_hist_bhlhe40_sites.txt', sep = "\t")


pdf(paste0(Sys.Date(),"_bhlhe40_site_density_DA_Peaks_GRZ.pdf"), height = 5, width = 6)
plot(grz.bhlhe40.dat.up[,1], grz.bhlhe40.dat.up[,2], type = 'l',
     ylim = c(0,1.5e-3), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "bhlhe40 Sites per bp per peak",
     main = "Killi Brain ATAC (GRZ)", las = 1)
points(grz.bhlhe40.dat.dwn[,1], grz.bhlhe40.dat.dwn[,2], type = 'l', col = "#333399")
points(grz.bhlhe40.dat.bcgd[,1], grz.bhlhe40.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()


pdf(paste0(Sys.Date(),"_bhlhe40_site_density_DA_Peaks_ZMZ.pdf"), height = 5, width = 6)
plot(zmz.bhlhe40.dat.up[,1], zmz.bhlhe40.dat.up[,2], type = 'l',
     ylim = c(0,1.5e-3), col = "#CC3333",
     xlab = "Distance from Peak Center (bp)",
     ylab = "bhlhe40 Sites per bp per peak",
     main = "Killi Brain ATAC (ZMZ)", las = 1)
points(zmz.bhlhe40.dat.dwn[,1], zmz.bhlhe40.dat.dwn[,2], type = 'l', col = "#333399")
points(zmz.bhlhe40.dat.bcgd[,1], zmz.bhlhe40.dat.bcgd[,2], type = 'l', col = "grey")
legend("topleft", c("Increased Accessibility", "Decreased Accessibility","Background"), pch = "_", col = c("#CC3333","#333399","grey"), bty = 'n')
dev.off()

#######################
sink(file = paste(Sys.Date(),"_Plot_Homer_motif_Densities_session_Info.txt", sep =""))
sessionInfo()
sink()


