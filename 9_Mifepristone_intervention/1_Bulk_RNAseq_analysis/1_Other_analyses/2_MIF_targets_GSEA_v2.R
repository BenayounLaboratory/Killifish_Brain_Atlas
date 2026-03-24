setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/GR_signaling/Mifepristone/bulk_Brain_RNAseq/DESeq2')
options(stringsAsFactors = F)

library("DESeq2")        #
library("sva")           #
library("limma")         #
library("pheatmap")      #
library("bitops")        #
library(phenoTest)

library(ggplot2) 
library(scales) 
theme_set(theme_bw())

library(Vennerable)

# 2025-07-15
# Compare aging and mif regulated genes

# 2025-07-28
# rerun with appropriate background size for Fisher

#####################################################################################################################
#### 1. Run GSEA

load('2025-07-15_Brain_MIF_DEseq2_results.RData')
ls()
# [1] "F.mif.res" "M.mif.res"


###############################################################################################################################################
# Geneset enrichment
F.mif.res$Mif           <- data.frame(F.mif.res$Mif)
F.mif.res$GeneName      <- rownames(F.mif.res$Mif )

M.mif.res$Mif           <- data.frame(M.mif.res$Mif)
M.mif.res$GeneName      <- rownames(M.mif.res$Mif )

######################## B. Prepare GeneLists using DEseq2 t-statistic to rank genes ########################
F.res.mif.geneList         <- F.mif.res$Mif$stat
names(F.res.mif.geneList)  <- rownames(F.mif.res$Mif )
F.res.mif.geneList         <- sort(F.res.mif.geneList , decreasing = TRUE)

M.res.mif.geneList         <- M.mif.res$Mif$stat
names(M.res.mif.geneList)  <- rownames(M.mif.res$Mif )
M.res.mif.geneList         <- sort(M.res.mif.geneList , decreasing = TRUE)

######################## C. Prep gene Sets ########################

my.age.gs.F <- list("Age_up"       = rownames(F.mif.res$Aging)[bitAnd(F.mif.res$Aging$padj < 0.05,F.mif.res$Aging$log2FoldChange > 0)>0 ],
                    "Age_dwn"      = rownames(F.mif.res$Aging)[bitAnd(F.mif.res$Aging$padj < 0.05,F.mif.res$Aging$log2FoldChange < 0)>0 ])


my.age.gs.M <- list("Age_up"       = rownames(M.mif.res$Aging)[bitAnd(M.mif.res$Aging$padj < 0.05,M.mif.res$Aging$log2FoldChange > 0)>0 ],
                    "Age_dwn"      = rownames(M.mif.res$Aging)[bitAnd(M.mif.res$Aging$padj < 0.05,M.mif.res$Aging$log2FoldChange < 0)>0 ])

######################## C. Gene Set Enrichment Analysis ########################
# set seed to stabilize output
set.seed(1234567890)

# run phenotest GSEA
gsea.res.F <- gsea( x         =  F.res.mif.geneList ,
                    gsets     =  my.age.gs.F,
                    mc.cores  =  1                 ,
                    logScale  =  FALSE             ,
                    B         =  10000              ,
                    minGenes  =  5                 ,
                    maxGenes  =  5000               )
my.summary.F <- data.frame(summary(gsea.res.F))
gsea.res.F$significance$summary
#            n         es       nes pval.es     pval.nes          fdr
# Age_up  1471 -0.5524104 -2.994689       0 0.000000e+00 0.000000e+00
# Age_dwn 1624  0.4826364  2.766773       0 2.220446e-16 2.220446e-16

pdf(paste(Sys.Date(), "Mif_responsive_GSEA_plot_F_aging.pdf", sep = "_"))
plot.gseaData(gsea.res.F, es.nes='nes', selGsets='Age_up', color = "purple")
plot.gseaData(gsea.res.F, es.nes='nes', selGsets='Age_dwn', color = "purple")
dev.off()



# run phenotest GSEA
gsea.res.M <- gsea( x         =  M.res.mif.geneList ,
                    gsets     =  my.age.gs.M,
                    mc.cores  =  1                 ,
                    logScale  =  FALSE             ,
                    B         =  10000              ,
                    minGenes  =  5                 ,
                    maxGenes  =  5000               )
my.summary.M <- data.frame(summary(gsea.res.M))
gsea.res.M$significance$summary
#            n         es       nes      pval.es     pval.nes          fdr
# Age_up  3002 -0.5832228 -3.550660 8.881784e-16 2.220446e-16 2.220446e-16
# Age_dwn 3357  0.5637614  3.365225 0.000000e+00 2.220446e-16 2.220446e-16

pdf(paste(Sys.Date(), "Mif_responsive_GSEA_plot_M_aging.pdf", sep = "_"))
plot.gseaData(gsea.res.M, es.nes='nes', selGsets='Age_up', color = "purple")
plot.gseaData(gsea.res.M, es.nes='nes', selGsets='Age_dwn', color = "purple")
dev.off()



############################################################################################
# Make bubble chart summary

############## Plot 
F.gsea.summary <- data.frame(gsea.res.F$significance$summary)
M.gsea.summary <- data.frame(gsea.res.M$significance$summary)

F.gsea.summary <- cbind(rownames(F.gsea.summary),F.gsea.summary)
M.gsea.summary <- cbind(rownames(M.gsea.summary),M.gsea.summary)

colnames(F.gsea.summary)[1] <- "GeneSet"
colnames(M.gsea.summary)[1] <- "GeneSet"

my.mif <- rbind(F.gsea.summary,M.gsea.summary)


# get merged datafame for ggplot
my.mif$minusLog10FDR <- -log10(my.mif$fdr + 1e-30)
my.mif$condition <- c(rep("Mif_F",2), rep("Mif_M",2))

my.max <- 4
my.min <- -4
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")

# to preserve the wanted order
my.mif$condition <- factor(my.mif$condition, levels = unique(my.mif$condition))
my.mif$GeneSet   <- factor(my.mif$GeneSet, levels = rev(unique(my.mif$GeneSet)))

pdf(paste0(Sys.Date(),"_F_Aging_GSEA_mif_GeneSets.pdf"),height = 3.5, width=6)
my.plot <- ggplot(my.mif,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme_bw()+ geom_point(shape = 16)
my.plot <- ggplot(my.mif,aes(x=condition,y=GeneSet,colour=nes,size=minusLog10FDR))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("GSEA") + labs(x = "Regulation in response to Mif", y = "Aging GeneSet")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(my.min,my.max))
my.plot <- my.plot + scale_size_area(limits = c(10,30))
print(my.plot)
dev.off()
#####################################################################################################################


#####################################################################################################################
## 2. Venn diagram


####

my.age.gs.F <- list("Age_up"       = rownames(F.mif.res$Aging)[bitAnd(F.mif.res$Aging$padj < 0.05,F.mif.res$Aging$log2FoldChange > 0)>0 ],
                    "Age_dwn"      = rownames(F.mif.res$Aging)[bitAnd(F.mif.res$Aging$padj < 0.05,F.mif.res$Aging$log2FoldChange < 0)>0 ])


my.age.gs.M <- list("Age_up"       = rownames(M.mif.res$Aging)[bitAnd(M.mif.res$Aging$padj < 0.05,M.mif.res$Aging$log2FoldChange > 0)>0 ],
                    "Age_dwn"      = rownames(M.mif.res$Aging)[bitAnd(M.mif.res$Aging$padj < 0.05,M.mif.res$Aging$log2FoldChange < 0)>0 ])


my.mif.gs.F <- list("Mif_up"       = rownames(F.mif.res$Mif)[bitAnd(F.mif.res$Mif$padj < 0.05,F.mif.res$Mif$log2FoldChange > 0)>0 ],
                    "Mif_dwn"      = rownames(F.mif.res$Mif)[bitAnd(F.mif.res$Mif$padj < 0.05,F.mif.res$Mif$log2FoldChange < 0)>0 ])


my.mif.gs.M <- list("Mif_up"       = rownames(M.mif.res$Mif)[bitAnd(M.mif.res$Mif$padj < 0.05,M.mif.res$Mif$log2FoldChange > 0)>0 ],
                    "Mif_dwn"      = rownames(M.mif.res$Mif)[bitAnd(M.mif.res$Mif$padj < 0.05,M.mif.res$Mif$log2FoldChange < 0)>0 ])




###########################
##### Female/C1 analysis

####### a. Female: Aging up/Mif down
F.AgeU.MifD <- list("Female Aging Up"   = my.age.gs.F$Age_up,
                    "Female Mif Down"   = my.mif.gs.F$Mif_dwn)
my.Venn <- Venn(F.AgeU.MifD)

pdf(paste0(Sys.Date(),"_Female_Aging_Up_MIf_Down.pdf"))
plot(my.Venn, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
dev.off()

test.F.AgeU.MifD <- fisher.test(matrix(c(96, 1375, 173, length(union(rownames(F.mif.res$Aging),rownames(F.mif.res$Mif)))-96 -1375 - 173),2,2))
test.F.AgeU.MifD$p.value
# [1] 2.230605e-48

F.AUMD_g <- intersect(my.age.gs.F$Age_up,my.mif.gs.F$Mif_dwn)


####### b. Female: Aging down/Mif up
F.AgeD.MifU <- list("Female Aging Down" = my.age.gs.F$Age_dwn,
                    "Female Mif Up"   = my.mif.gs.F$Mif_up)
my.Venn <- Venn(F.AgeD.MifU)

pdf(paste0(Sys.Date(),"_Female_Aging_Down_MIf_Up.pdf"))
plot(my.Venn, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
dev.off()

test.F.AgeD.MifU <- fisher.test(matrix(c(87, 1537, 116, length(union(rownames(F.mif.res$Aging),rownames(F.mif.res$Mif))) -87-1537-116),2,2))
test.F.AgeD.MifU$p.value
# [1] 4.016391e-48


F.ADMU_g <- intersect(my.age.gs.F$Age_dwn,my.mif.gs.F$Mif_up)


###########################
##### Male/C2 analysis

####### a. Male: Aging up/Mif down
M.AgeU.MifD <- list("Male Aging Up" = my.age.gs.M$Age_up,
                    "Male Mif Down" = my.mif.gs.M$Mif_dwn)
my.Venn <- Venn(M.AgeU.MifD)

pdf(paste0(Sys.Date(),"_Male_Aging_Up_MIf_Down.pdf"))
plot(my.Venn, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
dev.off()

test.M.AgeU.MifD <- fisher.test(matrix(c(54, 2948, 17, length(union(rownames(M.mif.res$Aging),rownames(M.mif.res$Mif))) -54-2948-17),2,2))
test.M.AgeU.MifD$p.value
# [1]  8.516741e-35

M.AUMD_g <- intersect(my.age.gs.M$Age_up,my.mif.gs.M$Mif_dwn)
# [1] "cd74"           "lemd3"          "LOC107385152"   "LOC107389042"   "myof"           "dync2li1"       "abcg8"          "LOC107376218"   "ints9"          "LOC107378257"  
# [11] "cnot10"         "invs"           "pbxip1"         "LOC107379864"   "rpl35a"         "LOC107380389"   "csgr07h15orf43" "csgr08h6orf120" "rrp36"          "LOC107382316"  
# [21] "palmd"          "myeov2"         "LOC107383058"   "LOC107383107"   "LOC107383647"   "LOC107384504"   "apbb3"          "crocc2"         "espnl"          "scly"          
# [31] "ush1c"          "LOC107389177"   "rftn2"          "LOC107390197"   "hnf4a"          "yipf2"          "LOC107392149"   "ushbp1"         "LOC107393170"   "LOC107394581"  
# [41] "tpmt"           "LOC107394987"   "LOC107395338"   "tsen54"         "tbc1d22b"       "mcoln1"         "LOC107373484"   "LOC107373826"   "LOC107373895"   "LOC107374198"  
# [51] "LOC107374592"   "Repeat_894"     "Repeat_1146"    "Repeat_2259"  

####### b. Male: Aging down/Mif up
M.AgeD.MifU <- list("Male Aging Down" = my.age.gs.M$Age_dwn,
                    "Male Mif Up"     = my.mif.gs.M$Mif_up)
my.Venn <- Venn(M.AgeD.MifU)

pdf(paste0(Sys.Date(),"_Male_Aging_Down_MIf_Up.pdf"))
plot(my.Venn, doWeights = TRUE, type = "circles",  show = list(Faces = FALSE))
dev.off()

test.M.AgeD.MifU <- fisher.test(matrix(c(49, 3308, 29, length(union(rownames(M.mif.res$Aging),rownames(M.mif.res$Mif)))-49-3308-29),2,2))
test.M.AgeD.MifU$p.value
# [1] 2.136057e-23

M.ADMU_g <- intersect(my.age.gs.M$Age_dwn,my.mif.gs.M$Mif_up)
# [1] "LOC107373705"  "LOC107376506"  "LOC107379789"  "abcd2"         "LOC107373279"  "prss12"        "cdh23"         "LOC107375660"  "rab23"         "prr18"        
# [11] "slc16a9"       "xpo6"          "rai1"          "LOC107376475"  "LOC107377236"  "pcdh7"         "robo2"         "tenm4"         "csgr06h3orf52" "LOC107380229" 
# [21] "fam168a"       "mfap1"         "LOC107382811"  "LOC107384204"  "cntn3"         "LOC107386508"  "LOC107386606"  "LOC107386615"  "cnih2"         "slc7a3"       
# [31] "gabrg2"        "LOC107386855"  "ndrg4"         "tspan7"        "LOC107389842"  "LOC107390240"  "LOC107391209"  "odc1"          "dtnb"          "LOC107392143" 
# [41] "klf11"         "disp2"         "mn1"           "LOC107394594"  "LOC107396829"  "LOC107396859"  "slc35f3"       "LOC107373001"  "Repeat_1080"  
################

intersect(F.AUMD_g, M.AUMD_g)
# [1] "LOC107389177" "LOC107394987" "LOC107373895"
# LOC107373895 heat shock 70 kDa protein 1 [ Nothobranchius furzeri (turquoise killifish) ]
# LOC107389177 tgm2l transglutaminase 2, like [ Nothobranchius furzeri (turquoise killifish) ]
# LOC107394987 tlcd3ba TLC domain containing 3Ba [ Nothobranchius furzeri (turquoise killifish) ]
# 

intersect(F.ADMU_g, M.ADMU_g)
# [1] "LOC107379789" "LOC107394594"
# LOC107379789 si:dkey-237h12.3 teneurin-3 [ Nothobranchius furzeri (turquoise killifish) ]
# LOC107394594 ictacalcin [ Nothobranchius furzeri (turquoise killifish) ]


#####################################################################


#######################
sink(file = paste(Sys.Date(),"_GSEA_Venn_bulk_killi_BrainAging_Mifepristone_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()
