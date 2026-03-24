

#### DEBUG
# mif.cts <- f.mif.cts
# sex     <- "Females"

# mif.cts <- m.mif.cts
# sex     <- "Males"

process_mif_data <- function (mif.cts, sex = "Females") {
  
  # remove column with gene names since it was put in rownames
  mif.cts <- mif.cts[,-1]
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_Bulk_GRZ_Brain_Aging_Mif_SVA_", sex)
  
  ###################################################################################################################
  # 0a. create covariates matrix for analysis
  
  # learn from column names
  my.group                                  <- rep(NA, ncol(mif.cts))
  my.group[grep("Y_CTL",colnames(mif.cts))] <- "Y_CTL"
  my.group[grep("O_CTL",colnames(mif.cts))] <- "O_CTL"
  my.group[grep("O_MIF",colnames(mif.cts))] <- "O_MIF"
  
  my.age                                 <- rep(NA, ncol(mif.cts))
  my.age[grep("Y_",colnames(mif.cts))]   <- "Young"
  my.age[grep("O_",colnames(mif.cts))]   <- "Old"
  
  ## make data frame of covariates
  mif.meta <- data.frame("Age"       = my.age,
                         "Group"     = my.group)
  
  rownames(mif.meta) <- colnames(mif.cts)
  ###################################################################################################################
  
  ###################################################################################################################
  # 0b. correct technical noise with SVA
  
  ###################################
  #######       Run SVA      #######
  
  # Set null and alternative models (ignore batch)
  mod1    = model.matrix(~ Group, data = mif.meta)
  n.sv.be = num.sv(mif.cts, mod1, method="be") # males is 2/Female is 1
  
  # apply SVAseq algortihm
  my.svseq = svaseq(as.matrix(mif.cts), mod1, n.sv=n.sv.be, constant = 0.1)
  
  # remove RIN and SV, preserve age and sex
  my.clean <- removeBatchEffect(log2(mif.cts + 0.1), 
                                covariates = cbind(my.svseq$sv),
                                design     = mod1)
  
  # delog and round data for DEseq2 processing
  my.filtered.sva <- round(2^my.clean-0.1)
  
  # output corrected counts
  write.table(my.filtered.sva, file = paste0(my.outprefix,"_SVA_corrected_counts_matrix.txt") , sep = "\t" , row.names = T, quote = F)
  
  # update for below
  mif.cts <- my.filtered.sva
  ###################################################################################################################
  
  
  #####################################################################################################################
  #### 1. Run DEseq2
  
  # see deseq2 vignette, remove genes without consistent expression (at least half samples)
  my.keep <- apply(mif.cts > 0, 1, sum) > ncol(mif.cts)/2
  
  # Now pull out the null/low expressed genes
  # Round values since multimapping was allowed
  my.filtered.matrix           <- round(mif.cts[my.keep,]) 

  # get matrix using age as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = my.filtered.matrix,
                                colData   = mif.meta,
                                design    = ~ Group)
  
  # run DESeq normalizations and export results
  dds.deseq <- DESeq(dds)
  
  # plot dispersion
  pdf(paste0(my.outprefix,"_dispersion_plot.pdf"))
  plotDispEsts(dds.deseq)
  dev.off()
  
  # get DESeq2 normalized expression value
  vst.cts.grz <- getVarianceStabilizedData(dds.deseq)
  
  # output counts
  my.out.ct.mat <- paste0(my.outprefix,"_VST_log2_counts_matrix.txt")
  write.table(vst.cts.grz, file = my.out.ct.mat , sep = "\t" , row.names = T, quote = F)
  
  ##### get color legend
  my.colors <- rep(NA, ncol(mif.cts))
  
  if (sex == "Females") {
    my.colors[my.group %in% "Y_CTL"] <- "deeppink"
    my.colors[my.group %in% "O_CTL"] <- "deeppink4"
    my.colors[my.group %in% "O_MIF"] <- "hotpink2"
    
  } else {
    my.colors[my.group %in% "Y_CTL"] <- "deepskyblue"
    my.colors[my.group %in% "O_CTL"] <- "deepskyblue4"
    my.colors[my.group %in% "O_MIF"] <- "lightskyblue"
    
  }
  

  
  # plot expression
  pdf(paste0(my.outprefix,"_expression_boxplot.pdf"))
  boxplot(vst.cts.grz, las = 2,col = my.colors)
  dev.off()
  
  # do MDS analysis
  mds.result <- cmdscale(1-cor(vst.cts.grz,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  x <- mds.result[, 1]
  y <- mds.result[, 2]
  
  pdf(paste0(my.outprefix,"_MDS_plot.pdf"), width = 5, height = 5)
  plot(x, y,
       xlab = "MDS dimension 1", ylab = "MDS dimension 2",
       main= paste(sex, "MIF brain aging"),
       cex=2, pch = 16, col = my.colors,
       las =1)
  points(x, y, cex = 2, pch = 1)
  legend("topleft",c("Y_CTL","O_CTL","O_MIF"), col = unique(my.colors), pch = 16, cex = 0.75, bty = 'n')
  dev.off()
  
  
  # PCA analysis
  my.pos.var <- apply(vst.cts.grz,1,var) > 0
  my.pca <- prcomp(t(vst.cts.grz[my.pos.var,]),scale = TRUE)
  x <- my.pca$x[,1]
  y <- my.pca$x[,2]
  
  my.summary <- summary(my.pca)
  
  pdf(paste0(my.outprefix,"_PCA_plot.pdf"), width = 5, height = 5)
  plot(x,y,
       cex=2, pch = 16, col = my.colors,
       xlab = paste('PC1 (', round(100*my.summary$importance[,1][2],1),"%)", sep=""),
       ylab = paste('PC2 (', round(100*my.summary$importance[,2][2],1),"%)", sep=""))
  points(x, y, cex = 2, pch = 1)
  legend("topleft",c("Y_CTL","O_CTL","O_MIF"), col = unique(my.colors), pch = 16, cex = 0.75, bty = 'n')
  dev.off()

  # extract gene significance by DEseq2
  res.mif <- results(dds.deseq, contrast = c("Group","O_MIF","O_CTL")) 
  res.age <- results(dds.deseq, contrast = c("Group","O_CTL","Y_CTL")) 
  
  # exclude genes with NA FDR value
  res.mif <- res.mif[!is.na(res.mif$padj),]
  res.age <- res.age[!is.na(res.age$padj),]
  
  ### get aging changes at FDR5
  genes.mif <- rownames(res.mif)[res.mif$padj < 0.05]
  genes.age <- rownames(res.age)[res.age$padj < 0.05]
  my.num.mif <- length(genes.mif) # 
  my.num.age <- length(genes.age) # 
  
  if (my.num.mif > 2) {
    pdf(paste0(my.outprefix,"_MIF_Heatmap_FDR5_GENES.pdf"), onefile = F, height = 10, width = 10)
    my.heatmap.title <- paste0("Brain Mif significant (FDR<5%), ", my.num.mif, " genes")
    pheatmap::pheatmap(vst.cts.grz[genes.mif,],
                       cluster_cols = F,
                       cluster_rows = T,
                       colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                       show_rownames = F, scale="row",
                       main = my.heatmap.title,
                       cellwidth = 15,
                       border    = NA,
                       cellheight = 0.1 )
    dev.off()
  }

  
  pdf(paste0(my.outprefix,"_AGING_Heatmap_FDR5_GENES.pdf"), onefile = F, height = 10, width = 10)
  my.heatmap.title <- paste0("Brain Aging significant (FDR<5%), ", my.num.age, " genes")
  pheatmap::pheatmap(vst.cts.grz[genes.age,],
                     cluster_cols = F,
                     cluster_rows = T,
                     colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                     show_rownames = F, scale="row",
                     main = my.heatmap.title,
                     cellwidth = 15,
                     border    = NA,
                     cellheight = 0.1 )
  dev.off()
  
  
  # output result tables of combined analysis to text files
  my.out.stats.mif <- paste0(my.outprefix,"_Mifepristone_all_genes_statistics.txt")
  write.table(res.mif, file = my.out.stats.mif , sep = "\t" , row.names = T, quote = F)
  
  my.out.stats.mif <- paste0(my.outprefix,"_AGING_all_genes_statistics.txt")
  write.table(res.age, file = my.out.stats.mif , sep = "\t" , row.names = T, quote = F)
  #####################################################################################################################
  
  #####################################################################################################################
  #### 2. Compare conditions pairwise
  
  merged.mif_age <- merge(data.frame(res.age), data.frame(res.mif), by = "row.names", suffixes = c(".age",".mif"))
  
  my.spear.cor <- cor.test(merged.mif_age$log2FoldChange.age,merged.mif_age$log2FoldChange.mif, method = 'spearman')
  my.rho       <- signif(my.spear.cor$estimate,3)
  
  pdf(paste0(my.outprefix,"AGING_vs_MIG_FC_scatterplot.pdf"))
  smoothScatter(merged.mif_age$log2FoldChange.age,merged.mif_age$log2FoldChange.mif, 
                xlim = c(-4,4), ylim = c(-4,4),
                xlab = paste("log2(FC) in ", sex, " with aging"),
                ylab = paste("log2(FC) in ", sex, " with Mif"  ),
                main = "Killifish Brain")
  abline(0,1, col = "grey", lty = "dashed")
  abline(h = 0, col = "red", lty = "dashed")
  abline(v = 0, col = "red", lty = "dashed")
  text(-3.9, 4, paste("Rho = ",my.rho), pos = 4)
  text(-3.9, 3.7, paste("p = ",signif(my.spear.cor$p.value,2)), pos = 4)
  dev.off()
  #####################################################################################
  
  return(list("Aging" = data.frame(res.age), "Mif" = data.frame(res.mif), "VST" = vst.cts.grz))
  
}
