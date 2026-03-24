setwd('/Volumes/BB_HQ_3/Brain_aging_Simons/scRNAseq/snRNAseq_Brain_Aging_Meta/Human_alignment/')
options(stringsAsFactors = FALSE)

# 2024-03-12
# Parse Killifish/Human BLAST alignments for homology mapping

#########################################################################################################################
# read in result table and add in column names
blast.rev.res <- read.table('2024-03-11_HumanGencode_Nfur2015_BestHits_1e-3_REVERSE.txt', sep = "\t", header = F)
colnames(blast.rev.res) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
blast.rev.res$Human_Symbol    <- unlist(lapply(strsplit(blast.rev.res$qseqid,"|", fixed = TRUE),'[',7))

# read in killifish gene annotation table
killi.annot    <- read.csv('GCF_001465895.1_Nfu_20140520_feature_table.txt', header = F, sep = "\t")
killi.annot.cl <- killi.annot[,c(1, 11, 14, 15)]
colnames(killi.annot.cl) <- c("Type", "Nfur_Accession", "Nfur_Description", "Nfur_Symbol")

# select on protein only
killi.annot.cl <- killi.annot.cl[killi.annot.cl$Type %in% "CDS",]

##### Merge Info from killi annot and BLAST results
killi.hsa.rev <- merge(killi.annot.cl, blast.rev.res, by.x = "Nfur_Accession", by.y = "sseqid")

# select columns and reorder
colnames(killi.hsa.rev)
my.col.order <- c("Nfur_Accession", "Nfur_Description", "Nfur_Symbol", "Human_Symbol", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore") 
killi.hsa.rev.cl <- killi.hsa.rev[,my.col.order]

write.table(killi.hsa.rev.cl, file = paste0(Sys.Date(),"_Human_Best_BLAST_hit_to_Killifish_Annotated_1e-3_withStats.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

##### Loop over table to retain only best hit
dim(killi.hsa.rev.cl) # [1] 104054     14

length(unique(killi.hsa.rev.cl$Nfur_Accession)) # 19621

killi.hsa.orth.rev <- unique(killi.hsa.rev.cl[,c("Nfur_Accession", "Nfur_Description", "Nfur_Symbol", "Human_Symbol") ])

write.table(killi.hsa.orth.rev, file = paste0(Sys.Date(),"_Human_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-3_Minimal_HOMOLOGY_TABLE_REV.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
