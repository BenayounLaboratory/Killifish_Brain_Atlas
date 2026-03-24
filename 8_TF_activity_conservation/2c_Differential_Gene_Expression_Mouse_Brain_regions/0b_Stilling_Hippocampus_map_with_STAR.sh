export GENIN="/Users/benayoun/Softwares/STAR-2.7.7a/genomes/MM39"

STAR --genomeDir $GENIN --readFilesIn ../FASTQ/GSM1517299_3m_Hippocampus_RNAseq_1.fastq.gz    ../FASTQ/GSM1517299_3m_Hippocampus_RNAseq_2.fastq.gz   --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix GSM1517299_3m_Hippocampus_RNAseq_STAR_    
STAR --genomeDir $GENIN --readFilesIn ../FASTQ/GSM1517300_3m_Hippocampus_RNAseq_1.fastq.gz    ../FASTQ/GSM1517300_3m_Hippocampus_RNAseq_2.fastq.gz   --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix GSM1517300_3m_Hippocampus_RNAseq_STAR_    
STAR --genomeDir $GENIN --readFilesIn ../FASTQ/GSM1517301_3m_Hippocampus_RNAseq_1.fastq.gz    ../FASTQ/GSM1517301_3m_Hippocampus_RNAseq_2.fastq.gz   --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix GSM1517301_3m_Hippocampus_RNAseq_STAR_    
STAR --genomeDir $GENIN --readFilesIn ../FASTQ/GSM1517302_29m_Hippocampus_RNAseq_1.fastq.gz   ../FASTQ/GSM1517302_29m_Hippocampus_RNAseq_2.fastq.gz  --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix GSM1517302_29m_Hippocampus_RNAseq_STAR_   
STAR --genomeDir $GENIN --readFilesIn ../FASTQ/GSM1517303_29m_Hippocampus_RNAseq_1.fastq.gz   ../FASTQ/GSM1517303_29m_Hippocampus_RNAseq_2.fastq.gz  --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix GSM1517303_29m_Hippocampus_RNAseq_STAR_   
STAR --genomeDir $GENIN --readFilesIn ../FASTQ/GSM1517304_29m_Hippocampus_RNAseq_1.fastq.gz   ../FASTQ/GSM1517304_29m_Hippocampus_RNAseq_2.fastq.gz  --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --readFilesCommand gzcat --runThreadN 1 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate  --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix GSM1517304_29m_Hippocampus_RNAseq_STAR_   



featureCounts -t exon -D 1500 -p --primary -T 1 -s 0 -a /Users/benayoun/Softwares/Genomes/mm39_refGene.gtf -o 2025-07-30_Stilling_Hippocampus_Aging_counts_MM39.txt \
 	GSM1517299_3m_Hippocampus_RNAseq_STAR_Aligned.sortedByCoord.out.bam    \
 	GSM1517300_3m_Hippocampus_RNAseq_STAR_Aligned.sortedByCoord.out.bam    \
 	GSM1517301_3m_Hippocampus_RNAseq_STAR_Aligned.sortedByCoord.out.bam    \
 	GSM1517302_29m_Hippocampus_RNAseq_STAR_Aligned.sortedByCoord.out.bam    \
 	GSM1517303_29m_Hippocampus_RNAseq_STAR_Aligned.sortedByCoord.out.bam    \
 	GSM1517304_29m_Hippocampus_RNAseq_STAR_Aligned.sortedByCoord.out.bam    \
