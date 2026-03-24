BOWTIE2_INDEX_KILLI="/Volumes/BB_Home_HQ2/Killifish_genome_variation/Genome/Process_data/GCF_001465895.1_Nfu_20140520_genomic"

bowtie2 --sensitive -k 2 -p 1 -x $BOWTIE2_INDEX_KILLI -U ../FASTQ/SRR1556403_GRZ_YM_Brain_H3K4me3_ChIPseq.fastq.gz | samtools view -b -S - > SRR1556403_GRZ_YM_Brain_H3K4me3_ChIPseq.bam
bowtie2 --sensitive -k 2 -p 1 -x $BOWTIE2_INDEX_KILLI -U ../FASTQ/SRR1556404_GRZ_YM_Brain_Input.fastq.gz           | samtools view -b -S - > SRR1556404_GRZ_YM_Brain_Input.bam       


## [samopen] SAM header is present: 5897 sequences.
## 34781631 reads; of these:
##   34781631 (100.00%) were unpaired; of these:
##     2778415 (7.99%) aligned 0 times
##     19474539 (55.99%) aligned exactly 1 time
##     12528677 (36.02%) aligned >1 times
## 92.01% overall alignment rate
## [samopen] SAM header is present: 5897 sequences.
## 45382183 reads; of these:
##   45382183 (100.00%) were unpaired; of these:
##     2996533 (6.60%) aligned 0 times
##     17985486 (39.63%) aligned exactly 1 time
##     24400164 (53.77%) aligned >1 times
## 93.40% overall alignment rate