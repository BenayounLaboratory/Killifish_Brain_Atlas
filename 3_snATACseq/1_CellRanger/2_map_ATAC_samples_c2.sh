~/Softwares/cellranger-atac-2.1.0/bin/cellranger-atac count --id Killi_Brain_YF2_FishTEDB_NR \
                 --reference=/home/benayoun/Softwares/cellranger-atac-2.1.0/Reference/GCF_001465895_FishTEDB_NR \
                 --fastqs /mnt/data3/Data/2023-01-18_Brain_killi_scATACseq_cohort2/FASTQ \
                 --sample Young_Female_2_CKDL220034403-1A_HMMFCDSX5  \
                 --localmem 64 \
                 --localcores 8

~/Softwares/cellranger-atac-2.1.0/bin/cellranger-atac count --id Killi_Brain_OF2_FishTEDB_NR \
                 --reference=/home/benayoun/Softwares/cellranger-atac-2.1.0/Reference/GCF_001465895_FishTEDB_NR \
                 --fastqs /mnt/data3/Data/2023-01-18_Brain_killi_scATACseq_cohort2/FASTQ \
                 --sample Old_Female_2_CKDL220034403-1A_HMMFCDSX5  \
                 --localmem 64 \
                 --localcores 8

~/Softwares/cellranger-atac-2.1.0/bin/cellranger-atac count --id Killi_Brain_YM2_FishTEDB_NR \
                 --reference=/home/benayoun/Softwares/cellranger-atac-2.1.0/Reference/GCF_001465895_FishTEDB_NR \
                 --fastqs /mnt/data3/Data/2023-01-18_Brain_killi_scATACseq_cohort2/FASTQ \
                 --sample Young_Male_2_CKDL220034403-1A_HMMFCDSX5_S2  \
                 --localmem 64 \
                 --localcores 8

~/Softwares/cellranger-atac-2.1.0/bin/cellranger-atac count --id Killi_Brain_OM2_FishTEDB_NR \
                 --reference=/home/benayoun/Softwares/cellranger-atac-2.1.0/Reference/GCF_001465895_FishTEDB_NR \
                 --fastqs /mnt/data3/Data/2023-01-18_Brain_killi_scATACseq_cohort2/FASTQ \
                 --sample Old_Male_2_CKDL220034403-1A_HMMFCDSX5  \
                 --localmem 64 \
                 --localcores 8


##########################################################################################
~/Softwares/cellranger-atac-2.1.0/bin/cellranger-atac aggr --id=Killi_Brain_ATAC_Set1 \
                 --reference=/home/benayoun/Softwares/cellranger-atac-2.1.0/Reference/GCF_001465895_FishTEDB_NR \
                 --csv=killi_Brain_2.csv \
                 --normalize=none        \
                 --localmem 64 \
                 --localcores 8

