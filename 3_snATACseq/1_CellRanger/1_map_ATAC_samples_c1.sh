~/Softwares/cellranger-atac-2.1.0/bin/cellranger-atac count --id Killi_Brain_YF1_FishTEDB_NR \
                 --reference=/home/benayoun/Softwares/cellranger-atac-2.1.0/Reference/GCF_001465895_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-04-19_Brain_killi_scATACseq/FASTQ \
                 --sample Young_Female_CKDL220008676-1a-SI_3A_A1_HK52VDSX3  \
                 --localmem 64 \
                 --localcores 8

~/Softwares/cellranger-atac-2.1.0/bin/cellranger-atac count --id Killi_Brain_OF1_FishTEDB_NR \
                 --reference=/home/benayoun/Softwares/cellranger-atac-2.1.0/Reference/GCF_001465895_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-04-19_Brain_killi_scATACseq/FASTQ \
                 --sample Old_Female_CKDL220008676-1a-SI_3A_B1_HK52VDSX3  \
                 --localmem 64 \
                 --localcores 8

~/Softwares/cellranger-atac-2.1.0/bin/cellranger-atac count --id Killi_Brain_YM1_FishTEDB_NR \
                 --reference=/home/benayoun/Softwares/cellranger-atac-2.1.0/Reference/GCF_001465895_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-04-19_Brain_killi_scATACseq/FASTQ \
                 --sample Young_Male_CKDL220008676-1a-SI_3A_C1_HK52VDSX3  \
                 --localmem 64 \
                 --localcores 8

~/Softwares/cellranger-atac-2.1.0/bin/cellranger-atac count --id Killi_Brain_OM1_FishTEDB_NR \
                 --reference=/home/benayoun/Softwares/cellranger-atac-2.1.0/Reference/GCF_001465895_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-04-19_Brain_killi_scATACseq/FASTQ \
                 --sample Old_Male_CKDL220008676-1a-SI_3A_D1_HK52VDSX3  \
                 --localmem 64 \
                 --localcores 8


##########################################################################################
~/Softwares/cellranger-atac-2.1.0/bin/cellranger-atac aggr --id=Killi_Brain_ATAC_Set1 \
                 --reference=/home/benayoun/Softwares/cellranger-atac-2.1.0/Reference/GCF_001465895_FishTEDB_NR \
                 --csv=killi_Brain_1.csv \
                 --normalize=none        \
                 --localmem 64 \
                 --localcores 8

