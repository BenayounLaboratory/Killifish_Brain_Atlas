~/Softwares/cellranger-6.1.2/bin/cellranger count --id GRZ_Brain_F_FishTEDB_NR \
                 --transcriptome ~/Softwares/cellranger-6.1.2/Reference/GCF_001465895.1_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-01-24_killifish_tissues_cohort3/FASTQ \
                 --sample GRZ_Brain_F_SIGAE8  \
                 --expect-cells 12000 \
                 --localmem 128 \
                 --localcores 16 \
                 --include-introns

~/Softwares/cellranger-6.1.2/bin/cellranger count --id GRZ_Brain_M_FishTEDB_NR \
                 --transcriptome ~/Softwares/cellranger-6.1.2/Reference/GCF_001465895.1_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-01-24_killifish_tissues_cohort3/FASTQ \
                 --sample GRZ_Brain_M_SIGAE7  \
                 --expect-cells 12000 \
                 --localmem 128 \
                 --localcores 16 \
                 --include-introns
