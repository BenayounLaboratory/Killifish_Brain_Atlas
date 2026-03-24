# GRZ samples (redo since issues on detected cells)
~/Softwares/cellranger-6.1.2/bin/cellranger count --id GRZ_5w_M_1_SIGAH1_FishTEDB_NR_v2 \
                 --transcriptome ~/Softwares/cellranger-6.1.2/Reference/GCF_001465895.1_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-06-02_Brain_killi_snRNAseq_cohort1 \
                 --sample GRZ_5w_M_1_SIGAH1  \
                 --force-cells 8000 \
                 --localmem 128 \
                 --localcores 16 \
                 --include-introns

~/Softwares/cellranger-6.1.2/bin/cellranger count --id GRZ_10w_M_1_SIGAH2_FishTEDB_NR_v2 \
                 --transcriptome ~/Softwares/cellranger-6.1.2/Reference/GCF_001465895.1_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-06-02_Brain_killi_snRNAseq_cohort1 \
                 --sample GRZ_10w_M_1_SIGAH2  \
                 --force-cells 8000 \
                 --localmem 128 \
                 --localcores 16 \
                 --include-introns

~/Softwares/cellranger-6.1.2/bin/cellranger count --id GRZ_15w_M_1_SIGAH3_FishTEDB_NR_v2 \
                 --transcriptome ~/Softwares/cellranger-6.1.2/Reference/GCF_001465895.1_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-06-02_Brain_killi_snRNAseq_cohort1 \
                 --sample GRZ_15w_M_1_SIGAH3  \
                 --force-cells 8000 \
                 --localmem 128 \
                 --localcores 16 \
                 --include-introns

~/Softwares/cellranger-6.1.2/bin/cellranger count --id GRZ_5w_F_1_SIGAH4_FishTEDB_NR_v2 \
                 --transcriptome ~/Softwares/cellranger-6.1.2/Reference/GCF_001465895.1_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-06-02_Brain_killi_snRNAseq_cohort1 \
                 --sample GRZ_5w_F_1_SIGAH4  \
                 --force-cells 8000 \
                 --localmem 128 \
                 --localcores 16 \
                 --include-introns

~/Softwares/cellranger-6.1.2/bin/cellranger count --id GRZ_10w_F_1_SIGAH5_FishTEDB_NR_v2 \
                 --transcriptome ~/Softwares/cellranger-6.1.2/Reference/GCF_001465895.1_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-06-02_Brain_killi_snRNAseq_cohort1 \
                 --sample GRZ_10w_F_1_SIGAH5  \
                 --force-cells 8000 \
                 --localmem 128 \
                 --localcores 16 \
                 --include-introns

~/Softwares/cellranger-6.1.2/bin/cellranger count --id GRZ_15w_F_1_SIGAH6_FishTEDB_NR_v2 \
                 --transcriptome ~/Softwares/cellranger-6.1.2/Reference/GCF_001465895.1_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-06-02_Brain_killi_snRNAseq_cohort1 \
                 --sample GRZ_15w_F_1_SIGAH6  \
                 --force-cells 8000 \
                 --localmem 128 \
                 --localcores 16 \
                 --include-introns


