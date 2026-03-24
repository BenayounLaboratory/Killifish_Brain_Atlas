# ### use trimgalore to recapitulate the behavior of fastx_trimmer (doesn't run on newer macos)

# ### cut 9 based on bias
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_F_Brain_c1_O_CTL1_1.fq.gz   ../FASTQ/GRZ_F_Brain_c1_O_CTL1_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_F_Brain_c1_O_CTL2_1.fq.gz   ../FASTQ/GRZ_F_Brain_c1_O_CTL2_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_F_Brain_c1_O_CTL3_1.fq.gz   ../FASTQ/GRZ_F_Brain_c1_O_CTL3_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_F_Brain_c1_O_MIF1_1.fq.gz   ../FASTQ/GRZ_F_Brain_c1_O_MIF1_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_F_Brain_c1_O_MIF2_1.fq.gz   ../FASTQ/GRZ_F_Brain_c1_O_MIF2_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_F_Brain_c1_Y_CTL1_1.fq.gz   ../FASTQ/GRZ_F_Brain_c1_Y_CTL1_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_F_Brain_c1_Y_CTL2_1.fq.gz   ../FASTQ/GRZ_F_Brain_c1_Y_CTL2_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_F_Brain_c1_Y_CTL3_1.fq.gz   ../FASTQ/GRZ_F_Brain_c1_Y_CTL3_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_M_Brain_c2_O_CTL1_1.fq.gz   ../FASTQ/GRZ_M_Brain_c2_O_CTL1_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_M_Brain_c2_O_CTL2_1.fq.gz   ../FASTQ/GRZ_M_Brain_c2_O_CTL2_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_M_Brain_c2_O_CTL3_1.fq.gz   ../FASTQ/GRZ_M_Brain_c2_O_CTL3_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_M_Brain_c2_O_MIF1_1.fq.gz   ../FASTQ/GRZ_M_Brain_c2_O_MIF1_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_M_Brain_c2_O_MIF2_1.fq.gz   ../FASTQ/GRZ_M_Brain_c2_O_MIF2_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_M_Brain_c2_O_MIF3_1.fq.gz   ../FASTQ/GRZ_M_Brain_c2_O_MIF3_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_M_Brain_c2_Y_CTL1_1.fq.gz   ../FASTQ/GRZ_M_Brain_c2_Y_CTL1_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_M_Brain_c2_Y_CTL2_1.fq.gz   ../FASTQ/GRZ_M_Brain_c2_Y_CTL2_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_M_Brain_c2_Y_CTL3_1.fq.gz   ../FASTQ/GRZ_M_Brain_c2_Y_CTL3_2.fq.gz
trim_galore --stringency 15 --clip_R1 9 --clip_R2 9 --paired ../FASTQ/GRZ_M_Brain_c2_Y_CTL4_1.fq.gz   ../FASTQ/GRZ_M_Brain_c2_Y_CTL4_2.fq.gz

