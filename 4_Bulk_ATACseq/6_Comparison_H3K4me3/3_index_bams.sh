# sort bams
samtools sort SRR1556403_GRZ_YM_Brain_H3K4me3_ChIPseq.bam SRR1556403_GRZ_YM_Brain_H3K4me3_ChIPseq.srt
samtools sort SRR1556404_GRZ_YM_Brain_Input.bam           SRR1556404_GRZ_YM_Brain_Input.srt

samtools index SRR1556403_GRZ_YM_Brain_H3K4me3_ChIPseq.srt.bam
samtools index SRR1556404_GRZ_YM_Brain_Input.srt.bam
