# index killi genome
###bwa index GCF_001465895.1_Nfu_20140520_genomic.fna

#########
# for each of them run qsub
#for f in $(find "../FASTQ" -name '*_1.fq')
#do
#       
#    f2=$(basename "${f}" | sed 's/_1\.fq/_2\.fq/g');
#	inf2="../FASTQ/${f2}"
#	
#	of1="../FASTQ/$(basename "${f}"  | sed 's/fq/sai/g')";
#	of2="../FASTQ/$(basename "${f2}" | sed 's/fq/sai/g')";
#    
#    bwa aln -t 4 ../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna $f > $of1
#    bwa aln -t 4 ../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna $inf2 > $of2
#    
#    of=$(basename "${f}" | sed 's/_1\.fq/\.sam/g');
#    oFname="./${of}"
#    
#    bwa sampe ../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna $of1 $of2 $f $inf2 > $oFname
#    
#done

#samtools view -bS G2_CSFP200009003-1a_HNV2GDSXY_L3.sam > G2_CSFP200009003-1a_HNV2GDSXY_L3.bam
#samtools view -bS G1_CSFP200009002-1a_HNV2GDSXY_L3.sam > G1_CSFP200009002-1a_HNV2GDSXY_L3.bam
#samtools view -bS Z1_CSFP200009004-1a_HNV2GDSXY_L3.sam > Z1_CSFP200009004-1a_HNV2GDSXY_L3.bam
#samtools view -bS Z2_CSFP200009005-1a_HNV2GDSXY_L3.sam > Z2_CSFP200009005-1a_HNV2GDSXY_L3.bam

#samtools sort G2_CSFP200009003-1a_HNV2GDSXY_L3.bam G2_CSFP200009003-1a_HNV2GDSXY_L3.srt
#samtools sort G1_CSFP200009002-1a_HNV2GDSXY_L3.bam G1_CSFP200009002-1a_HNV2GDSXY_L3.srt
#samtools sort Z1_CSFP200009004-1a_HNV2GDSXY_L3.bam Z1_CSFP200009004-1a_HNV2GDSXY_L3.srt
#samtools sort Z2_CSFP200009005-1a_HNV2GDSXY_L3.bam Z2_CSFP200009005-1a_HNV2GDSXY_L3.srt

# samtools index G2_CSFP200009003-1a_HNV2GDSXY_L3.srt.bam
# samtools index G1_CSFP200009002-1a_HNV2GDSXY_L3.srt.bam
# samtools index Z1_CSFP200009004-1a_HNV2GDSXY_L3.srt.bam
# samtools index Z2_CSFP200009005-1a_HNV2GDSXY_L3.srt.bam

gatk MarkDuplicates -I G2_CSFP200009003-1a_HNV2GDSXY_L3.srt.bam -M G2_CSFP200009003-1a_HNV2GDSXY_L3.srt.bam.DupMetrics -O G2_CSFP200009003-1a_HNV2GDSXY_L3.srt.markDup.bam --CREATE_INDEX true
gatk MarkDuplicates -I G1_CSFP200009002-1a_HNV2GDSXY_L3.srt.bam -M G1_CSFP200009002-1a_HNV2GDSXY_L3.srt.bam.DupMetrics -O G1_CSFP200009002-1a_HNV2GDSXY_L3.srt.markDup.bam --CREATE_INDEX true
gatk MarkDuplicates -I Z1_CSFP200009004-1a_HNV2GDSXY_L3.srt.bam -M Z1_CSFP200009004-1a_HNV2GDSXY_L3.srt.bam.DupMetrics -O Z1_CSFP200009004-1a_HNV2GDSXY_L3.srt.markDup.bam --CREATE_INDEX true
gatk MarkDuplicates -I Z2_CSFP200009005-1a_HNV2GDSXY_L3.srt.bam -M Z2_CSFP200009005-1a_HNV2GDSXY_L3.srt.bam.DupMetrics -O Z2_CSFP200009005-1a_HNV2GDSXY_L3.srt.markDup.bam --CREATE_INDEX true
