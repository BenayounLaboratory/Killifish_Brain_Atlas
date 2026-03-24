for f in $(find "../FLI_data/" -name '*_1.fastq\.gz')
do
       
    f2=$(basename "${f}" | sed 's/_1\.fastq\.gz/_2\.fastq\.gz/g');
	inf2="../FLI_data/${f2}"
	
	of1="./$(basename "${f}"  | sed 's/fastq\.gz/sai/g')";
	of2="./$(basename "${f2}" | sed 's/fastq\.gz/sai/g')";
    
    bwa aln -t 4 ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna <(gzcat $f)    > $of1
    bwa aln -t 4 ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna <(gzcat $inf2) > $of2
    
    of=$(basename "${f}" | sed 's/_1\.fastq\.gz/\.sam/g');
    oFname="./${of}"
    bwa sampe ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna $of1 $of2 $f $inf2 > $oFname
     
    obam=$(basename "${f}" | sed 's/_1\.fastq\.gz/\.bam/g');
    samtools view -bS -F 4 $oFname   > $obam  # mapped


done

for f in $(find "../Stanford/" -name '*_1.fastq\.gz')
do
       
    f2=$(basename "${f}" | sed 's/_1\.fastq\.gz/_2\.fastq\.gz/g');
	inf2="../Stanford/${f2}"
	
	of1="./$(basename "${f}"  | sed 's/fastq\.gz/sai/g')";
	of2="./$(basename "${f2}" | sed 's/fastq\.gz/sai/g')";
    
    bwa aln -t 4 ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna <(gzcat $f)    > $of1
    bwa aln -t 4 ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna <(gzcat $inf2) > $of2
    
    of=$(basename "${f}" | sed 's/_1\.fastq\.gz/\.sam/g');
    oFname="./${of}"
    bwa sampe ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna $of1 $of2 $f $inf2 > $oFname
    
    obam=$(basename "${f}" | sed 's/_1\.fastq\.gz/\.bam/g');
    samtools view -bS -F 4 $oFname   > $obam  # mapped
    
done

#### bwa aln -t 4 ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna <(gunzip /Volumes/BB_Home_HQ2/Killifish_genome_variation/Genome/FASTQ_public/FLI_data/PRJEB5837_ERR878522_MZM-0403_M006_genomic_1.fastq.gz) > test.sai
### https://tldp.org/LDP/abs/html/process-sub.html

#samtools index SRR1718758_MZM0403.srt.bam
#samtools index SRR1718759_MZM0403.srt.bam
#samtools index SRR1261482_MZM0703.srt.bam
#
#gatk MarkDuplicates -I SRR1718758_MZM0403.srt.mp.bam -M SRR1718758_MZM0403.srt.mp.bam.DupMetrics -O SRR1718758_MZM0403.srt.mp.markDup.bam --CREATE_INDEX true
#gatk MarkDuplicates -I SRR1718759_MZM0403.srt.mp.bam -M SRR1718759_MZM0403.srt.mp.bam.DupMetrics -O SRR1718759_MZM0403.srt.mp.markDup.bam --CREATE_INDEX true
#gatk MarkDuplicates -I SRR1261482_MZM0703.srt.mp.bam -M SRR1261482_MZM0703.srt.mp.bam.DupMetrics -O SRR1261482_MZM0703.srt.mp.markDup.bam --CREATE_INDEX true
#
#
#samtools merge POOLED_SRR171875x_MZM0403.mp.bam SRR1718758_MZM0403.srt.mp.bam SRR1718759_MZM0403.srt.mp.bam
#
#gatk MarkDuplicates -I POOLED_SRR171875x_MZM0403.srt.mp.bam -M POOLED_SRR171875x_MZM0403.srt.mp.bam.DupMetrics -O POOLED_SRR171875x_MZM0403.srt.mp.markDup.bam --CREATE_INDEX true
#
#

