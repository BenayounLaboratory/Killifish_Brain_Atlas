for f in $(find "." -name 'PRJEB5837_*mergedLanes.bam')
do
       
    srtp=$(basename "${f}" | sed 's/\.bam/\.srt/g');
	srtbam=$(basename "${f}" | sed 's/\.bam/\.srt.bam/g');
	
	samtools sort $f $srtp
	samtools index $srtbam
	
	dupMet=$(basename "${f}" | sed 's/\.bam/\.DupMetrics/g');
	nodupbam=$(basename "${f}" | sed 's/\.bam/\.markDup\.bam/g');

	gatk MarkDuplicates -I $srtbam -M $dupMet -O $nodupbam --CREATE_INDEX true

done