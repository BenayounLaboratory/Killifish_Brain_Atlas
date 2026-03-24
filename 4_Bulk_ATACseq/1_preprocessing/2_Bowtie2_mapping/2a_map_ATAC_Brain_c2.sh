FqDir="/home/benayoun/Projects/Killi_ATAC/FQ_c2"
BowtieDir="/home/benayoun/Projects/Killi_ATAC/BOWTIE2"
BOWTIE2_INDEX_KILLI="/home/benayoun/Projects/Killi_ATAC/BOWTIE2_IX/GCF_001465895.1_Nfu_20140520_genomic"

for f in $(find $FqDir -name '*_NGmerge_1.fastq.gz')
do
    f2=$(basename "${f}" | sed 's/_NGmerge_1\.fastq\.gz/_NGmerge_2\.fastq\.gz/g');
    inf2="${FqDir}/${f2}"

    of=$(basename "${f}" | sed 's/_NGmerge_1\.fastq\.gz/\.bam/g');
    oFname="${BowtieDir}/${of}"
    
    logf=$(basename "${f}" | sed 's/NGmerge_1\.fastq\.gz/stats\.file\.txt/g');
    logFname="${BowtieDir}/${logf}"

    bowtie2 --sensitive --no-discordant --no-mixed  -k 2 -p 6 -X 2500 -x $BOWTIE2_INDEX_KILLI -1 $f -2 $inf2 2> $logFname | samtools view -b -S - > $oFname

done

