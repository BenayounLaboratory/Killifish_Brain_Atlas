FqDir="/Volumes/BB_Home_HQ2/SIngle_Cell_analysis/ATAC/2022-04-15_Killi_Brain_Aging_Bulk_ATACseq/NGMerge"
BowtieDir="/Volumes/BB_Home_HQ2/SIngle_Cell_analysis/ATAC/2022-04-15_Killi_Brain_Aging_Bulk_ATACseq/BOWTIE2"
BOWTIE2_INDEX_KILLI="/Volumes/BB_Home_HQ2/Killifish_genome_variation/Genome/Process_data/GCF_001465895.1_Nfu_20140520_genomic"

for f in $(find $FqDir -name '*_NGmerge_1.fastq.gz')
do
    f2=$(basename "${f}" | sed 's/_NGmerge_1\.fastq\.gz/_NGmerge_2\.fastq\.gz/g');
    inf2="${FqDir}/${f2}"

    of=$(basename "${f}" | sed 's/_NGmerge_1\.fastq\.gz/\.bam/g');
    oFname="${BowtieDir}/${of}"

    bowtie2 --sensitive --no-discordant --no-mixed  -k 2 -p 2 -X 2500 -x $BOWTIE2_INDEX_KILLI -1 $f -2 $inf2 | samtools view -b -S - > $oFname

done



###clean_BAM_PE_v2.sh Killi2015 . ../..