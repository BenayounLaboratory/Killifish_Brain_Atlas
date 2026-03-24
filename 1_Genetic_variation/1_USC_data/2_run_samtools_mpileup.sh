#http://www.htslib.org/doc/samtools-mpileup.html
#https://samtools.github.io/bcftools/howtos/variant-calling.html //https://www.biostars.org/p/378305/

### samtools mpileup -uf ref.fa aln.bam | bcftools call -mv > var.raw.vcf
### bcftools filter -s LowQual -e '%QUAL<20 || DP>100' var.raw.vcf  > var.flt.vcf


#-C, --adjust-MQ INT
# Coefficient for downgrading mapping quality for reads containing excessive mismatches. Given a read with a phred-scaled probability q
# of being generated from the mapped position, the new mapping quality is about sqrt((INT-q)/INT)*INT. 
# A zero value disables this functionality; if enabled, the recommended value for BWA is 50. [0]

bcftools mpileup -Ou -C 50 -q 20 -Q 25 -d 8000 -f ../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna G2_CSFP200009003-1a_HNV2GDSXY_L3.srt.markDup.bam | bcftools call -mv -Oz -o G2_CSFP200009003-1a_HNV2GDSXY_L3.srt.markDup.bcftools.calls_v2.vcf.gz
bcftools mpileup -Ou -C 50 -q 20 -Q 25 -d 8000 -f ../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna G1_CSFP200009002-1a_HNV2GDSXY_L3.srt.markDup.bam | bcftools call -mv -Oz -o G1_CSFP200009002-1a_HNV2GDSXY_L3.srt.markDup.bcftools.calls_v2.vcf.gz
bcftools mpileup -Ou -C 50 -q 20 -Q 25 -d 8000 -f ../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna Z1_CSFP200009004-1a_HNV2GDSXY_L3.srt.markDup.bam | bcftools call -mv -Oz -o Z1_CSFP200009004-1a_HNV2GDSXY_L3.srt.markDup.bcftools.calls_v2.vcf.gz
bcftools mpileup -Ou -C 50 -q 20 -Q 25 -d 8000 -f ../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna Z2_CSFP200009005-1a_HNV2GDSXY_L3.srt.markDup.bam | bcftools call -mv -Oz -o Z2_CSFP200009005-1a_HNV2GDSXY_L3.srt.markDup.bcftools.calls_v2.vcf.gz
#### bcftools view -i '%QUAL>=20' calls.bcf
