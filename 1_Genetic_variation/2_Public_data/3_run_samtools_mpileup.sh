#http://www.htslib.org/doc/samtools-mpileup.html
#https://samtools.github.io/bcftools/howtos/variant-calling.html //https://www.biostars.org/p/378305/

### samtools mpileup -uf ref.fa aln.bam | bcftools call -mv > var.raw.vcf
### bcftools filter -s LowQual -e '%QUAL<20 || DP>100' var.raw.vcf  > var.flt.vcf


#-C, --adjust-MQ INT
# Coefficient for downgrading mapping quality for reads containing excessive mismatches. Given a read with a phred-scaled probability q
# of being generated from the mapped position, the new mapping quality is about sqrt((INT-q)/INT)*INT. 
# A zero value disables this functionality; if enabled, the recommended value for BWA is 50. [0]

bcftools mpileup -Ou -C 50 -q 20 -Q 25 -d 8000 -f ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna PRJEB5837_GRZ_M004_genomic.mergedLanes.markDup.bam      |  bcftools call -mv -Oz -o PRJEB5837_GRZ_M004_genomic.mergedLanes.bcftools.calls.vcf.gz  
bcftools mpileup -Ou -C 50 -q 20 -Q 25 -d 8000 -f ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna PRJEB5837_GRZ_M005_genomic.mergedLanes.markDup.bam      |  bcftools call -mv -Oz -o PRJEB5837_GRZ_M005_genomic.mergedLanes.bcftools.calls.vcf.gz  
bcftools mpileup -Ou -C 50 -q 20 -Q 25 -d 8000 -f ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna PRJEB5837_MZM0403_M006_genomic.mergedLanes.markDup.bam  |  bcftools call -mv -Oz -o PRJEB5837_MZM0403_M006_genomic.mergedLanes.bcftools.calls.vcf.gz  
bcftools mpileup -Ou -C 50 -q 20 -Q 25 -d 8000 -f ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna PRJEB5837_MZM0403_M007_genomic.mergedLanes.markDup.bam  |  bcftools call -mv -Oz -o PRJEB5837_MZM0403_M007_genomic.mergedLanes.bcftools.calls.vcf.gz 
bcftools mpileup -Ou -C 50 -q 20 -Q 25 -d 8000 -f ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna PRJEB5837_MZM0410_M012_genomic.mergedLanes.markDup.bam  |  bcftools call -mv -Oz -o PRJEB5837_MZM0410_M012_genomic.mergedLanes.bcftools.calls.vcf.gz 
bcftools mpileup -Ou -C 50 -q 20 -Q 25 -d 8000 -f ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna PRJEB5837_MZM0410_M013_genomic.mergedLanes.markDup.bam  |  bcftools call -mv -Oz -o PRJEB5837_MZM0410_M013_genomic.mergedLanes.bcftools.calls.vcf.gz  
bcftools mpileup -Ou -C 50 -q 20 -Q 25 -d 8000 -f ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna PRJEB5837_MZZW0701_M001_genomic.mergedLanes.markDup.bam |  bcftools call -mv -Oz -o PRJEB5837_MZZW0701_M001_genomic.mergedLanes.bcftools.calls.vcf.gz  
bcftools mpileup -Ou -C 50 -q 20 -Q 25 -d 8000 -f ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna PRJEB5837_MZZW0701_M002_genomic.mergedLanes.markDup.bam |  bcftools call -mv -Oz -o PRJEB5837_MZZW0701_M002_genomic.mergedLanes.bcftools.calls.vcf.gz  

bcftools mpileup -Ou -C 50 -q 20 -Q 25 -d 8000 -f ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna PRJNA221357_SRR1261472_GRZ_GRZ300.markDup.bam              |  bcftools call -mv -Oz -o PRJNA221357_SRR1261472_GRZ_GRZ300.bcftools.calls.vcf.gz          
bcftools mpileup -Ou -C 50 -q 20 -Q 25 -d 8000 -f ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna PRJNA221357_SRR1261474_GRZ_GRZ340.markDup.bam              |  bcftools call -mv -Oz -o PRJNA221357_SRR1261474_GRZ_GRZ340.bcftools.calls.vcf.gz          
bcftools mpileup -Ou -C 50 -q 20 -Q 25 -d 8000 -f ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna PRJNA221357_SRR1261476_GRZ_BGI170.markDup.bam              |  bcftools call -mv -Oz -o PRJNA221357_SRR1261476_GRZ_BGI170.bcftools.calls.vcf.gz          
bcftools mpileup -Ou -C 50 -q 20 -Q 25 -d 8000 -f ../../Process_data/GCF_001465895.1_Nfu_20140520_genomic.fna PRJNA221357_SRR1261480_GRZ_P0FA.markDup.bam                |  bcftools call -mv -Oz -o PRJNA221357_SRR1261480_GRZ_P0FA.bcftools.calls.vcf.gz            
