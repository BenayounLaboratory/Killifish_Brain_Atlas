#cp ../BWA/G1_CSFP200009002-1a_HNV2GDSXY_L3.srt.markDup.bcftools.calls_v2.vcf.gz .
#cp ../BWA/G2_CSFP200009003-1a_HNV2GDSXY_L3.srt.markDup.bcftools.calls_v2.vcf.gz .
#cp ../BWA/Z1_CSFP200009004-1a_HNV2GDSXY_L3.srt.markDup.bcftools.calls_v2.vcf.gz .
#cp ../BWA/Z2_CSFP200009005-1a_HNV2GDSXY_L3.srt.markDup.bcftools.calls_v2.vcf.gz .
#cp ../Stanford_LL/SRR1261482_MZM0703.bcftools.calls.vcf.gz  .
#cp ../Stanford_LL/POOLED_SRR171875x_MZM0403.bcftools.calls.vcf.gz .

cp ../FASTQ_public/BWA_public/PRJNA221357_SRR1261480_GRZ_P0FA.bcftools.calls.vcf.gz .
cp ../FASTQ_public/BWA_public/PRJNA221357_SRR1261476_GRZ_BGI170.bcftools.calls.vcf.gz .
cp ../FASTQ_public/BWA_public/PRJNA221357_SRR1261474_GRZ_GRZ340.bcftools.calls.vcf.gz .
cp ../FASTQ_public/BWA_public/PRJNA221357_SRR1261472_GRZ_GRZ300.bcftools.calls.vcf.gz .
cp ../FASTQ_public/BWA_public/PRJEB5837_MZZW0701_M002_genomic.mergedLanes.bcftools.calls.vcf.gz .
cp ../FASTQ_public/BWA_public/PRJEB5837_MZZW0701_M001_genomic.mergedLanes.bcftools.calls.vcf.gz .
cp ../FASTQ_public/BWA_public/PRJEB5837_MZM0410_M013_genomic.mergedLanes.bcftools.calls.vcf.gz .
cp ../FASTQ_public/BWA_public/PRJEB5837_MZM0410_M012_genomic.mergedLanes.bcftools.calls.vcf.gz .
cp ../FASTQ_public/BWA_public/PRJEB5837_MZM0403_M007_genomic.mergedLanes.bcftools.calls.vcf.gz .
cp ../FASTQ_public/BWA_public/PRJEB5837_MZM0403_M006_genomic.mergedLanes.bcftools.calls.vcf.gz .
cp ../FASTQ_public/BWA_public/PRJEB5837_GRZ_M005_genomic.mergedLanes.bcftools.calls.vcf.gz .
cp ../FASTQ_public/BWA_public/PRJEB5837_GRZ_M004_genomic.mergedLanes.bcftools.calls.vcf.gz .

# make index
#tabix -fp vcf  G1_CSFP200009002-1a_HNV2GDSXY_L3.srt.markDup.bcftools.calls_v2.vcf.gz
#tabix -fp vcf  G2_CSFP200009003-1a_HNV2GDSXY_L3.srt.markDup.bcftools.calls_v2.vcf.gz 
#tabix -fp vcf  Z1_CSFP200009004-1a_HNV2GDSXY_L3.srt.markDup.bcftools.calls_v2.vcf.gz 
#tabix -fp vcf  Z2_CSFP200009005-1a_HNV2GDSXY_L3.srt.markDup.bcftools.calls_v2.vcf.gz 
#tabix -fp vcf  SRR1261482_MZM0703.bcftools.calls.vcf.gz 
#tabix -fp vcf  POOLED_SRR171875x_MZM0403.bcftools.calls.vcf.gz

tabix -fp vcf PRJNA221357_SRR1261480_GRZ_P0FA.bcftools.calls.vcf.gz 
tabix -fp vcf PRJNA221357_SRR1261476_GRZ_BGI170.bcftools.calls.vcf.gz 
tabix -fp vcf PRJNA221357_SRR1261474_GRZ_GRZ340.bcftools.calls.vcf.gz 
tabix -fp vcf PRJNA221357_SRR1261472_GRZ_GRZ300.bcftools.calls.vcf.gz 
tabix -fp vcf PRJEB5837_MZZW0701_M002_genomic.mergedLanes.bcftools.calls.vcf.gz 
tabix -fp vcf PRJEB5837_MZZW0701_M001_genomic.mergedLanes.bcftools.calls.vcf.gz 
tabix -fp vcf PRJEB5837_MZM0410_M013_genomic.mergedLanes.bcftools.calls.vcf.gz 
tabix -fp vcf PRJEB5837_MZM0410_M012_genomic.mergedLanes.bcftools.calls.vcf.gz 
tabix -fp vcf PRJEB5837_MZM0403_M007_genomic.mergedLanes.bcftools.calls.vcf.gz 
tabix -fp vcf PRJEB5837_MZM0403_M006_genomic.mergedLanes.bcftools.calls.vcf.gz 
tabix -fp vcf PRJEB5837_GRZ_M005_genomic.mergedLanes.bcftools.calls.vcf.gz 
tabix -fp vcf PRJEB5837_GRZ_M004_genomic.mergedLanes.bcftools.calls.vcf.gz

####### combine all calls into one multi sample vcf
# https://www.biostars.org/p/469300/
bcftools merge --missing-to-ref -o 2021-09-28_Killifish_Strains_genotypes.merged.vcf G1_CSFP200009002-1a_HNV2GDSXY_L3.srt.markDup.bcftools.calls_v2.vcf.gz G2_CSFP200009003-1a_HNV2GDSXY_L3.srt.markDup.bcftools.calls_v2.vcf.gz Z1_CSFP200009004-1a_HNV2GDSXY_L3.srt.markDup.bcftools.calls_v2.vcf.gz Z2_CSFP200009005-1a_HNV2GDSXY_L3.srt.markDup.bcftools.calls_v2.vcf.gz SRR1261482_MZM0703.bcftools.calls.vcf.gz POOLED_SRR171875x_MZM0403.bcftools.calls.vcf.gz PRJNA221357_SRR1261480_GRZ_P0FA.bcftools.calls.vcf.gz PRJNA221357_SRR1261476_GRZ_BGI170.bcftools.calls.vcf.gz  PRJNA221357_SRR1261474_GRZ_GRZ340.bcftools.calls.vcf.gz  PRJNA221357_SRR1261472_GRZ_GRZ300.bcftools.calls.vcf.gz  PRJEB5837_MZZW0701_M002_genomic.mergedLanes.bcftools.calls.vcf.gz  PRJEB5837_MZZW0701_M001_genomic.mergedLanes.bcftools.calls.vcf.gz  PRJEB5837_MZM0410_M013_genomic.mergedLanes.bcftools.calls.vcf.gz  PRJEB5837_MZM0410_M012_genomic.mergedLanes.bcftools.calls.vcf.gz  PRJEB5837_MZM0403_M007_genomic.mergedLanes.bcftools.calls.vcf.gz  PRJEB5837_MZM0403_M006_genomic.mergedLanes.bcftools.calls.vcf.gz  PRJEB5837_GRZ_M005_genomic.mergedLanes.bcftools.calls.vcf.gz  PRJEB5837_GRZ_M004_genomic.mergedLanes.bcftools.calls.vcf.gz 

gzip 2021-09-28_Killifish_Strains_genotypes.merged.vcf


#prune the data
plink --vcf  2021-09-28_Killifish_Strains_genotypes.merged.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:#_\$1_\$2 --indep-pairwise 100 50 0.9 --out  2021-09-28_Killifish_Strains_genotypes.PRUNED

#Calculate PCs from genotypes:
plink --vcf  2021-09-28_Killifish_Strains_genotypes.merged.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:#_\$1_\$2 --extract 2021-09-28_Killifish_Strains_genotypes.PRUNED.prune.in --pca var-wts --out 2021-09-28_Killifish_Strains_genotypes_PRUNED_PLINK_pca

##### https://www.biostars.org/p/44735/