# Get JASPAR core/vertebrate, 2022-04-20
# JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt

# filtering is differnt for gtf on ATAC
# make concatenated versions for cellranger genome generate step
cat GCF_001465895.1_Nfu_20140520_genomic_CLEAN_MT_exon.gtf 2021-09-06_FishTEDB_Nothobranchius_furzeri.cdhit90perc.gtf > 2022-04-20_GCF_001465895.1_genes_withFishTEDB_forATAC.gtf


~/Softwares/cellranger-atac-2.1.0/bin/cellranger-atac  mkref --config=~/Softwares/cellranger-atac-2.1.0/Reference/killi.genome.config