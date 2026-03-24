
# 2022-02-11

# get protein sequences from NCBI for killifish 2015 assembly
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/465/895/GCF_001465895.1_Nfu_20140520/

# get Zebrafish proteins from Ensembl Biomart 105

		Dataset 30368 / 37241 Genes
		Zebrafish genes (GRCz11)

		Filters
			Gene type: IG_C_gene , IG_C_pseudogene , IG_J_pseudogene , IG_pseudogene , IG_V_pseudogene , processed_pseudogene , protein_coding , pseudogene , TR_D_gene , TR_J_gene , TR_V_gene
		
		Attributes
			Peptide
			Gene name
			Gene stable ID
			Protein stable ID
			Transcript stable ID
		
				2022-02-11_Zabrafish_Proteome_Ens105.fa

# Zebrafish proteins from uniprot (download Zebrafish proteome)
	 https://www.uniprot.org/proteomes/?query=zebrafish&sort=score
	 UP000000437	Danio rerio (Zebrafish) (Brachydanio rerio) (Strain: Tuebingen)	7955	46844	 n:3640
C:95.5% (S:52.9% D:42.6%) F:1.7% M:2.8% actinopterygii_odb10	Close to standard (high value)	full
		UniProtKB (46,844)
		Swiss-Prot (3,216)
		TrEMBL (43,628)

			uniprot-proteome_UP000000437.fasta

makeblastdb -in 2022-02-11_Zebrafish_Proteome_Ens105.fa -out Zebrafish_Proteome_Ens105 -dbtype prot -title "GRCz11 Database"

blastp -query GCF_001465895.1_Nfu_20140520_protein.faa  -db Zebrafish_Proteome_Ens105 -evalue 1e-5 -out 2022-02-11_Nfur2015_Zebrafish_BestHits_1e-5.txt -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -outfmt 6 -max_target_seqs 1


BLASTn tabular output format 6
Column headers: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore


#### 2022-02-22
# reverse blast
makeblastdb -in GCF_001465895.1_Nfu_20140520_protein.faa -out Nfu_20140520_protein -dbtype prot -title "Killi 2015 Database"
blastp -query 2022-02-11_Zebrafish_Proteome_Ens105.fa  -db Nfu_20140520_protein -evalue 1e-5 -out 2022-02-22_Zebrafish_Nfur2015_BestHits_1e-5_REVERSE.txt -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -outfmt 6 -max_target_seqs 1
