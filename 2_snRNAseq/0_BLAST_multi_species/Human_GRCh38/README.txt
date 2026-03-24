# 2024-03-11

# get protein sequences from NCBI for killifish 2015 assembly
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/465/895/GCF_001465895.1_Nfu_20140520/

# get Human proteins from Ensembl Biomart 111
#
#		Dataset 23820 / 70711 Genes
#		Human genes (GRCh38.p14)
#
#		Filters
#			Gene type: IG_C_gene , IG_D_gene , IG_J_gene , IG_V_gene , protein_coding , TR_C_gene , TR_D_gene , TR_J_gene , TR_V_gene
#		
#		Attributes
#			Peptide
#			Gene name
#			Gene stable ID
#			Protein stable ID
#			Transcript stable ID
#		
#				2024-03-11_Human_Proteome_Ens111.fa
#

# Get human proteome from Uniprot
     https://www.uniprot.org/proteomes/UP000005640
     only swissprot reviewed sequences [Download only reviewed (Swiss-Prot) canonical proteins (20,418)]

		uniprotkb_proteome_UP000005640_AND_revi_2024_03_12.fasta
		
		
	1 protein sequence per gene
	UP000005640_9606.fasta

Get uniprot to gene symbol conversion table from biomart
	
BLASTn tabular output format 6
Column headers: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

# reverse blast
# makeblastdb -in GCF_001465895.1_Nfu_20140520_protein.faa -out Nfu_20140520_protein -dbtype prot -title "Killi 2015 Database"
# blastp -query 2024-03-11_Human_Proteome_Ens111.fa -db Nfu_20140520_protein -evalue 1e-3 -out 2024-03-11_Human_Nfur2015_BestHits_1e-3_REVERSE.txt -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -outfmt 6 -max_target_seqs 1
blastp -query UP000005640_9606.fasta -db Nfu_20140520_protein -evalue 1e-3 -out 2024-03-11_HumanUNIPROT_Nfur2015_BestHits_1e-3_REVERSE.txt -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -outfmt 6 -max_target_seqs 1

blastp -query gencode.v45.pc_translations.fa -db Nfu_20140520_protein -evalue 1e-3 -out 2024-03-11_HumanGencode_Nfur2015_BestHits_1e-3_REVERSE.txt -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -outfmt 6 -max_target_seqs 1
