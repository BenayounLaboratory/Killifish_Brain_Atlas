# 2024-05-03
# footprint analysis in bulk DA peaks only
# by cell type, pseudobulked, F and M merged

# get clean peaks for footprint analysis
# ln -s /Volumes/BB_HQ_3/Brain_aging_Simons/ATAC/Bulk_ATAC_Analysis_Folder/TOBIAS/POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed
###    26154 POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed


# run a loop over files
for f in $(find . -name '*scATAC_PB.srt.bam')
do
    of=$(basename "${f}" | sed 's/\.srt\.bam/ATACorrect_MACS2_DApeaks/g');
    oFname="./${of}"
    
    # correct Tn5 bias
	TOBIAS ATACorrect --bam $f --prefix $oFname --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed  --cores 1

    of2=$(basename "${f}" | sed 's/\.srt\.bam/ATACorrect_MACS2_DApeaks_corrected\.bw/g');
    oFname2="./${of2}"
    
    of3=$(basename "${f}" | sed 's/\.srt\.bam/ATACorrect_MACS2_DApeaks_footprints\.bw/g');
    oFname3="./${of3}"
    
	# score footprints
	TOBIAS FootprintScores --signal $oFname2 --output $oFname3 --regions POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed --cores 1 --score footprint

done
ulimit -n 1000000

# detect differential binding at JASPAR motifs across conditions
TOBIAS BINDetect --signals Astrocytes_Radial_Glia_YF_YM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw Astrocytes_Radial_Glia_OF_OM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw         --outdir Astrocytes_Radial_Glia_Y_vs_O_scATAC_BINDetect_JASPAR      --cond_names Y O --motifs JASPAR2022_CORE_non-redundant_pfms_jaspar.txt  --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed --cores 1
TOBIAS BINDetect --signals Ependymal_cells_YF_YM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw Ependymal_cells_OF_OM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw                       --outdir Ependymal_cells_Y_vs_O_scATAC_BINDetect_JASPAR             --cond_names Y O --motifs JASPAR2022_CORE_non-redundant_pfms_jaspar.txt  --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed --cores 1
TOBIAS BINDetect --signals GABAergic_neurons_YF_YM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw GABAergic_neurons_OF_OM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw                   --outdir GABAergic_neurons_Y_vs_O_scATAC_BINDetect_JASPAR           --cond_names Y O --motifs JASPAR2022_CORE_non-redundant_pfms_jaspar.txt  --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed --cores 1
TOBIAS BINDetect --signals Granule_Excitatory_Neurons_YF_YM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw Granule_Excitatory_Neurons_OF_OM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw --outdir Granule_Excitatory_Neurons_Y_vs_O_scATAC_BINDetect_JASPAR  --cond_names Y O --motifs JASPAR2022_CORE_non-redundant_pfms_jaspar.txt  --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed --cores 1
TOBIAS BINDetect --signals Microglia_YF_YM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw Microglia_OF_OM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw                                   --outdir Microglia_Y_vs_O_scATAC_BINDetect_JASPAR                   --cond_names Y O --motifs JASPAR2022_CORE_non-redundant_pfms_jaspar.txt  --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed --cores 1
TOBIAS BINDetect --signals NSPCs_YF_YM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw NSPCs_OF_OM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw                                           --outdir NSPCs_Y_vs_O_scATAC_BINDetect_JASPAR                       --cond_names Y O --motifs JASPAR2022_CORE_non-redundant_pfms_jaspar.txt  --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed --cores 1
TOBIAS BINDetect --signals Oligodendrocytes_YF_YM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw Oligodendrocytes_OF_OM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw                     --outdir Oligodendrocytes_Y_vs_O_scATAC_BINDetect_JASPAR            --cond_names Y O --motifs JASPAR2022_CORE_non-redundant_pfms_jaspar.txt  --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed --cores 1
TOBIAS BINDetect --signals OPCs_YF_YM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw OPCs_OF_OM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw                                             --outdir OPCs_Y_vs_O_scATAC_BINDetect_JASPAR                        --cond_names Y O --motifs JASPAR2022_CORE_non-redundant_pfms_jaspar.txt  --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed --cores 1
TOBIAS BINDetect --signals PV_interneurons_YF_YM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw PV_interneurons_OF_OM_scATAC_PBATACorrect_MACS2_DApeaks_footprints.bw                       --outdir PV_interneurons_Y_vs_O_scATAC_BINDetect_JASPAR             --cond_names Y O --motifs JASPAR2022_CORE_non-redundant_pfms_jaspar.txt  --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed --cores 1
