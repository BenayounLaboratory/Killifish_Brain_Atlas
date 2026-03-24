# # 2024-04-09
# # footprint analysis in DA peaks only
# # get clean peaks for footprint analysis
#  cat ../DESeq2/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age.bed \
#      ../DESeq2/2024-04-09_GRZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age.bed  \
#      ../DESeq2/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_DWN_with_Age.bed \
#      ../DESeq2/2024-04-09_ZMZ_brain_Aging_ATAC_DESeq2_Analysis_AGING_FDR5_UP_with_Age.bed  \
#      | perl -lane 'print if ($F[1]>0)' |sortBed -i - | mergeBed -i -  | cut -f 1,2,3,4 > POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed
# ###    26154 POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed
# 
# # # get pooled by condition bams
# ln -s ../NucleoATAC/Killi_Brain_ATAC_GRZ_OF_Merged.noMT.srt.bam
# ln -s ../NucleoATAC/Killi_Brain_ATAC_GRZ_OF_Merged.noMT.srt.bam.bai
# ln -s ../NucleoATAC/Killi_Brain_ATAC_GRZ_OM_Merged.noMT.srt.bam
# ln -s ../NucleoATAC/Killi_Brain_ATAC_GRZ_OM_Merged.noMT.srt.bam.bai
# ln -s ../NucleoATAC/Killi_Brain_ATAC_GRZ_YF_Merged.noMT.srt.bam
# ln -s ../NucleoATAC/Killi_Brain_ATAC_GRZ_YF_Merged.noMT.srt.bam.bai
# ln -s ../NucleoATAC/Killi_Brain_ATAC_GRZ_YM_Merged.noMT.srt.bam
# ln -s ../NucleoATAC/Killi_Brain_ATAC_GRZ_YM_Merged.noMT.srt.bam.bai
# ln -s ../NucleoATAC/Killi_Brain_ATAC_ZMZ_GF_Merged.noMT.srt.bam
# ln -s ../NucleoATAC/Killi_Brain_ATAC_ZMZ_GF_Merged.noMT.srt.bam.bai
# ln -s ../NucleoATAC/Killi_Brain_ATAC_ZMZ_GM_Merged.noMT.srt.bam
# ln -s ../NucleoATAC/Killi_Brain_ATAC_ZMZ_GM_Merged.noMT.srt.bam.bai
# ln -s ../NucleoATAC/Killi_Brain_ATAC_ZMZ_OF_Merged.noMT.srt.bam
# ln -s ../NucleoATAC/Killi_Brain_ATAC_ZMZ_OF_Merged.noMT.srt.bam.bai
# ln -s ../NucleoATAC/Killi_Brain_ATAC_ZMZ_OM_Merged.noMT.srt.bam
# ln -s ../NucleoATAC/Killi_Brain_ATAC_ZMZ_OM_Merged.noMT.srt.bam.bai
# ln -s ../NucleoATAC/Killi_Brain_ATAC_ZMZ_YF_Merged.noMT.srt.bam
# ln -s ../NucleoATAC/Killi_Brain_ATAC_ZMZ_YF_Merged.noMT.srt.bam.bai
# ln -s ../NucleoATAC/Killi_Brain_ATAC_ZMZ_YM_Merged.noMT.srt.bam
# ln -s ../NucleoATAC/Killi_Brain_ATAC_ZMZ_YM_Merged.noMT.srt.bam.bai
# 
# 
# # correct Tn5 bias
# TOBIAS ATACorrect --bam Killi_Brain_ATAC_GRZ_YF_Merged.noMT.srt.bam --prefix Killi_Brain_ATAC_GRZ_YF_ATACorrect_MACS2_DApeaks --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed  --cores 1
# TOBIAS ATACorrect --bam Killi_Brain_ATAC_GRZ_OF_Merged.noMT.srt.bam --prefix Killi_Brain_ATAC_GRZ_OF_ATACorrect_MACS2_DApeaks --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed  --cores 1
# TOBIAS ATACorrect --bam Killi_Brain_ATAC_GRZ_YM_Merged.noMT.srt.bam --prefix Killi_Brain_ATAC_GRZ_YM_ATACorrect_MACS2_DApeaks --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed  --cores 1
# TOBIAS ATACorrect --bam Killi_Brain_ATAC_GRZ_OM_Merged.noMT.srt.bam --prefix Killi_Brain_ATAC_GRZ_OM_ATACorrect_MACS2_DApeaks --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed  --cores 1
# TOBIAS ATACorrect --bam Killi_Brain_ATAC_ZMZ_YF_Merged.noMT.srt.bam --prefix Killi_Brain_ATAC_ZMZ_YF_ATACorrect_MACS2_DApeaks --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed  --cores 1
# TOBIAS ATACorrect --bam Killi_Brain_ATAC_ZMZ_OF_Merged.noMT.srt.bam --prefix Killi_Brain_ATAC_ZMZ_OF_ATACorrect_MACS2_DApeaks --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed  --cores 1
# TOBIAS ATACorrect --bam Killi_Brain_ATAC_ZMZ_GF_Merged.noMT.srt.bam --prefix Killi_Brain_ATAC_ZMZ_GF_ATACorrect_MACS2_DApeaks --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed  --cores 1
# TOBIAS ATACorrect --bam Killi_Brain_ATAC_ZMZ_YM_Merged.noMT.srt.bam --prefix Killi_Brain_ATAC_ZMZ_YM_ATACorrect_MACS2_DApeaks --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed  --cores 1
# TOBIAS ATACorrect --bam Killi_Brain_ATAC_ZMZ_OM_Merged.noMT.srt.bam --prefix Killi_Brain_ATAC_ZMZ_OM_ATACorrect_MACS2_DApeaks --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed  --cores 1
# TOBIAS ATACorrect --bam Killi_Brain_ATAC_ZMZ_GM_Merged.noMT.srt.bam --prefix Killi_Brain_ATAC_ZMZ_GM_ATACorrect_MACS2_DApeaks --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed  --cores 1
# 
# # score footprints
# # The main task in footprinting is to identify regions of protein binding across the genome. 
# # Using single basepair cutsite tracks (as produced by ATACorrect), TOBIAS ScoreBigwig is used to calculate a continuous footprinting score across regions. 
# TOBIAS FootprintScores --signal Killi_Brain_ATAC_GRZ_YF_ATACorrect_MACS2_DApeaks_corrected.bw --output Killi_Brain_ATAC_GRZ_YF_ATACorrect_MACS2_DApeaks_footprints.bw --regions POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed --cores 1 --score footprint
# TOBIAS FootprintScores --signal Killi_Brain_ATAC_GRZ_OF_ATACorrect_MACS2_DApeaks_corrected.bw --output Killi_Brain_ATAC_GRZ_OF_ATACorrect_MACS2_DApeaks_footprints.bw --regions POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed --cores 1 --score footprint
# TOBIAS FootprintScores --signal Killi_Brain_ATAC_GRZ_YM_ATACorrect_MACS2_DApeaks_corrected.bw --output Killi_Brain_ATAC_GRZ_YM_ATACorrect_MACS2_DApeaks_footprints.bw --regions POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed --cores 1 --score footprint
# TOBIAS FootprintScores --signal Killi_Brain_ATAC_GRZ_OM_ATACorrect_MACS2_DApeaks_corrected.bw --output Killi_Brain_ATAC_GRZ_OM_ATACorrect_MACS2_DApeaks_footprints.bw --regions POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed --cores 1 --score footprint
# TOBIAS FootprintScores --signal Killi_Brain_ATAC_ZMZ_YF_ATACorrect_MACS2_DApeaks_corrected.bw --output Killi_Brain_ATAC_ZMZ_YF_ATACorrect_MACS2_DApeaks_footprints.bw --regions POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed --cores 1 --score footprint
# TOBIAS FootprintScores --signal Killi_Brain_ATAC_ZMZ_OF_ATACorrect_MACS2_DApeaks_corrected.bw --output Killi_Brain_ATAC_ZMZ_OF_ATACorrect_MACS2_DApeaks_footprints.bw --regions POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed --cores 1 --score footprint
# TOBIAS FootprintScores --signal Killi_Brain_ATAC_ZMZ_GF_ATACorrect_MACS2_DApeaks_corrected.bw --output Killi_Brain_ATAC_ZMZ_GF_ATACorrect_MACS2_DApeaks_footprints.bw --regions POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed --cores 1 --score footprint
# TOBIAS FootprintScores --signal Killi_Brain_ATAC_ZMZ_YM_ATACorrect_MACS2_DApeaks_corrected.bw --output Killi_Brain_ATAC_ZMZ_YM_ATACorrect_MACS2_DApeaks_footprints.bw --regions POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed --cores 1 --score footprint
# TOBIAS FootprintScores --signal Killi_Brain_ATAC_ZMZ_OM_ATACorrect_MACS2_DApeaks_corrected.bw --output Killi_Brain_ATAC_ZMZ_OM_ATACorrect_MACS2_DApeaks_footprints.bw --regions POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed --cores 1 --score footprint
# TOBIAS FootprintScores --signal Killi_Brain_ATAC_ZMZ_GM_ATACorrect_MACS2_DApeaks_corrected.bw --output Killi_Brain_ATAC_ZMZ_GM_ATACorrect_MACS2_DApeaks_footprints.bw --regions POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed --cores 1 --score footprint
# 
#removeOutOfBoundsReads.pl POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.bed Killi2015 > POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed
###ulimit -n 1000000
# detect differential binding at JASPAR motifs across conditions
TOBIAS BINDetect --signals Killi_Brain_ATAC_GRZ_YF_ATACorrect_MACS2_DApeaks_footprints.bw Killi_Brain_ATAC_GRZ_OF_ATACorrect_MACS2_DApeaks_footprints.bw --outdir GRZ_YF_OF_BINDetect_JASPAR --cond_names YF OF --motifs JASPAR2022_CORE_non-redundant_pfms_jaspar.txt  --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed --cores 1
TOBIAS BINDetect --signals Killi_Brain_ATAC_GRZ_YM_ATACorrect_MACS2_DApeaks_footprints.bw Killi_Brain_ATAC_GRZ_OM_ATACorrect_MACS2_DApeaks_footprints.bw --outdir GRZ_YM_OM_BINDetect_JASPAR --cond_names YM OM --motifs JASPAR2022_CORE_non-redundant_pfms_jaspar.txt  --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed --cores 1
TOBIAS BINDetect --signals Killi_Brain_ATAC_ZMZ_YF_ATACorrect_MACS2_DApeaks_footprints.bw Killi_Brain_ATAC_ZMZ_OF_ATACorrect_MACS2_DApeaks_footprints.bw --outdir ZMZ_YF_OF_BINDetect_JASPAR --cond_names YF OF --motifs JASPAR2022_CORE_non-redundant_pfms_jaspar.txt  --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed --cores 1
TOBIAS BINDetect --signals Killi_Brain_ATAC_ZMZ_YF_ATACorrect_MACS2_DApeaks_footprints.bw Killi_Brain_ATAC_ZMZ_GF_ATACorrect_MACS2_DApeaks_footprints.bw --outdir ZMZ_YF_GF_BINDetect_JASPAR --cond_names YF GF --motifs JASPAR2022_CORE_non-redundant_pfms_jaspar.txt  --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed --cores 1
TOBIAS BINDetect --signals Killi_Brain_ATAC_ZMZ_YM_ATACorrect_MACS2_DApeaks_footprints.bw Killi_Brain_ATAC_ZMZ_OM_ATACorrect_MACS2_DApeaks_footprints.bw --outdir ZMZ_YM_OM_BINDetect_JASPAR --cond_names YM OM --motifs JASPAR2022_CORE_non-redundant_pfms_jaspar.txt  --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed --cores 1
TOBIAS BINDetect --signals Killi_Brain_ATAC_ZMZ_YM_ATACorrect_MACS2_DApeaks_footprints.bw Killi_Brain_ATAC_ZMZ_GM_ATACorrect_MACS2_DApeaks_footprints.bw --outdir ZMZ_YM_OG_BINDetect_JASPAR --cond_names YM GO --motifs JASPAR2022_CORE_non-redundant_pfms_jaspar.txt  --genome /Volumes/BB_HQ_3/Brain_aging_Simons/GENOME_REF/2021-12-02_GCF_001465895.1_MASKED_withFishTEDB_NR.fa --peaks POOLED_KilliBrainAging_ATACseq_MACS2_peaks_Diffbind_FDR5peaks_ONLY.FIX.bed --cores 1
