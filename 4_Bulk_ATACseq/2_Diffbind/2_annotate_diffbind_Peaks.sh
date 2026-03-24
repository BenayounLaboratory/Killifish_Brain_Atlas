annotatePeaks.pl 2024-04-08_GRZ_Aging_Brains_MACS2_MSPC_diffbind_peaks.bed   Killi2015g > HOMER_2024-04-08_GRZ_Aging_Brains_MACS2_MSPC_diffbind_peaks.xls  
annotatePeaks.pl 2024-04-08_ZMZ_Aging_Brains_MACS2_MSPC_diffbind_peaks.bed   Killi2015g > HOMER_2024-04-08_ZMZ_Aging_Brains_MACS2_MSPC_diffbind_peaks.xls  


### intersect with H3K4me3 peaks
# remove negative bed entries
cat 2024-04-08_GRZ_Aging_Brains_MACS2_MSPC_diffbind_peaks.bed | perl -lane 'if ($F[1] < 0) { print "$F[0]\t0\t$F[2]";} else {print;}' > 2024-04-08_GRZ_Aging_Brains_MACS2_MSPC_diffbind_peaks.CLEAN.bed
cat 2024-04-08_ZMZ_Aging_Brains_MACS2_MSPC_diffbind_peaks.bed | perl -lane 'if ($F[1] < 0) { print "$F[0]\t0\t$F[2]";} else {print;}' > 2024-04-08_ZMZ_Aging_Brains_MACS2_MSPC_diffbind_peaks.CLEAN.bed


intersectBed -wao -a 2024-04-08_GRZ_Aging_Brains_MACS2_MSPC_diffbind_peaks.CLEAN.bed -b /Volumes/BB_HQ_3/Brain_aging_Simons/ATAC/Bulk_ATAC_Analysis_Folder/H3K4me3_ChIPseq/MACS2/SRR1556403_GRZ_YM_Brain_H3K4me3_ChIPseq_peaks.broadPeak > GRZ_Aging_Brains_MACS2_MSPC_diffbind_peaks_H3K4me3_overlap.bed
intersectBed -wao -a 2024-04-08_ZMZ_Aging_Brains_MACS2_MSPC_diffbind_peaks.CLEAN.bed -b /Volumes/BB_HQ_3/Brain_aging_Simons/ATAC/Bulk_ATAC_Analysis_Folder/H3K4me3_ChIPseq/MACS2/SRR1556403_GRZ_YM_Brain_H3K4me3_ChIPseq_peaks.broadPeak > ZMZ_Aging_Brains_MACS2_MSPC_diffbind_peaks_H3K4me3_overlap.bed