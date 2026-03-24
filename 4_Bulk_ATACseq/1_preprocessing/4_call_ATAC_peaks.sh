# https://informatics.fas.harvard.edu/atac-seq-guidelines-old-version.html#peak

MAPPEDBAM="/Volumes/BB_HQ_3/Brain_aging_Simons/ATAC/Bulk_ATAC_Analysis_Folder/BAM_SORT"

/Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c1_Brain_ATACseq_OF1 -t $MAPPEDBAM/GRZ_c1_Brain_ATACseq_OF1.srt.bam -f "BAMPE" -g 1.25e9
/Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c1_Brain_ATACseq_OF2 -t $MAPPEDBAM/GRZ_c1_Brain_ATACseq_OF2.srt.bam -f "BAMPE" -g 1.25e9
/Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c1_Brain_ATACseq_OF3 -t $MAPPEDBAM/GRZ_c1_Brain_ATACseq_OF3.srt.bam -f "BAMPE" -g 1.25e9
/Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c1_Brain_ATACseq_OM1 -t $MAPPEDBAM/GRZ_c1_Brain_ATACseq_OM1.srt.bam -f "BAMPE" -g 1.25e9
/Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c1_Brain_ATACseq_OM2 -t $MAPPEDBAM/GRZ_c1_Brain_ATACseq_OM2.srt.bam -f "BAMPE" -g 1.25e9
/Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c1_Brain_ATACseq_OM3 -t $MAPPEDBAM/GRZ_c1_Brain_ATACseq_OM3.srt.bam -f "BAMPE" -g 1.25e9
/Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c1_Brain_ATACseq_YF1 -t $MAPPEDBAM/GRZ_c1_Brain_ATACseq_YF1.srt.bam -f "BAMPE" -g 1.25e9
/Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c1_Brain_ATACseq_YF2 -t $MAPPEDBAM/GRZ_c1_Brain_ATACseq_YF2.srt.bam -f "BAMPE" -g 1.25e9
/Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c1_Brain_ATACseq_YF3 -t $MAPPEDBAM/GRZ_c1_Brain_ATACseq_YF3.srt.bam -f "BAMPE" -g 1.25e9
/Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c1_Brain_ATACseq_YM1 -t $MAPPEDBAM/GRZ_c1_Brain_ATACseq_YM1.srt.bam -f "BAMPE" -g 1.25e9
/Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c1_Brain_ATACseq_YM2 -t $MAPPEDBAM/GRZ_c1_Brain_ATACseq_YM2.srt.bam -f "BAMPE" -g 1.25e9
/Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c1_Brain_ATACseq_YM3 -t $MAPPEDBAM/GRZ_c1_Brain_ATACseq_YM3.srt.bam -f "BAMPE" -g 1.25e9


# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c2_Brain_ATACseq_OF1 -t $MAPPEDBAM/GRZ_c2_Brain_ATACseq_OF1.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c2_Brain_ATACseq_OF2 -t $MAPPEDBAM/GRZ_c2_Brain_ATACseq_OF2.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c2_Brain_ATACseq_OF3 -t $MAPPEDBAM/GRZ_c2_Brain_ATACseq_OF3.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c2_Brain_ATACseq_OM1 -t $MAPPEDBAM/GRZ_c2_Brain_ATACseq_OM1.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c2_Brain_ATACseq_OM2 -t $MAPPEDBAM/GRZ_c2_Brain_ATACseq_OM2.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c2_Brain_ATACseq_OM3 -t $MAPPEDBAM/GRZ_c2_Brain_ATACseq_OM3.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c2_Brain_ATACseq_YF1 -t $MAPPEDBAM/GRZ_c2_Brain_ATACseq_YF1.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c2_Brain_ATACseq_YF2 -t $MAPPEDBAM/GRZ_c2_Brain_ATACseq_YF2.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c2_Brain_ATACseq_YF3 -t $MAPPEDBAM/GRZ_c2_Brain_ATACseq_YF3.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c2_Brain_ATACseq_YM1 -t $MAPPEDBAM/GRZ_c2_Brain_ATACseq_YM1.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c2_Brain_ATACseq_YM2 -t $MAPPEDBAM/GRZ_c2_Brain_ATACseq_YM2.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n GRZ_c2_Brain_ATACseq_YM3 -t $MAPPEDBAM/GRZ_c2_Brain_ATACseq_YM3.srt.bam -f "BAMPE" -g 1.25e9

# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c2_Brain_ATACseq_GF1 -t $MAPPEDBAM/ZMZ_c2_Brain_ATACseq_GF1.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c2_Brain_ATACseq_GF2 -t $MAPPEDBAM/ZMZ_c2_Brain_ATACseq_GF2.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c2_Brain_ATACseq_GM1 -t $MAPPEDBAM/ZMZ_c2_Brain_ATACseq_GM1.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c2_Brain_ATACseq_GM2 -t $MAPPEDBAM/ZMZ_c2_Brain_ATACseq_GM2.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c2_Brain_ATACseq_OF1 -t $MAPPEDBAM/ZMZ_c2_Brain_ATACseq_OF1.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c2_Brain_ATACseq_OM1 -t $MAPPEDBAM/ZMZ_c2_Brain_ATACseq_OM1.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c2_Brain_ATACseq_OM2 -t $MAPPEDBAM/ZMZ_c2_Brain_ATACseq_OM2.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c2_Brain_ATACseq_YF1 -t $MAPPEDBAM/ZMZ_c2_Brain_ATACseq_YF1.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c2_Brain_ATACseq_YF2 -t $MAPPEDBAM/ZMZ_c2_Brain_ATACseq_YF2.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c2_Brain_ATACseq_YM1 -t $MAPPEDBAM/ZMZ_c2_Brain_ATACseq_YM1.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c2_Brain_ATACseq_YM2 -t $MAPPEDBAM/ZMZ_c2_Brain_ATACseq_YM2.srt.bam -f "BAMPE" -g 1.25e9


# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c3_Brain_ATACseq_GF3 -t $MAPPEDBAM/ZMZ_c3_Brain_ATACseq_GF3.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c3_Brain_ATACseq_GF4 -t $MAPPEDBAM/ZMZ_c3_Brain_ATACseq_GF4.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c3_Brain_ATACseq_GM3 -t $MAPPEDBAM/ZMZ_c3_Brain_ATACseq_GM3.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c3_Brain_ATACseq_GM4 -t $MAPPEDBAM/ZMZ_c3_Brain_ATACseq_GM4.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c3_Brain_ATACseq_OM3 -t $MAPPEDBAM/ZMZ_c3_Brain_ATACseq_OM3.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c3_Brain_ATACseq_OM4 -t $MAPPEDBAM/ZMZ_c3_Brain_ATACseq_OM4.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c3_Brain_ATACseq_OF3 -t $MAPPEDBAM/ZMZ_c3_Brain_ATACseq_OF3.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c3_Brain_ATACseq_OF4 -t $MAPPEDBAM/ZMZ_c3_Brain_ATACseq_OF4.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c3_Brain_ATACseq_YF3 -t $MAPPEDBAM/ZMZ_c3_Brain_ATACseq_YF3.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c3_Brain_ATACseq_YF4 -t $MAPPEDBAM/ZMZ_c3_Brain_ATACseq_YF4.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c3_Brain_ATACseq_YM3 -t $MAPPEDBAM/ZMZ_c3_Brain_ATACseq_YM3.srt.bam -f "BAMPE" -g 1.25e9
# /Library/Frameworks/Python.framework/Versions/3.7/bin/macs2 callpeak -n ZMZ_c3_Brain_ATACseq_YM4 -t $MAPPEDBAM/ZMZ_c3_Brain_ATACseq_YM4.srt.bam -f "BAMPE" -g 1.25e9
