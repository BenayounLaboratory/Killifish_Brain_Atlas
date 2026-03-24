cp /Volumes/BB_HQ_3/Brain_aging_Simons/ATAC/Bulk_ATAC_Analysis_Folder/MACS2/*narrowPeak ./Brain_ATAC/

mspc -f ./Brain_ATAC/ -c 5  -r bio -w 1e-5 -s 1e-9 -m Lowest --excludeHeader -o MSPC_Killifish_Brain_ATAC_05minimum  # out of 47 samples
mspc -f ./Brain_ATAC/ -c 10 -r bio -w 1e-5 -s 1e-9 -m Lowest --excludeHeader -o MSPC_Killifish_Brain_ATAC_10minimum  # out of 47 samples
mspc -f ./Brain_ATAC/ -c 15 -r bio -w 1e-5 -s 1e-9 -m Lowest --excludeHeader -o MSPC_Killifish_Brain_ATAC_15minimum  # out of 47 samples
mspc -f ./Brain_ATAC/ -c 20 -r bio -w 1e-5 -s 1e-9 -m Lowest --excludeHeader -o MSPC_Killifish_Brain_ATAC_20minimum  # out of 47 samples


cp ./MSPC_Killifish_Brain_ATAC_05minimum/ConsensusPeaks.bed  Killifish_Brain_ATAC_05minimum_ATAC_MSPC_ConsensusPeaks.bed  
cp ./MSPC_Killifish_Brain_ATAC_10minimum/ConsensusPeaks.bed  Killifish_Brain_ATAC_10minimum_ATAC_MSPC_ConsensusPeaks.bed  
cp ./MSPC_Killifish_Brain_ATAC_15minimum/ConsensusPeaks.bed  Killifish_Brain_ATAC_15minimum_ATAC_MSPC_ConsensusPeaks.bed  
cp ./MSPC_Killifish_Brain_ATAC_20minimum/ConsensusPeaks.bed  Killifish_Brain_ATAC_20minimum_ATAC_MSPC_ConsensusPeaks.bed  

wc -l *_ConsensusPeaks.bed  
#  105053 Killifish_Brain_ATAC_05minimum_ATAC_MSPC_ConsensusPeaks.bed
#   89980 Killifish_Brain_ATAC_10minimum_ATAC_MSPC_ConsensusPeaks.bed
#   80727 Killifish_Brain_ATAC_15minimum_ATAC_MSPC_ConsensusPeaks.bed
#   73472 Killifish_Brain_ATAC_20minimum_ATAC_MSPC_ConsensusPeaks.bed
#  349232 total
  
annotatePeaks.pl Killifish_Brain_ATAC_05minimum_ATAC_MSPC_ConsensusPeaks.bed   Killi2015g > HOMER_Killifish_Brain_ATAC_05minimum_ATAC_MSPC_ConsensusPeaks.xls  
annotatePeaks.pl Killifish_Brain_ATAC_10minimum_ATAC_MSPC_ConsensusPeaks.bed   Killi2015g > HOMER_Killifish_Brain_ATAC_10minimum_ATAC_MSPC_ConsensusPeaks.xls  
annotatePeaks.pl Killifish_Brain_ATAC_15minimum_ATAC_MSPC_ConsensusPeaks.bed   Killi2015g > HOMER_Killifish_Brain_ATAC_15minimum_ATAC_MSPC_ConsensusPeaks.xls  
annotatePeaks.pl Killifish_Brain_ATAC_20minimum_ATAC_MSPC_ConsensusPeaks.bed   Killi2015g > HOMER_Killifish_Brain_ATAC_20minimum_ATAC_MSPC_ConsensusPeaks.xls  

### Will run diffbind with peaks from 20 or more samples