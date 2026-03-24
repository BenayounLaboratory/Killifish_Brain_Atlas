export CRDIR1="/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2022-04-15_Killi_Brain_Aging_scATACseq_cohort1"
export CRDIR2="/Volumes/BB_HQ_3/Brain_aging_Simons/scATACseq/2023-01-05_Killi_Brain_Aging_scATACseq_cohort2"


### 
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Astrocytes_Radial_Glia_OF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Astrocytes_Radial_Glia_OF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Astrocytes_Radial_Glia_OF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Astrocytes_Radial_Glia_OF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Astrocytes_Radial_Glia_OM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Astrocytes_Radial_Glia_OM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Astrocytes_Radial_Glia_OM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Astrocytes_Radial_Glia_OM2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Astrocytes_Radial_Glia_YF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Astrocytes_Radial_Glia_YF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Astrocytes_Radial_Glia_YF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Astrocytes_Radial_Glia_YF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Astrocytes_Radial_Glia_YM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Astrocytes_Radial_Glia_YM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Astrocytes_Radial_Glia_YM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Astrocytes_Radial_Glia_YM2.bam

### 
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Ependymal_cells_OF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Ependymal_cells_OF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Ependymal_cells_OF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Ependymal_cells_OF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Ependymal_cells_OM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Ependymal_cells_OM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Ependymal_cells_OM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Ependymal_cells_OM2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Ependymal_cells_YF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Ependymal_cells_YF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Ependymal_cells_YF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Ependymal_cells_YF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Ependymal_cells_YM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Ependymal_cells_YM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Ependymal_cells_YM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Ependymal_cells_YM2.bam

### 
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_GABAergic_neurons_OF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_GABAergic_neurons_OF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_GABAergic_neurons_OF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_GABAergic_neurons_OF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_GABAergic_neurons_OM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_GABAergic_neurons_OM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_GABAergic_neurons_OM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_GABAergic_neurons_OM2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_GABAergic_neurons_YF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_GABAergic_neurons_YF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_GABAergic_neurons_YF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_GABAergic_neurons_YF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_GABAergic_neurons_YM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_GABAergic_neurons_YM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_GABAergic_neurons_YM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_GABAergic_neurons_YM2.bam

### 
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Granule_Excitatory_Neurons_OF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Granule_Excitatory_Neurons_OF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Granule_Excitatory_Neurons_OF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Granule_Excitatory_Neurons_OF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Granule_Excitatory_Neurons_OM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Granule_Excitatory_Neurons_OM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Granule_Excitatory_Neurons_OM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Granule_Excitatory_Neurons_OM2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Granule_Excitatory_Neurons_YF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Granule_Excitatory_Neurons_YF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Granule_Excitatory_Neurons_YF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Granule_Excitatory_Neurons_YF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Granule_Excitatory_Neurons_YM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Granule_Excitatory_Neurons_YM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Granule_Excitatory_Neurons_YM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Granule_Excitatory_Neurons_YM2.bam

### 
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Microglia_OF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Microglia_OF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Microglia_OF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Microglia_OF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Microglia_OM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Microglia_OM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Microglia_OM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Microglia_OM2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Microglia_YF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Microglia_YF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Microglia_YF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Microglia_YF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Microglia_YM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Microglia_YM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Microglia_YM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Microglia_YM2.bam

### 
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_NSPCs_OF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_NSPCs_OF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_NSPCs_OF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_NSPCs_OF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_NSPCs_OM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_NSPCs_OM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_NSPCs_OM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_NSPCs_OM2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_NSPCs_YF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_NSPCs_YF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_NSPCs_YF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_NSPCs_YF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_NSPCs_YM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_NSPCs_YM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_NSPCs_YM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_NSPCs_YM2.bam

### 
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Oligodendrocytes_OF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Oligodendrocytes_OF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Oligodendrocytes_OF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Oligodendrocytes_OF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Oligodendrocytes_OM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Oligodendrocytes_OM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Oligodendrocytes_OM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Oligodendrocytes_OM2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Oligodendrocytes_YF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Oligodendrocytes_YF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Oligodendrocytes_YF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Oligodendrocytes_YF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Oligodendrocytes_YM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_Oligodendrocytes_YM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_Oligodendrocytes_YM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_Oligodendrocytes_YM2.bam

### 
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_OPCs_OF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_OPCs_OF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_OPCs_OF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_OPCs_OF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_OPCs_OM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_OPCs_OM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_OPCs_OM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_OPCs_OM2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_OPCs_YF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_OPCs_YF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_OPCs_YF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_OPCs_YF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_OPCs_YM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_OPCs_YM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_OPCs_YM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_OPCs_YM2.bam

### 
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_PV_interneurons_OF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_PV_interneurons_OF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_PV_interneurons_OF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_PV_interneurons_OF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_OM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_PV_interneurons_OM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_PV_interneurons_OM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_OM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_PV_interneurons_OM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_PV_interneurons_OM2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YF1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_PV_interneurons_YF1_barcode_list.txt --cores 1 --out-bam 2024-04-24_PV_interneurons_YF1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YF2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_PV_interneurons_YF2_barcode_list.txt --cores 1 --out-bam 2024-04-24_PV_interneurons_YF2.bam
subset-bam --bam $CRDIR1/CellRanger/Killi_Brain_YM1_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_PV_interneurons_YM1_barcode_list.txt --cores 1 --out-bam 2024-04-24_PV_interneurons_YM1.bam
subset-bam --bam $CRDIR2/CellRanger/Killi_Brain_YM2_FishTEDB_NR/outs/possorted_bam.bam --cell-barcodes 2024-04-24_PV_interneurons_YM2_barcode_list.txt --cores 1 --out-bam 2024-04-24_PV_interneurons_YM2.bam

