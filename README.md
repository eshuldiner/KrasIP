# KrasIP Tuba-seq

This repository contains code used to analyze the tumor barcode sequencing ("Tuba-seq") data from the manuscript *Multiplexed identification of RAS paralog imbalance as a driver of lung cancer growth*.

## Contents

The following scripts analyze the raw data (fastq files produced from sequencing of the sgID-barcode amplicon) sequentially to produce estimates of the size of each tumor in each mouse:
1. **count_reads.py**: counts reads that map to each unique sgID-barcode combination in each sample.
2. **filter_tumors.py**: removes spurious tumors and likely contaminants.
3. **convert_to_cells.py**: converts read counts to estimates of tumor size (cell counts)

The script **process_tumors.py** collects the processed data across samples to produce a single dataset that is the starting point for most statistical analyses and visualizations.

The script **statistical_analysis.py** contains the code used to perform statistical analyses of the effects of different tumor genotypes on tumor fitness (i.e. relative tumor size, relative tumor number and relative tumor burden).
