#!/bin/bash
#########################
# ChIPmentation project #
#########################

SAMPLES_FILE=/home/arendeiro/projects/chipmentation/samples.txt
PEAKS_FILE=/home/arendeiro/projects/chipmentation/samples_peaks.txt


### PREPROCESSING
# Trim, map, mark duplicates, shift reads, qc, create homer dir
while read SAMPLE_NAME SAMPLE_FILE; do
    sbatch /home/arendeiro/projects/chipmentation/src/pipelines/mapping_pipeline.sh $SAMPLE_NAME $SAMPLE_FILE
done < $SAMPLES_FILE

# Merge technical replicates
sh scripts/merge_replicates.sh # <- for two technical replicates
sh scripts/joining_replicates.sh # <- by hand for more than two replicates


### CORRELATIONS
while read SAMPLE_NAME CONTROL_NAME; do
    sbatch /home/arendeiro/projects/chipmentation/src/scripts/counts_on_windows.sh $SAMPLE_NAME
    sbatch /home/arendeiro/projects/chipmentation/src/scripts/counts_on_windows.sh $CONTROL_NAME
done < $PEAKS_FILE
# concatenate files
sh scripts/concatenate_counts_on_windows.sh
# plot correlations
R scripts/plot_correlations.R # needs work to accept arguments


### PEAKS
# Call peaks with MACS2 and homer
while read SAMPLE_NAME CONTROL_FILE; do
sbatch /home/arendeiro/projects/chipmentation/src/scripts/mapping_pipeline.sh $SAMPLE_NAME $CONTROL_FILE
done < $PEAKS_FILE

# Derive consensus peaks between replicates, count #'s
R scripts/diffBind_analysis.R
sh scripts/peak_intersection.sh #<- by hand, very bad

# Find motifs in TF samples
sh # missing!!!


### FOOTPRINTING

## H3K4me3 on TSSs
# Get CAGE data
# Filter for high-scoring CAGE peaks
# Get window around CAGE summits with orientation
# Count reads, GC content and conservation on windows around CAGE summits
# Plot
# Plot xlim(-100,100)
# Filter on TATA-containing promoters

# Count reads on windows centered on peaks' summits
#run_counts_on_windows_peaks.sh # deprecated

# Center peaks on motif, get 4kb window around 
# Count reads around motifs
sh scripts/footprint_analysis.sh
R scripts/plot_footprints.R

# Run CENTIPEDE

# Run pyDNASE


# Autocorrelations in H3K4me3 peaks
python scripts/distance_reads_correlations.py