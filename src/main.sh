#!/bin/bash
#########################
# ChIPmentation project #
#########################

SAMPLES_FILE=/home/arendeiro/projects/chipmentation/samples.txt
PEAKS_FILE=/home/arendeiro/projects/chipmentation/samples_peaks.txt
REPLICATE_FILE=/home/arendeiro/projects/chipmentation/samples_to_merge.txt
PEAKS_FILE_MERGED=/home/arendeiro/projects/chipmentation/samples_peaks_merged.txt


### PREPROCESSING
# Trim, map, mark duplicates, shift reads, qc, create homer dir
while read SAMPLE_NAME SAMPLE_FILE; do
    sbatch /home/arendeiro/projects/chipmentation/src/pipelines/mapping_pipeline.sh $SAMPLE_NAME $SAMPLE_FILE
done < $SAMPLES_FILE

# Merge technical replicates
sh scripts/merge_replicates.sh # <- for two technical replicates
sh scripts/joining_replicates.sh # <- by hand for more than two replicates


### CORRELATIONS
# biological replicates
while read SAMPLE_NAME CONTROL_NAME; do
    sbatch /home/arendeiro/projects/chipmentation/src/scripts/counts_on_windows.sh $SAMPLE_NAME
    sbatch /home/arendeiro/projects/chipmentation/src/scripts/counts_on_windows.sh $CONTROL_NAME
done < $PEAKS_FILE

# merged replicates
while read SAMPLE_NAME CONTROL_NAME; do
    sbatch /home/arendeiro/projects/chipmentation/src/scripts/counts_on_windows_mergedSamples.sh $SAMPLE_NAME
    sbatch /home/arendeiro/projects/chipmentation/src/scripts/counts_on_windows_mergedSamples.sh $CONTROL_NAME
done < $PEAKS_FILE_MERGED

# concatenate files
sbatch /home/arendeiro/projects/chipmentation/src/scripts/concatenate_counts_on_windows.sh
# plot correlations
R scripts/plot_correlations.R # needs more memory, ask with salloc


### MERGE SAMPLES
# work on this...
join_technical_replicates.sh

while read SAMPLE_NAME_R1 SAMPLE_NAME_R2; do
    sbatch /home/arendeiro/projects/chipmentation/src/scripts/merge_biological_replicates.sh $SAMPLE_NAME_R1 $SAMPLE_NAME_R2
done < $REPLICATE_FILE

### PEAKS
# Call peaks with MACS2 and homer
while read SAMPLE_NAME CONTROL_FILE; do
    sbatch /home/arendeiro/projects/chipmentation/src/scripts/peak_calling.sh $SAMPLE_NAME $CONTROL_FILE
done < $PEAKS_FILE_MERGED

# with spp
while read SAMPLE_NAME CONTROL_NAME; do
    sbatch /home/arendeiro/projects/chipmentation/src/scripts/submit_call_peaks_spp.sh $SAMPLE_NAME $CONTROL_NAME
done < $PEAKS_FILE_MERGED



# OLD:
# Derive consensus peaks between replicates, count #'s
# R scripts/diffBind_analysis.R

# Compare peaks between samples
sh scripts/peak_intersection.sh


### FOOTPRINTING
## H3K4me3 on TSSs
# Get CAGE data
# Filter for high-scoring CAGE peaks
# Get window around CAGE summits with orientation
# Count reads, GC content and conservation on windows around CAGE summits
sbatch /home/arendeiro/projects/chipmentation/src/scripts/footprint_analysis_tss.sh H3K4me3_K562_500k_CM IgG_K562_500k_CM

source /home/arendeiro/venv/bin/activate
python /home/arendeiro/projects/chipmentation/src/scripts/plot_coverage_on_tss.py H3K4me3_K562_500k_CM
deactivate

# Filter on TATA-containing promoters
# Plot
# Plot xlim(-100,100)

## TFs
# Count reads on windows centered on peaks' summits
#run_counts_on_windows_peaks.sh # deprecated


# Find motifs
# Center peaks on motif, get 4kb window around 
# Count reads around motifs
sbatch /home/arendeiro/projects/chipmentation/src/scripts/footprint_analysis.sh PU1_K562_10mio_CM IgG_K562_10mio_CM
sbatch /home/arendeiro/projects/chipmentation/src/scripts/footprint_analysis.sh CTCF_K562_10mio_CM IgG_K562_10mio_CM

source /home/arendeiro/venv/bin/activate
python /home/arendeiro/projects/chipmentation/src/scripts/plot_coverage_on_tss.py PU1_K562_10mio_CM
deactivate

# Run CENTIPEDE

# Run pyDNASE


# Autocorrelations in H3K4me3 peaks
python scripts/distance_reads_correlations.py