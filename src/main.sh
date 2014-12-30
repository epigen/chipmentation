#!/bin/bash
#########################
# ChIPmentation project #
#########################

SRCDIR=/home/arendeiro/projects/chipmentation/
DATADIR=/home/arendeiro/data/human/chipmentation/

SAMPLES_FILE=$SRCDIR/samples.txt
REPLICATE_FILE=$SRCDIR/samples_to_merge.txt

PEAKS_FILE=$SRCDIR/samples_peaks.txt
PEAKS_FILE_MERGED=$SRCDIR/samples_peaks_merged.txt

ALL_SAMPLES=$SRCDIR/all_samples_for_correlation.txt


### PREPROCESSING
# External data fetched as in:
# $SRCDIR/scripts/getEncodeData.sh
# $SRCDIR/scripts/getCageData.sh

# Fastqc
while read SAMPLE_NAME SAMPLE_FILE; do
    sbatch $SRCDIR/scripts/fastqc.sh $SAMPLE_FILE
done < $SAMPLES_FILE

# Trim, map, mark duplicates, shift reads, qc
while read SAMPLE_NAME SAMPLE_FILE; do
    sbatch $SRCDIR/pipelines/mapping_pipeline.sh $SAMPLE_NAME $SAMPLE_FILE
done < $SAMPLES_FILE


### MERGE REPLICATES
# Merge technical replicates
sh $SRCDIR/scripts/merge_technical_replicates.sh # <- by hand - needs work

# Merge biological replicates
while read REP1 REP2; do
    sbatch $SRCDIR/scripts/merge_biological_replicates.sh $REP1 $REP2
done < $REPLICATE_FILE


### UCSC TRACKS
# make bigwig tracks for merged replicates
while read SAMPLE_NAME CONTROL_NAME; do
    sbatch $SRCDIR/scripts/bamToBigWig.sh $DATADIR/mapped/merged/$SAMPLE_NAME
    sbatch $SRCDIR/scripts/bamToBigWig.sh $DATADIR/mapped/merged/$CONTROL_NAME
done < $PEAKS_FILE_MERGED

# bigwig files:     http://www.biomedical-sequencing.at/bocklab/arendeiro/chipmentation2014/bigWigTracks/
# trackHub file:    http://www.biomedical-sequencing.at/bocklab/arendeiro/chipmentation2014/bigWigTracks/tracks_hub.txt
# UCSC link:        http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&position=chr22&hgt.customText=http://www.biomedical-sequencing.at/bocklab/arendeiro/chipmentation2014/bigWigTracks/tracks_hub.txt


### CORRELATIONS
## make bams without duplicates
while read SAMPLE_NAME CONTROL_NAME; do
    sbatch $SRCDIR/scripts/remove_duplicates.sh $SAMPLE_NAME
    sbatch $SRCDIR/scripts/remove_duplicates.sh $CONTROL_NAME
done < $PEAKS_FILE

while read SAMPLE_NAME CONTROL_NAME; do
    sbatch $SRCDIR/scripts/remove_duplicates_mergedSamples.sh $SAMPLE_NAME
    sbatch $SRCDIR/scripts/remove_duplicates_mergedSamples.sh $CONTROL_NAME
done < $PEAKS_FILE_MERGED

# make genome-wide 1kb windows:
GENOMESIZE=/fhgfs/prod/ngs_resources/genomes/hg19/hg19_chromLengths_sorted.txt
WINDOWS=/home/arendeiro/reference/Homo_sapiens/1kb_windows.bed
bedtools makewindows -g $GENOMESIZE -w 1000 > $WINDOWS

# count reads for biological replicates
while read SAMPLE_NAME CONTROL_NAME; do
    sbatch $SRCDIR/scripts/counts_on_windows.sh $SAMPLE_NAME $WINDOWS
    sbatch $SRCDIR/scripts/counts_on_windows.sh $CONTROL_NAME $WINDOWS
done < $PEAKS_FILE

# count reads for merged replicates
while read SAMPLE_NAME CONTROL_NAME; do
    sbatch $SRCDIR/scripts/counts_on_windows_mergedSamples.sh $SAMPLE_NAME $WINDOWS
    sbatch $SRCDIR/scripts/counts_on_windows_mergedSamples.sh $CONTROL_NAME $WINDOWS
done < $PEAKS_FILE_MERGED

# concatenate files
cat $SAMPLES_FILE $PEAKS_FILE_MERGED > $ALL_SAMPLES
sbatch $SRCDIR/scripts/concatenate_counts_on_windows.sh $ALL_SAMPLES

# plot correlations
# needs more memory: ask with salloc:
# salloc --partition=develop --time=10:00:00 --job-name=R --nodes=1 --mem=20000 --nodelist=n002 srun R --no-save
R scripts/plot_correlations.R 


### CALL PEAKS
# Call peaks with MACS2 on merged replicates
while read SAMPLE_NAME CONTROL_FILE; do
    sbatch $SRCDIR/scripts/macs2_peak_calling.sh $SAMPLE_NAME $CONTROL_FILE
done < $PEAKS_FILE_MERGED

# Call peaks with SPP on merged replicates
while read SAMPLE_FILE CONTROL_FILE; do
    sbatch $SRCDIR/scripts/submit_call_peaks_spp.sh \
    $DATADIR/mapped/merged/${SAMPLE_FILE}.bam \
    $DATADIR/mapped/merged/${CONTROL_FILE}.bam
done < $PEAKS_FILE_MERGED

# Compare peaks between samples
# terrible brace yourself!
sh peak_intersection_histone.sh 
sh peak_intersection_TFs.sh
sh peak_intersection_H3K27me3.sh # <- Terrible +1


### H3K4me3 on TSSs
# Get signal coverage over CAGE TSSs (120bp) at 1bp resolution
SAMPLE=H3K4me3_K562_500k
CAGE=hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed
for SIGNAL in CM IgG ChIP DNase
do
    sbatch $SRCDIR/scripts/cage_tss_coverage_job.sh $SAMPLE $CAGE $SIGNAL
done

# Parse signal (also inverts negative strand signal)
for SIGNAL in CM IgG ChIP DNase
do
    sbatch $SRCDIR/scripts/cage_tss_coverage-pythonParse_job.sh $CAGE $SIGNAL
done

# Plot average profiles around CAGE TSSs
R TSS_plots.R

# Cluster TSSs based on signals from several techniques
R TSS_heatmap_square.R

# Heatmaps were exported from javatreeview manually (Christian)

# TSS motifs enrichment p-values
# add code from dell ~/


### TF FOOTPRINTING
# Find motifs
for SAMPLE in PU1_K562_10mio_CM CTCF_K562_10mio_CM
    sbatch $SRCDIR/scripts/findMotifs.sh $SAMPLE
done

for SAMPLE_NAME,CONTROL_NAME in PU1_K562_10mio_CM CTCF_K562_10mio_CM
    sbatch $SRCDIR/scripts/footprint_analysis.sh $SAMPLE_NAME
done

# Plot average profiles, footprint with Centipede, export CDTs
R TF_plots.R

