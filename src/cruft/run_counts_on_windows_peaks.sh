#!/bin/bash

SAMPLES_FILE=/home/arendeiro/projects/chipmentation/samples_peaks.txt
PROJECTDIR=/home/arendeiro/data/human/chipmentation
GENOMESIZE=/fhgfs/prod/ngs_resources/genomes/hg19/hg19_chromLengths_sorted.txt
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19.fa
BEDTOOLSDIR=/home/arendeiro/.local/software/bedtools2/bin

while read SAMPLE_NAME CONTROL_NAME; do
	# filter for qvalue
	awk '$9 >= 10' $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_bw200/${SAMPLE_NAME}_bw200_peaks.narrowPeak > $PROJECTDIR/peaks/${SAMPLE_NAME}.narrowPeak
	# get summits
	awk -v OFS='\t' '{print $1, $2 + $10, $2 + $10 + 1, $4}' $PROJECTDIR/peaks/${SAMPLE_NAME}.narrowPeak > $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.bed
	# get 4kb window around summits
	$BEDTOOLSDIR/slopBed -b 2000 -i $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.bed -g $GENOMESIZE > $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.2kb.bed
	# get read count of PU1 in windows
	$BEDTOOLSDIR/coverageBed -d -abam $PROJECTDIR/mapped/${SAMPLE_NAME}.bam -b $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.2kb.bed > $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage.bed
	# get read count of IgG in windows
	$BEDTOOLSDIR/coverageBed -d -abam $PROJECTDIR/mapped/${CONTROL_NAME}.bam -b $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.2kb.bed > $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_IgG.bed
	# get nucleotide composition in windows
	## center at single nucleotide resolution
	## get 12bp windows around each nucleotide position in the 2kb windows
	## see nucleotide composition in those 12bp windows
	awk -v OFS='\t' '{print $1, $2 + $5, $2 + $5 + 1, $4}' $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage.bed | \
	$BEDTOOLSDIR/slopBed -b 6 -i stdin -g $GENOMESIZE | \
	$BEDTOOLSDIR/nucBed -fi $GENOMEREF -bed stdin > $PROJECTDIR/bed/${SAMPLE_NAME}_peak_nucleotide_compos.bed

	# For each position in each peak put all toghether (paste relevant columns from each file, remove headers when existing)
	cut -f 6 $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_IgG.bed | paste $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage.bed - > tmp
	cut -f 5,6,7,8,9,10 $PROJECTDIR/bed/${SAMPLE_NAME}_peak_nucleotide_compos.bed | tail -n +2 | paste tmp - > $PROJECTDIR/bed/${SAMPLE_NAME}_peak.2kb_coverage.bed
done < $SAMPLES_FILE
