#!/bin/bash
export SCRATCH_PATH=/fhgfs/scratch/users/arendeiro
export TMPDIR=$SCRATCH_PATH/tmp

SAMPLES_FILE=/home/arendeiro/projects/chipmentation/samples_peaks.txt
PROJECTDIR=/home/arendeiro/data/human/chipmentation
GENOMESIZE=/fhgfs/prod/ngs_resources/genomes/hg19/hg19_chromLengths_sorted.txt
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19.fa
BEDTOOLSDIR=/home/arendeiro/.local/software/bedtools2/bin
HOMERDIR=/home/arendeiro/.local/software/homer-4.6/bin
CONSERVATION=/home/arendeiro/reference/Homo_sapiens/phyloP/placentalMammals/

### CAGE
# get data
cd /home/arendeiro/reference/Homo_sapiens/
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeRikenCageK562CellPapTssHmm.txt.gz
gzip -d wgEncodeRikenCageK562CellPapTssHmm.txt.gz
# filter on ~score, make bed
awk -v OFS='\t' '{(if $8 > 10) print $2, $3, $4, $5, $6, $7}' /home/arendeiro/reference/Homo_sapiens/wgEncodeRikenCageK562CellPapTssHmm.txt > /home/arendeiro/reference/Homo_sapiens/cage_K562_cell_tss_f10.bed

while read SAMPLE_NAME CONTROL_NAME; do
	# center on motifs
	annotatePeaks.pl $PROJECTDIR/peaks/${SAMPLE_NAME}.narrowPeak hg19 -size 4000 -center $PROJECTDIR/motifs/$SAMPLE_NAME/homerResults/motif1.motif | \
	awk -v OFS='\t' '{print $2, $3, $4, $1}' | \
	sortBed > $PROJECTDIR/bed/${SAMPLE_NAME}.motif.bed

	# get summits
	#awk -v OFS='\t' '{print $1, $2 + $10, $2 + $10 + 1, $4}' $PROJECTDIR/peaks/${SAMPLE_NAME}.narrowPeak > $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.bed
	# get 4kb window around summits
	#bedtools slopBed -b 2000 -i $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.bed -g $GENOMESIZE > $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.2kb.bed
	
	# get read count of PU1 in windows
	bedtools coverage -d -abam $PROJECTDIR/mapped/${SAMPLE_NAME}.bam -b $PROJECTDIR/bed/${SAMPLE_NAME}.motif.bed > $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage.bed
	# get read count of IgG in windows
	bedtools coverage -d -abam $PROJECTDIR/mapped/${CONTROL_NAME}.bam -b $PROJECTDIR/bed/${SAMPLE_NAME}.motif.bed > $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_IgG.bed
	
	# get nucleotide composition in windows
	## center at single nucleotide resolution
	## get 12bp windows around each nucleotide position in the 2kb windows
	## see nucleotide composition in those 12bp windows
	awk -v OFS='\t' '{print $1, $2 + $5, $2 + $5 + 1, $4}' $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage.bed | \
	bedtools slop -b 6 -i stdin -g $GENOMESIZE | \
	bedtools nuc -fi $GENOMEREF -bed stdin > $PROJECTDIR/bed/${SAMPLE_NAME}_peak_nucleotide_compos.bed

	# get conservation score
	# independently by chromossome
	for chr in `ls $CONSERVATION/*.bed`; \ 
    do bedmap --echo --echo-map-score $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage.bed $chr \
        > $TMPDIR/${SAMPLE_NAME}_peak_conservation.${chr}.bed; \
    done
    # concatenate and sort
	cat $TMPDIR/${SAMPLE_NAME}_peak_conservation.* | bedtools sort > $PROJECTDIR/bed/${SAMPLE_NAME}_peak_conservation.bed

	# For each position in each peak put all toghether (paste relevant columns from each file, remove headers when existing)
	cut -f 6 $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_IgG.bed | paste $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage.bed - > $TMPDIR/tmp
	cut -f 5,6,7,8,9,10 $PROJECTDIR/bed/${SAMPLE_NAME}_peak_nucleotide_compos.bed | tail -n +2 | paste $TMPDIR/tmp - > $TMPDIR/tmp2
	cut -f 5 $PROJECTDIR/bed/${SAMPLE_NAME}_peak_conservation.bed | paste $TMPDIR/tmp2 - > $PROJECTDIR/bed/${SAMPLE_NAME}_peak.2kb_coverage.bed
done < $SAMPLES_FILE



