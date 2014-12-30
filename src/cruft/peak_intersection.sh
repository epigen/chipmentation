#!/bin/bash


PROJECTDIR=/home/arendeiro/data/human/chipmentation
INTERSECTIONS=/home/arendeiro/projects/chipmentation/results/peak_intersections.txt

SAMPLES_FILE=/home/arendeiro/projects/chipmentation/samples_peaks.txt

SAMPLE_CM=$1

for SAMPLE_NAME in 
    
    SAMPLE=${SAMPLE_NAME/_CM/}
    SAMPLE_CM=${SAMPLE}_CM
    SAMPLE_CHIP=${SAMPLE}_ChIP

    # Merge peaks within samples that are <100bp apart
    echo "Merging peaks in" $SAMPLE_CM
    bedtools merge -d 100 \
    -i $PROJECTDIR/peaks/${SAMPLE_CM}_peaks/${SAMPLE_CM}_peaks.narrowPeak \
    > $PROJECTDIR/peaks/${SAMPLE_CM}_peaks/${SAMPLE_CM}_peaks.merged.bed

    echo "Merging peaks in" $SAMPLE_CHIP
    bedtools merge -d 100 \
    -i $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak \
    > $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.merged.bed

    # Count total number of peaks
    echo "Counting peaks in" $SAMPLE_CM
    CM_COUNT="$(wc -l $PROJECTDIR/peaks/${SAMPLE_CM}_peaks/${SAMPLE_CM}_peaks.merged.bed)"
    echo "Counting peaks in" $SAMPLE_CHIP
    CHIP_COUNT="$(wc -l $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.merged.bed)"

    # Intersect and count
    echo "Intersecting" $SAMPLE_CM "with" $SAMPLE_CHIP
    CM_CHIP="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM}_peaks/${SAMPLE_CM}_peaks.merged.bed \
    -b $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.merged.bed \
    | wc -l)"

    echo "Intersecting" $SAMPLE_CHIP "with" $SAMPLE_CM
    CHIP_CM="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.merged.bed \
    -b $PROJECTDIR/peaks/${SAMPLE_CM}_peaks/${SAMPLE_CM}_peaks.merged.bed \
    | wc -l)"

    echo $SAMPLE_CM $CM_COUNT $CHIP_COUNT $CM_CHIP $CHIP_CM >> $INTERSECTIONS

done < $SAMPLES_FILE

column $INTERSECTIONS > tmp
mv tmp $INTERSECTIONS