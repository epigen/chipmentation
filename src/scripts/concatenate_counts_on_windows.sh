#!/bin/bash

SAMPLES_FILE=/home/arendeiro/projects/chipmentation/samples_peaks.txt
PROJECTDIR=/home/arendeiro/data/human/chipmentation
COUNTS=$PROJECTDIR/counts_1kb_windows.tsv

COUNTER=0
while read SAMPLE_NAME CONTROL_NAME; do
    echo "Doing sample: " $SAMPLE_NAME "and control: " $CONTROL_NAME
    cut -f 4 $PROJECTDIR/bed/${SAMPLE_NAME}_1kb_windows.bed | sed "1s/^/$SAMPLE_NAME\n/" > $PROJECTDIR/bed/${SAMPLE_NAME}_1kb_windows_1col.bed
    cut -f 4 $PROJECTDIR/bed/${CONTROL_NAME}_1kb_windows.bed | sed "1s/^/$CONTROL_NAME\n/" > $PROJECTDIR/bed/${CONTROL_NAME}_1kb_windows_1col.bed
    if [[ $COUNTER == 0 ]]
        then
        cp $PROJECTDIR/bed/${SAMPLE_NAME}_1kb_windows_1col.bed $COUNTS
        paste $COUNTS $PROJECTDIR/bed/${CONTROL_NAME}_1kb_windows_1col.bed > tmp
        mv tmp $COUNTS
    else
        paste $COUNTS $PROJECTDIR/bed/${SAMPLE_NAME}_1kb_windows_1col.bed > tmp
        mv tmp $COUNTS
        paste $COUNTS $PROJECTDIR/bed/${CONTROL_NAME}_1kb_windows_1col.bed > tmp
        mv tmp $COUNTS
    fi
    COUNTER=$[COUNTER + 1]
done < $SAMPLES_FILE
