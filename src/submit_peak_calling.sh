#!/bin/bash

SAMPLES_FILE=/home/arendeiro/projects/chipmentation/samples_peaks.txt

while read SAMPLE_NAME CONTROL_NAME; do
    sbatch /home/arendeiro/projects/chipmentation/src/peak_calling.sh \
    $SAMPLE_NAME \
    $CONTROL_NAME \
    --job-name="peakcalling_${SAMPLE_NAME}" \
    --error="/home/arendeiro/logs/peakcalling_${SAMPLE_NAME}_%j.err" \
    --output="/home/arendeiro/logs/peakcalling_${SAMPLE_NAME}_%j.log"
done < $SAMPLES_FILE