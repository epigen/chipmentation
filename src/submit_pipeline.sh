#!/bin/bash

SAMPLES_FILE=/home/arendeiro/projects/chipmentation/samples.txt

while read SAMPLE_NAME SAMPLE_FILE; do
sbatch /home/arendeiro/projects/chipmentation/src/mapping_pipeline.sh $SAMPLE_NAME $SAMPLE_FILE
done < $SAMPLES_FILE