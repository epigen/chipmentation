#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=48000
#SBATCH --nodes=1

#SBATCH --job-name=concatenateCountsOnWindows
#SBATCH --output=/home/arendeiro/logs/concatenateCountsOnWindows_%j.log

# *** setup environment ***
# load the required environmental modules that your script depends upon

# set the temporary folder to your scratch location (avoid using node local /tmp)
export SCRATCH_PATH=/fhgfs/scratch/users/arendeiro
export TMPDIR=$SCRATCH_PATH/tmp

# *** run the job ***
hostname
date

### Specify paths
PROJECTDIR=/home/arendeiro/data/human/chipmentation

# 1kb windows
COUNTS=$PROJECTDIR/counts_1kb_windows.nodups.tsv

# Samples
ALL_SAMPLES=$1

COUNTER=0
while read SAMPLE_NAME CONTROL_NAME; do
    echo "Doing sample: " $SAMPLE_NAME "and control: " $CONTROL_NAME
    if [[ $COUNTER == 0 ]]
        then
        cp $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows_1col.bed $COUNTS
    else
        paste $COUNTS $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows_1col.bed > tmp
        mv tmp $COUNTS
    fi
    COUNTER=$[COUNTER + 1]
done < $ALL_SAMPLES

date


