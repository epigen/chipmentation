#!/bin/bash
#SBATCH --partition=develop
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16000
#SBATCH --nodes=1

#SBATCH --job-name=homer_findPeaks_%j
#SBATCH --output=/home/arendeiro/logs/homer_findPeaks_%j.out

# *** setup environment ***
# load the required environmental modules that your script depends upon
module load python/2.7.6

# set the temporary folder to your scratch location (avoid using node local /tmp)
export SCRATCH_PATH=/fhgfs/scratch/users/arendeiro
export TMPDIR=$SCRATCH_PATH/tmp

# *** run the job ***
hostname
date

### Get sample info from arguments
SAMPLE_NAME=$1
CONTROL_NAME=$2

### Specify paths
PROJECTDIR=/home/arendeiro/data/human/chipmentation
HOMERDIR=/home/arendeiro/.local/software/homer-4.6/bin

# HOMER
# find peaks
if [[ $SAMPLE_NAME == *H3K4* ]]
    then
    echo "Running with settings for H3K4me3 for sample: " $SAMPLE_NAME
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style histone -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer
elif [[ $SAMPLE_NAME == *H3K27* ]]
    then
    echo "Running with settings for H3K27me3 for sample: " $SAMPLE_NAME
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style histone -region -size 1000 -minDist 10000 -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer
else
    echo "Running with settings for TFs for sample: " $SAMPLE_NAME
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style factor -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer
fi

# make bed file
sed '/^#/d' $PROJECTDIR/homer/${SAMPLE_NAME}_homer/regions.txt | \
cut -f 2,3,4,11,12,13 > $PROJECTDIR/homer/${SAMPLE_NAME}_homer/${SAMPLE_NAME}_homer_peaks.bed

date