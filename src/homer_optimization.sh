#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --nodes=1

#SBATCH --job-name=homer_findPeaks
#SBATCH --output=/home/arendeiro/logs/homer_findPeaks_%j.out
#SBATCH --error=/home/arendeiro/logs/homer_findPeaks_%j.err

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
if [[ $SAMPLE_NAME == *H3* ]]
    then
    echo "Detected histone data"
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style histone -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer_std
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style histone -region -size 150 -minDist 370 -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer_region_150_370
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style histone -region -size 150 -minDist 1000 -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer_region_150_1000
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style histone -region -size 1000 -minDist 2500 -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer_region_1000_2500
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style histone -region -size 1000 -minDist 10000 -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer_region_1000_10000
else
    echo "Detected TF data"
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style factor -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer_std
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style factor -C 1 -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer_c1
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style factor -C 0 -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer_c0
fi

date