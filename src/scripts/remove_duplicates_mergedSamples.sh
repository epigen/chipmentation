#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000
#SBATCH --nodes=1

#SBATCH --job-name=countsOnWindows1kb
#SBATCH --output=/home/arendeiro/logs/countsOnWindows1kb_%j.log

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

### Get sample info from arguments
SAMPLE_NAME=$1

echo "Removing duplicates for sample: $SAMPLE_NAME"
sambamba markdup -t 16 -r $PROJECTDIR/mapped/merged/${SAMPLE_NAME}.bam $PROJECTDIR/mapped/merged/${SAMPLE_NAME}.nodups.bam
sambamba index -t 16 $PROJECTDIR/mapped/merged/${SAMPLE_NAME}.nodups.bam

date
