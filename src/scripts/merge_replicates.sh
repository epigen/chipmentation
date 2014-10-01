#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

# Optional parameters
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000
#SBATCH --nodes=1

#SBATCH --job-name=mergeReplicates
#SBATCH --output=/home/arendeiro/logs/mergeReplicates_%j.out
#SBATCH --error=/home/arendeiro/logs/mergeReplicates_%j.err


# *** setup environment ***
# load the required environmental modules that your script depends upon
module load samtools

# set the temporary folder to your scratch location (avoid using node local /tmp)
export SCRATCH_PATH=/fhgfs/scratch/users/arendeiro
export TMPDIR=$SCRATCH_PATH/tmp

# *** run the job ***
hostname
date

### Get sample info from arguments
SAMPLE_NAME=$1

### Specify paths
PROJECTDIR=/home/arendeiro/data/human/chipmentation

if [[ $SAMPLE_NAME == *ATAC-seq* ]] || [[ $SAMPLE_NAME == *_CM_* ]]
    then
    samtools merge -f $PROJECTDIR/mapped/${SAMPLE_NAME}.bam $PROJECTDIR/mapped/${SAMPLE_NAME}-1.trimmed.bowtie2.sorted.shifted.dup.bam $PROJECTDIR/mapped/${SAMPLE_NAME}-2.trimmed.bowtie2.sorted.shifted.dup.bam
else
    samtools merge -f $PROJECTDIR/mapped/${SAMPLE_NAME}.bam $PROJECTDIR/mapped/${SAMPLE_NAME}-1.trimmed.bowtie2.sorted.dup.bam $PROJECTDIR/mapped/${SAMPLE_NAME}-2.trimmed.bowtie2.sorted.dup.bam
fi
date