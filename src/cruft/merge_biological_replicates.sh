#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

# Optional parameters
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10000
#SBATCH --nodes=1

#SBATCH --job-name=mergeReplicates
#SBATCH --output=/home/arendeiro/logs/mergeReplicates_%j.out


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
SAMPLE_NAME_R1=$1
SAMPLE_NAME_R2=$2

### Specify paths
PROJECTDIR=/home/arendeiro/data/human/chipmentation

echo $SAMPLE_NAME_R1
echo $SAMPLE_NAME_R2

# Get sample name
if [[ $SAMPLE_NAME_R1 == *_CM_* ]]
	then
	SAMPLE_NAME=`echo $SAMPLE_NAME_R1 | sed -e "s/_CM_.*/_CM/g"`
elif [[ $SAMPLE_NAME_R1 == *Encode* ]]
	then
	SAMPLE_NAME=`echo $SAMPLE_NAME_R1 | sed -e "s/_ChIP_Encode.*/_ChIP_Encode/g"`
else
	SAMPLE_NAME=`echo $SAMPLE_NAME_R1 | sed -e "s/_ChIP_.*/_ChIP/g"`
fi

echo $SAMPLE_NAME

if [[ $SAMPLE_NAME == *_CM* ]]
    then
    samtools merge -f $PROJECTDIR/mapped/merged/${SAMPLE_NAME}.bam $PROJECTDIR/mapped/${SAMPLE_NAME_R1}.trimmed.bowtie2.sorted.shifted.dup.bam $PROJECTDIR/mapped/${SAMPLE_NAME_R2}.trimmed.bowtie2.sorted.shifted.dup.bam
else
    samtools merge -f $PROJECTDIR/mapped/merged/${SAMPLE_NAME}.bam $PROJECTDIR/mapped/${SAMPLE_NAME_R1}.trimmed.bowtie2.sorted.dup.bam $PROJECTDIR/mapped/${SAMPLE_NAME_R2}.trimmed.bowtie2.sorted.dup.bam
fi

samtools index $PROJECTDIR/mapped/merged/${SAMPLE_NAME}.bam
samtools flagstat $PROJECTDIR/mapped/merged/${SAMPLE_NAME}.bam > $PROJECTDIR/mapped/merged/${SAMPLE_NAME}.flagstat

date


