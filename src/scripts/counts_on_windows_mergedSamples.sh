#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=48000
#SBATCH --nodes=1

#SBATCH --job-name=countsOnWindows10kb
#SBATCH --output=/home/arendeiro/logs/countsOnWindows10kb_%j.log

# *** setup environment ***
# load the required environmental modules that your script depends upon

# set the temporary folder to your scratch location (avoid using node local /tmp)
export SCRATCH_PATH=/fhgfs/scratch/users/arendeiro
export TMPDIR=$SCRATCH_PATH/tmp

# *** run the job ***
hostname
date

### Get sample info from arguments
SAMPLE_NAME=$1

### Get window info from arguments
WINDOWS=$2

### Specify paths
PROJECTDIR=/home/arendeiro/data/human/chipmentation

mkdir -p $PROJECTDIR/bed
mkdir -p $PROJECTDIR/bed/correlations

# needs a lot of memory with 1kb windows
if [[ $SAMPLE_NAME == *Encode* ]]; then
	echo "Counting reads on Encode sample: " $SAMPLE_NAME
	bedtools bamtobed -i $PROJECTDIR/../encode/chip-seq/mapped/merged/${SAMPLE_NAME}.nodups.bam | \
	bedtools intersect -c -a $WINDOWS -b stdin | \
	bedtools sort \
	> $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows.bed
	
	cut -f 4 $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows.bed | \
	sed "1s/^/$SAMPLE_NAME\n/" \
	> $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows_1col.bed
else
	echo "Counting reads on chipmentation sample: " $SAMPLE_NAME
	bedtools bamtobed -i $PROJECTDIR/mapped/merged/${SAMPLE_NAME}.nodups.bam | \
	bedtools intersect -c -a $WINDOWS -b stdin | \
	bedtools sort \
	> $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows.bed
	
	cut -f 4 $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows.bed | \
	sed "1s/^/$SAMPLE_NAME\n/" \
	> $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows_1col.bed
fi


date
