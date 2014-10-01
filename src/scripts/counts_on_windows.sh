#!/bin/bash
#SBATCH --partition=develop
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000
#SBATCH --nodes=1

#SBATCH --job-name=countsOnWindows
#SBATCH --output=/home/arendeiro/logs/countsOnWindows_%j.log

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

### Specify paths
PROJECTDIR=/home/arendeiro/data/human/chipmentation
BEDTOOLSDIR=/home/arendeiro/.local/software/bedtools2/bin
GENOMESIZE=/fhgfs/prod/ngs_resources/genomes/hg19/hg19_chromLengths_sorted.txt

mkdir -p $PROJECTDIR/bed/

# $BEDTOOLSDIR/windowMaker -g $GENOMESIZE -w 1000 > /home/arendeiro/reference/Homo_sapiens/1kb_windows.bed
# $BEDTOOLSDIR/windowMaker -g $GENOMESIZE -w 10000 > /home/arendeiro/reference/Homo_sapiens/10kb_windows.bed
WINDOWS=/home/arendeiro/reference/Homo_sapiens/1kb_windows.bed

# needs a lot of memory with 1kb windows
if [[ $SAMPLE_NAME == *Encode* ]]; then
	$BEDTOOLSDIR/bamToBed -i $PROJECTDIR/../encode/chip-seq/mapped/$SAMPLE_NAME.bam | $BEDTOOLSDIR/intersectBed -c -a $WINDOWS -b stdin | $BEDTOOLSDIR/sortBed > $PROJECTDIR/bed/${SAMPLE_NAME}_1kb_windows.bed
	cut -f 4 $PROJECTDIR/bed/${SAMPLE_NAME}_1kb_windows.bed | sed "1s/^/$SAMPLE_NAME\n/" > $PROJECTDIR/bed/${SAMPLE_NAME}_1kb_windows_1col.bed
else
	$BEDTOOLSDIR/bamToBed -i $PROJECTDIR/mapped/$SAMPLE_NAME.bam | $BEDTOOLSDIR/intersectBed -c -a $WINDOWS -b stdin | $BEDTOOLSDIR/sortBed > $PROJECTDIR/bed/${SAMPLE_NAME}_10kb_windows.bed
	cut -f 4 $PROJECTDIR/bed/${SAMPLE_NAME}_10kb_windows.bed | sed "1s/^/$SAMPLE_NAME\n/" > $PROJECTDIR/bed/${SAMPLE_NAME}_10kb_windows_1col.bed
fi


date