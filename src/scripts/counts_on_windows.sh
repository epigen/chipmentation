#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=48000
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

### Get sample info from arguments
SAMPLE_NAME=$1

### Specify paths
PROJECTDIR=/home/arendeiro/data/human/chipmentation
BEDTOOLSDIR=/home/arendeiro/.local/software/bedtools2/bin
GENOMESIZE=/fhgfs/prod/ngs_resources/genomes/hg19/hg19_chromLengths_sorted.txt

mkdir -p $PROJECTDIR/bed
mkdir -p $PROJECTDIR/bed/correlations

# to make windows run:
# bedtools makewindows -g $GENOMESIZE -w 1000 > /home/arendeiro/reference/Homo_sapiens/1kb_windows.bed
# bedtools makewindows -g $GENOMESIZE -w 10000 > /home/arendeiro/reference/Homo_sapiens/10kb_windows.bed
WINDOWS=/home/arendeiro/reference/Homo_sapiens/1kb_windows.bed

# needs a lot of memory with 1kb windows
if [[ $SAMPLE_NAME == *Encode* ]]; then
	echo "Counting reads on Encode sample: " $SAMPLE_NAME
	bedtools bamtobed -i $PROJECTDIR/../encode/chip-seq/mapped/$SAMPLE_NAME.bam | \
	bedtools intersect -c -a $WINDOWS -b stdin | \
	bedtools sort \
	> $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows.bed
	
	cut -f 4 $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows.bed | \
	sed "1s/^/$SAMPLE_NAME\n/" \
	> $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows_1col.bed

else
	if [[ $SAMPLE_NAME == *_CM_* ]]; then
		echo "Counting reads on chipmentation sample: " $SAMPLE_NAME
		bedtools bamtobed -i $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam | \
		bedtools intersect -c -a $WINDOWS -b stdin | \
		bedtools sort \
		> $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows.bed
		
		cut -f 4 $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows.bed | \
		sed "1s/^/$SAMPLE_NAME\n/" \
		> $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows_1col.bed

	elif [[ $SAMPLE_NAME == *_ChIP_* ]]; then
		echo "Counting reads on ChIP sample: " $SAMPLE_NAME
		bedtools bamtobed -i $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.bam | \
		bedtools intersect -c -a $WINDOWS -b stdin | \
		bedtools sort \
		> $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows.bed
		
		cut -f 4 $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows.bed | \
		sed "1s/^/$SAMPLE_NAME\n/" \
		> $PROJECTDIR/bed/correlations/${SAMPLE_NAME}_1kb_windows_1col.bed

	fi
fi


date