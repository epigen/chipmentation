#!/bin/bash
#SBATCH --partition=longq
#SBATCH --ntasks=1
#SBATCH --time=72:00:00

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32000
#SBATCH --nodes=1

#SBATCH --job-name=footprintingPeaks
#SBATCH --output=/home/arendeiro/logs/footprintingPeaks_%j.out

PROJECTDIR=/home/arendeiro/data/human/chipmentation/

SAMPLE_NAME=$1

findMotifsGenome.pl \
$PROJECTDIR/peaks/${SAMPLE_NAME}_peaks/${SAMPLE_NAME}_peaks.narrowPeak \
hg19 \
$PROJECTDIR/motifs/${SAMPLE_NAME}/ \
-mask -size 150 -len 8,10,12

date