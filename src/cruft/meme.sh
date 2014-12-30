#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=12:00:00

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32000
#SBATCH --nodes=1

#SBATCH --job-name=meme
#SBATCH --output=/home/arendeiro/logs/meme%j.out

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
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19.fa

mkdir -p $PROJECTDIR/motifs_meme
mkdir -p $PROJECTDIR/motifs_meme/$SAMPLE_NAME

# meme de novo find motifs
bedtools getfasta \
-fi $GENOMEREF \
-bed $PROJECTDIR/bed/${SAMPLE_NAME}.motif.bed \
-fo $PROJECTDIR/motifs_meme/$SAMPLE_NAME/${SAMPLE_NAME}_motifs.fa \
-name


meme \
$PROJECTDIR/motifs_meme/$SAMPLE_NAME/${SAMPLE_NAME}_motifs.fa \
-o $PROJECTDIR/motifs_meme/$SAMPLE_NAME/ \
-text \
-dna \
-maxsize 1000000000


date
