#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --nodes=1

#SBATCH --job-name=fastqc
#SBATCH --output=/home/arendeiro/logs/fastqc_%j.out
#SBATCH --error=/home/arendeiro/logs/fastqc_%j.err

# *** setup environment ***
# load the required environmental modules that your script depends upon
module load java/jdk/1.7.0_65 
module load FastQC/0.11.2

# set the temporary folder to your scratch location (avoid using node local /tmp)
export SCRATCH_PATH=/fhgfs/scratch/users/arendeiro
export TMPDIR=$SCRATCH_PATH/tmp

# *** run the job ***
hostname
date

### Get sample info from arguments
SAMPLE_FILE=$1

### Specify paths
PROJECTDIR=/home/arendeiro/data/human/chipmentation
ADAPTERS=/home/arendeiro/adapters/chipmentation.fa

### Start work on samples 
# Fastqc on raw data
mkdir -p $PROJECTDIR/raw/fastqc

fastqc $SAMPLE_FILE \
-format bam \
--noextract \
--contaminants $ADAPTERS \
--outdir $PROJECTDIR/raw/fastqc

date
