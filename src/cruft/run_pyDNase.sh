#!/bin/bash
#SBATCH --partition=longq
#SBATCH --ntasks=1
#SBATCH --time=42:00:00

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16000
#SBATCH --nodes=1

#SBATCH --job-name=pyDNase
#SBATCH --output=/home/arendeiro/logs/pyDNase_%j.out

# *** setup environment ***
# load the required environmental modules that your script depends upon
module load FastQC/0.11.2
module load trimmomatic/0.32
module load samtools
module load bamtools/bamtools
module load python/2.7.6

# set the temporary folder to your scratch location (avoid using node local /tmp)
export SCRATCH_PATH=/fhgfs/scratch/users/arendeiro
export TMPDIR=$SCRATCH_PATH/tmp

# *** run the job ***
hostname
date

### Get sample info from arguments
SAMPLE_NAME=$1
SS=$2
SE=$3
FS=$4
FE=$5

### Specify paths
PROJECTDIR=/home/arendeiro/data/human/chipmentation

echo "Doing sample" $SAMPLE_NAME ", with options: ${SS}_${SE}_${FS}_${FE}"

# make output dir
mkdir -p $PROJECTDIR/pyDNase/${SAMPLE_NAME}_${SS}_${SE}_${FS}_${FE}

# run
wellington_footprints.py -b -fp ${FS},${FE},2 -sh ${SS},${SE},1 \
$PROJECTDIR/peaks/${SAMPLE_NAME}_peaks/${SAMPLE_NAME}_peaks.bed \
$PROJECTDIR/mapped/merged/${SAMPLE_NAME}.bam \
$PROJECTDIR/pyDNase/${SAMPLE_NAME}_${SS}_${SE}_${FS}_${FE}

date
