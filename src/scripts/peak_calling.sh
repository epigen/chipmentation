#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

# Optional parameters
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8000
#SBATCH --nodes=1
#SBATCH --job-name=call_peaks
#SBATCH --output=/home/arendeiro/logs/call_peaks_%j.out

# *** setup environment ***
# load the required environmental modules that your script depends upon
module load FastQC/0.11.2
module load trimmomatic/0.32
module load samtools
module load bamtools/bamtools
module load bowtie/2.2.0
module load picard-tools/1.118

# set the temporary folder to your scratch location (avoid using node local /tmp)
export SCRATCH_PATH=/fhgfs/scratch/users/arendeiro
export TMPDIR=$SCRATCH_PATH/tmp

# *** run the job ***
hostname
date

### Get sample info from arguments
SAMPLE_NAME=$1
CONTROL_NAME=$2

### Activate virtual environment
source /home/arendeiro/venv/bin/activate


### Specify paths
PROJECTDIR=/home/arendeiro/data/human/chipmentation
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19

### Call Peaks
mkdir -p $PROJECTDIR/peaks/

# MACS
if [[ $SAMPLE_NAME != *H3K27me3* ]]
    then
    echo "Calling broad peaks for sample " $SAMPLE_NAME

    macs2 callpeak -t $PROJECTDIR/mapped/merged/$SAMPLE_NAME.bam \
    -c $PROJECTDIR/mapped/merged/$CONTROL_NAME.bam \
    --bw 200 \
    -g hs -n ${SAMPLE_NAME} --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks

    awk '$9 >= 10' $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks/${SAMPLE_NAME}_peaks.narrowPeak > $PROJECTDIR/peaks/${SAMPLE_NAME}.narrowPeak
    awk '{print $1, $2 + $10, $2 + $10 + 1, $4}' $PROJECTDIR/peaks/${SAMPLE_NAME}.narrowPeak > $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.bed
    $BEDTOOLSDIR/slopBed -b 2000 -i $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.bed -g $GENOMESIZE > $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.2kb.bed
else
    macs2 callpeak -B -t $PROJECTDIR/mapped/merged/$SAMPLE_NAME.bam \
    -c $PROJECTDIR/mapped/merged/$CONTROL_NAME.bam \
    --broad --nomodel --extsize 200 --pvalue 1e-3 \
    -g hs -n ${SAMPLE_NAME} --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks
fi

deactivate
date


