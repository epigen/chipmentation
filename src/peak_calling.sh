#!/bin/bash
#SBATCH --partition=develop
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

# Optional parameters
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2000
#SBATCH --nodes=1
#SBATCH --mail-type=end
#SBATCH --mail-user=arendeiro@cemm.oeaw.ac.at


####### EXPLANATION OF USAGE ##########

# Submit jobs for this script by doing:
# $ sbatch <path>/pipeline_with_trimming.sh SAMPLE_NAME SAMPLE_FILE
# 
# You can encapsulate this in a loop as such to run multiple samples from a file:
#
# SAMPLES_FILE=/home/arendeiro/projects/chipmentation/samples_test_4.txt
# $ while read SAMPLE_NAME SAMPLE_FILE; do
# $     sbatch /home/arendeiro/projects/chipmentation/src/pipeline_with_trimming.sh $SAMPLE_NAME $SAMPLE_FILE
# $ done < $SAMPLES_FILE
# 
# You can supply additional options to the sbatch command, so that the jobname and logs
# reflect the current sample:

# SAMPLES_FILE=/home/arendeiro/projects/chipmentation/samples_test_4.txt
#
# $ while read SAMPLE_NAME SAMPLE_FILE; do
# $    sbatch /home/arendeiro/projects/chipmentation/src/pipeline_with_trimming.sh \
# $    $SAMPLE_NAME \
# $    $SAMPLE_FILE \
# $    --job-name="pipeline_trim_${SAMPLE_NAME}" \
# $    --error="/home/arendeiro/logs/pipeline_trim_${SAMPLE_NAME}_%j.err" \
# $    --output="/home/arendeiro/logs/pipeline_trim_${SAMPLE_NAME}_%j.log" \
# $ done < $SAMPLES_FILE


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
SAMPLE_FILE=$2
CONTROL_NAME=$3
CONTROL_FILE=$4

### Specify paths
PROJECTDIR=/home/arendeiro/data/human/chipmentation
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19
PICARDDIR=/cm/shared/apps/picard-tools/1.118
PRESEQ=/home/arendeiro/.local/software/preseq-0.1.0/preseq
MACS2="python2.7 .local/software/bin/macs2"
HOMERDIR=/home/arendeiro/.local/software/homer-4.6/bin

### Call Peaks
mkdir -p $PROJECTDIR/trimmed/peaks/

# MACS
$MACS2 callpeak -t $SAMPLE_FILE \
-c $CONTROL_FILE \
-g hs -n $SAMPLE_NAME --outdir $PROJECTDIR/trimmed/peaks/${SAMPLE_NAME}_peaks_MACS2

# HOMER
# make tag dirs
mkdir -p $PROJECTDIR/homer/${SAMPLE_NAME}_homer
$HOMERDIR/makeTagDirectory $PROJECTDIR/homer/${SAMPLE_NAME}_homer $SAMPLE_FILE
$HOMERDIR/makeTagDirectory $PROJECTDIR/homer/${CONTROL_NAME}_homer $CONTROL_FILE
# make UCSC tracks
$HOMERDIR/makeUCSCfile $PROJECTDIR/homer/${SAMPLE_NAME}_homer -name $SAMPLE_NAME -o auto
$HOMERDIR/makeUCSCfile $PROJECTDIR/homer/${CONTROL_NAME}_homer -name $CONTROL_NAME -o auto

# find peaks
if [[ $SAMPLE_NAME == *H3* ]]
    then
    echo "Detected histone data"
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style histone -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer
else
    echo "Detected TF data"
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style factor -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer
fi

date