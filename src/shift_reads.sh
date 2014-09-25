#!/bin/bash
#SBATCH --partition=develop
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --nodes=1

#SBATCH --job-name=shiftReads
#SBATCH --output=/home/arendeiro/logs/shiftReads_%j.out
#SBATCH --error=/home/arendeiro/logs/shiftReads_%j.err


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
# $    --job-name="pipeline_trim_${SAMPLE_NAME}" \
# $    --error="/home/arendeiro/logs/pipeline_trim_${SAMPLE_NAME}_%j.err" \
# $    --output="/home/arendeiro/logs/pipeline_trim_${SAMPLE_NAME}_%j.log" \
# $    $SAMPLE_NAME \
# $    $SAMPLE_FILE
# $ done < $SAMPLES_FILE


# *** setup environment ***
# load the required environmental modules that your script depends upon
module load FastQC/0.11.2
module load trimmomatic/0.32
module load samtools
module load bamtools/bamtools
module load bowtie/2.2.0

# set the temporary folder to your scratch location (avoid using node local /tmp)
export SCRATCH_PATH=/fhgfs/scratch/users/arendeiro
export TMPDIR=$SCRATCH_PATH/tmp

# *** run the job ***
hostname
date

### Get sample info from arguments
SAMPLE_NAME=$1
SAMPLE_FILE=$2

### Specify paths
PROJECTDIR=/home/arendeiro/data/human/chipmentation
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19
PRESEQ=/home/arendeiro/.local/software/preseq-0.1.0/preseq

# Shift reads from tagmentation samples
 if [[ $SAMPLE_NAME == *ATAC-seq* ]] || [[ $SAMPLE_NAME == *CM* ]]
    then
    echo "Shifting reads for sample:" $SAMPLE_NAME
    samtools view $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.bam | awk -F "\t" 'BEGIN{OFS="\t"}($2==16){ offset=$4-5; print $1,$2,$3, offset,$5,$6,$7,$8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}' > $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted_minus.sam
    samtools view $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.bam | awk -F "\t" 'BEGIN{OFS="\t"}($2==0){ offset=$4+4; print $1,$2,$3, offset,$5,$6,$7,$8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}' > $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted_plus.sam
    cat $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted_minus.sam $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted_plus.sam | samtools view -S -b - | samtools sort - $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.bam
    if [ -s $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.bam ]
    	then
        rm $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted_minus.sam $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted_plus.sam
    fi
fi

date