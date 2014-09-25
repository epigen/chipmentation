#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --nodes=1

#SBATCH --job-name=shiftreads_markdups
#SBATCH --output=/home/arendeiro/logs/shiftreads_markdups_%j.out
#SBATCH --error=/home/arendeiro/logs/shiftreads_markdups_%j.err


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
module load python/2.7.6

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
PICARDDIR=/cm/shared/apps/picard-tools/1.118/

# Shift reads from tagmentation samples
if [[ $SAMPLE_NAME == *ATAC-seq* ]] || [[ $SAMPLE_NAME == *_CM_* ]]
    then
    echo "Shifting reads for sample:" $SAMPLE_NAME
    samtools view -h $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.bam | python2.7 projects/chipmentation/src/shift_reads.py | samtools view -S -b - | samtools sort - $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted
fi

# Mark & remove duplicates
if [[ $SAMPLE_NAME == *ATAC-seq* ]] || [[ $SAMPLE_NAME == *_CM_* ]]
    then
    echo "Marking & removing duplicates: " $SAMPLE_NAME
    java -Xmx4g -jar $PICARDDIR/MarkDuplicates.jar \
    INPUT=$PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.bam \
    OUTPUT=$PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.notsorted.bam \
    METRICS_FILE=$PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dups.txt \
    VALIDATION_STRINGENCY=LENIENT \
    TMP_DIR=$TMPDIR

    samtools sort $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.notsorted.bam $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup

    if [[ -s $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.notsorted.bam ]]
        then
        rm $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.notsorted.bam
    fi
else
    echo "Marking & removing duplicates: " $SAMPLE_NAME
    java -Xmx4g -jar $PICARDDIR/MarkDuplicates.jar \
    INPUT=$PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.bam \
    OUTPUT=$PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.bam \
    METRICS_FILE=$PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dups.txt \
    VALIDATION_STRINGENCY=LENIENT \
    TMP_DIR=$TMPDIR

    samtools sort $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.notsorted.bam $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup
    
    if [[ -s $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.notsorted.bam ]]
        then
        rm $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.notsorted.bam
    fi
fi

# Preseq
if [[ $SAMPLE_NAME == *ATAC-seq* ]] || [[ $SAMPLE_NAME == *_CM_* ]]
    then
    echo "QC: " $SAMPLE_NAME
    $PRESEQ c_curve -B $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam -o $PROJECTDIR/mapped/$SAMPLE_NAME_qc_c_curve.txt
    $PRESEQ lc_extrap -e 1e8 -s 2e6 -B $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam -o $PROJECTDIR/mapped/$SAMPLE_NAME_qc_lc_extrap.txt
else
    echo "QC: " $SAMPLE_NAME
    $PRESEQ c_curve -B $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.bam -o $PROJECTDIR/mapped/$SAMPLE_NAME_qc_c_curve.txt
    $PRESEQ lc_extrap -e 1e8 -s 2e6 -B $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.bam -o $PROJECTDIR/mapped/$SAMPLE_NAME_qc_lc_extrap.txt
fi


date