#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --nodes=1

#SBATCH --job-name=macs2_callPeaks
#SBATCH --output=/home/arendeiro/logs/macs2_callPeaks_%j.out
#SBATCH --error=/home/arendeiro/logs/macs2_callPeaks_%j.err

# *** setup environment ***
# load the required environmental modules that your script depends upon
module load python/2.7.6

# set the temporary folder to your scratch location (avoid using node local /tmp)
export SCRATCH_PATH=/fhgfs/scratch/users/arendeiro
export TMPDIR=$SCRATCH_PATH/tmp

# *** run the job ***
hostname
date

### Get sample info from arguments
SAMPLE_NAME=$1
CONTROL_NAME=$2

### Specify paths
PROJECTDIR=/home/arendeiro/data/human/chipmentation
MACS2="python2.7 .local/software/bin/macs2"

### Call Peaks
mkdir -p $PROJECTDIR/peaks/

# MACS
if [[ $SAMPLE_NAME == *ATAC-seq* ]] || [[ $SAMPLE_NAME == *_CM_* ]]
    then
    if [[ $SAMPLE_NAME == *H3K27me3* ]]
        then
        echo "Calling broad peaks for sample " $SAMPLE_NAME
        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        --broad -g hs \
        -n ${SAMPLE_NAME}_std -B --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_std

        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        --broad -g hs --bw 200 \
        -n ${SAMPLE_NAME}_bw200 -B --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_bw200

        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        --broad -g hs --keep-dup \
        -n ${SAMPLE_NAME}_keepdup -B --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_keepdup

        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        --broad -g hs --call-summits \
        -n ${SAMPLE_NAME}_summits -B --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_summits
    else
        echo "Calling normal peaks for sample " $SAMPLE_NAME
        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        -g hs \
        -n ${SAMPLE_NAME}_std -B --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_std

        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        -g hs --bw 200 \
        -n ${SAMPLE_NAME}_bw200 -B --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_bw200

        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        -g hs --keep-dup \
        -n ${SAMPLE_NAME}_keepdup -B --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_keepdup

        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        -g hs --call-summits \
        -n ${SAMPLE_NAME}_summits -B --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_summits
    fi
else
    if [[ $SAMPLE_NAME == *H3K27me3* ]]
        then
        echo "Calling broad peaks for sample " $SAMPLE_NAME
        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.dup.bam \
        --broad -g hs \
        -n ${SAMPLE_NAME}_std -B --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_std

        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.dup.bam \
        --broad -g hs --bw 200 \
        -n ${SAMPLE_NAME}_bw200 -B --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_bw200

        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.dup.bam \
        --broad -g hs --keep-dup \
        -n ${SAMPLE_NAME}_keepdup -B --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_keepdup

        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.dup.bam \
        --broad -g hs --call-summits \
        -n ${SAMPLE_NAME}_summits -B --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_summits
    else
        echo "Calling normal peaks for sample " $SAMPLE_NAME
        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.dup.bam \
        -g hs \
        -n ${SAMPLE_NAME}_std -B --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_std

        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.dup.bam \
        -g hs --bw 200 \
        -n ${SAMPLE_NAME}_bw200 -B --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_bw200

        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.dup.bam \
        -g hs --keep-dup \
        -n ${SAMPLE_NAME}_keepdup -B --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_keepdup

        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.dup.bam \
        -g hs --call-summits \
        -n ${SAMPLE_NAME}_summits -B --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_summits
    fi
fi

date