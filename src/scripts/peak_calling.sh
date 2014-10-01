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

### Specify paths
PROJECTDIR=/home/arendeiro/data/human/chipmentation
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19
PICARDDIR=/cm/shared/apps/picard-tools/1.118
PRESEQ=/home/arendeiro/.local/software/preseq-0.1.0/preseq
MACS2="python2.7 /home/arendeiro/.local/software/bin/macs2"
HOMERDIR=/home/arendeiro/.local/software/homer-4.6/bin

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
        --broad -g hs -n $SAMPLE_NAME --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2
    else
        echo "Calling normal peaks for sample " $SAMPLE_NAME
        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.shifted.dup.bam \
        -g hs -n $SAMPLE_NAME --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2
    fi
else
    if [[ $SAMPLE_NAME == *H3K27me3* ]]
        then
        echo "Calling broad peaks for sample " $SAMPLE_NAME
        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.dup.bam \
        --broad -g hs -n $SAMPLE_NAME --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2
    else
        echo "Calling normal peaks for sample " $SAMPLE_NAME
        $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.bam \
        -c $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.dup.bam \
        -g hs -n $SAMPLE_NAME --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2
    fi
fi


# HOMER
# make tag dirs
mkdir -p $PROJECTDIR/homer/${SAMPLE_NAME}_homer
if [[ $SAMPLE_NAME == *ATAC-seq* ]] || [[ $SAMPLE_NAME == *_CM_* ]]
    then
    $HOMERDIR/makeTagDirectory $PROJECTDIR/homer/${SAMPLE_NAME}_homer $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam
    $HOMERDIR/makeTagDirectory $PROJECTDIR/homer/${CONTROL_NAME}_homer $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.shifted.dup.bam
    # make UCSC tracks
    $HOMERDIR/makeUCSCfile $PROJECTDIR/homer/${SAMPLE_NAME}_homer -name $SAMPLE_NAME -o auto
    $HOMERDIR/makeUCSCfile $PROJECTDIR/homer/${CONTROL_NAME}_homer -name $CONTROL_NAME -o auto
else
    $HOMERDIR/makeTagDirectory $PROJECTDIR/homer/${SAMPLE_NAME}_homer $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.bam
    $HOMERDIR/makeTagDirectory $PROJECTDIR/homer/${CONTROL_NAME}_homer $PROJECTDIR/mapped/$CONTROL_NAME.trimmed.bowtie2.sorted.dup.bam
    # make UCSC tracks
    $HOMERDIR/makeUCSCfile $PROJECTDIR/homer/${SAMPLE_NAME}_homer -name $SAMPLE_NAME -o auto
    $HOMERDIR/makeUCSCfile $PROJECTDIR/homer/${CONTROL_NAME}_homer -name $CONTROL_NAME -o auto
fi

# find peaks
if [[ $SAMPLE_NAME == *H3* ]]
    then
    echo "Detected histone data"
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style histone -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer_std
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style histone -region -size 150 -minDist 370 -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer_region_150_370
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style histone -region -size 150 -minDist 1000 -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer_region_150_1000
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style histone -region -size 1000 -minDist 2500 -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer_region_1000_2500
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style histone -region -size 1000 -minDist 10000 -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer_region_1000_10000
else
    echo "Detected TF data"
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style factor -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer_std
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style factor -C 1 -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer_c1
    $HOMERDIR/findPeaks $PROJECTDIR/homer/${SAMPLE_NAME}_homer -style factor -C 0 -o auto -i $PROJECTDIR/homer/${CONTROL_NAME}_homer_c0
fi

date