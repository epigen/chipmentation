#!/bin/bash
#SBATCH --partition=develop
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

# Optional parameters
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --nodes=1

# *** setup environment ***
# load the required environmental modules that your script depends upon
module load FastQC/0.11.2
module load trimmomatic/0.32
module load samtools
module load bamtools/bamtools
module load bowtie/2.2.0
module load picard-tools/1.118
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
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19
GENOMESIZE=/fhgfs/prod/ngs_resources/genomes/hg19/hg19_chromLengths_sorted.txt
PICARDDIR=/cm/shared/apps/picard-tools/1.118
PRESEQ=/home/arendeiro/.local/software/preseq-0.1.0/preseq
HOMERDIR=/home/arendeiro/.local/software/homer-4.6/bin
MACS2="python2.7 /home/arendeiro/.local/software/bin/macs2"
BEDTOOLSDIR=/home/arendeiro/.local/software/bedtools2/bin

### Call Peaks
mkdir -p $PROJECTDIR/peaks

if [[ $SAMPLE_NAME != *H3K27me3* ]]
    then
    #$MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.bam \
    #-c $PROJECTDIR/mapped/$CONTROL_NAME.bam \
    #-g hs -n ${SAMPLE_NAME}_std --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_std

    $MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.bam \
    -c $PROJECTDIR/mapped/$CONTROL_NAME.bam \
    --bw 200 \
    -g hs -n ${SAMPLE_NAME}_bw200 --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_bw200

    #$MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.bam \
    #-c $PROJECTDIR/mapped/$CONTROL_NAME.bam \
    # --summits \
    #-g hs -n ${SAMPLE_NAME}_summits --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_summits

    awk '$9 >= 10' $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_bw200/${SAMPLE_NAME}_bw200_peaks.narrowPeak > $PROJECTDIR/peaks/${SAMPLE_NAME}.narrowPeak
    awk '{print $1, $2 + $10, $2 + $10 + 1, $4}' $PROJECTDIR/peaks/${SAMPLE_NAME}.narrowPeak > $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.bed
    $BEDTOOLSDIR/slopBed -b 2000 -i $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.bed -g $GENOMESIZE > $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.2kb.bed

else
    $MACS2 callpeak -B -t $PROJECTDIR/mapped/$SAMPLE_NAME.bam \
    -c $PROJECTDIR/mapped/$CONTROL_NAME.bam \
    --broad --nomodel --extsize 200 \
    -g hs -n ${SAMPLE_NAME}_std --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_broad_nomodel_extsize200

    $MACS2 callpeak -B -t $PROJECTDIR/mapped/$SAMPLE_NAME.bam \
    -c $PROJECTDIR/mapped/$CONTROL_NAME.bam \
    --broad --nomodel --extsize 200 --pvalue 1e-3 \
    -g hs -n ${SAMPLE_NAME}_std --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3

    #$MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.bam \
    #--format BED \
    #-c $PROJECTDIR/mapped/$CONTROL_NAME.bam \
    #--broad --nomodel --extsize 200 \ # provide fragment length if MACS has problems figuring it out
    #-g hs -n ${SAMPLE_NAME}_std --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_std

    #$MACS2 callpeak -t $PROJECTDIR/mapped/$SAMPLE_NAME.bam \
    #-c $PROJECTDIR/mapped/$CONTROL_NAME.bam \
    # --broad --bw 200 \
    #-g hs -n ${SAMPLE_NAME}_bw200 --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_bw200

    #awk '$9 >= 10' $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_bw200/${SAMPLE_NAME}_bw200_peaks.broadPeak > $PROJECTDIR/peaks/${SAMPLE_NAME}.broadPeak
fi

date