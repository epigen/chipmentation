#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2000
#SBATCH --nodes=1

#SBATCH --job-name=mapping
#SBATCH --output=/home/arendeiro/logs/mapping_%j.out
#SBATCH --error=/home/arendeiro/logs/mapping_%j.err

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
HOMERDIR=/home/arendeiro/.local/software/homer-4.6/bin

### Start work on samples 
# Fastqc on raw data
mkdir -p $PROJECTDIR/raw/fastqc
fastqc -o $PROJECTDIR/raw/fastqc $SAMPLE_FILE

# Convert to fastq
echo "Converting to fastq: " $SAMPLE_NAME
mkdir -p $PROJECTDIR/raw/fastq
bamtools convert -in $SAMPLE_FILE -format fastq > $PROJECTDIR/raw/fastq/$SAMPLE_NAME.fastq

# Trimming 
# Remove adapters, trim reads when needed
echo "Trimming reads: " $SAMPLE_NAME
java -jar `which trimmomatic-0.32.jar` SE \
$PROJECTDIR/raw/fastq/$SAMPLE_NAME.fastq ${PROJECTDIR}/raw/${SAMPLE_NAME}.trimmed.fastq \
ILLUMINACLIP:/home/arendeiro/adapters/chipmentation.fa:1:40:15 \
LEADING:3 TRAILING:3 \
SLIDINGWINDOW:4:10 \
MINLEN:36

# Map with bowtie2 end-to-end sensitive settings
echo "Mapping: " $SAMPLE_NAME
mkdir -p $PROJECTDIR/mapped/
bowtie2 --very-sensitive -p 16 -x $GENOMEREF $PROJECTDIR/raw/$SAMPLE_NAME.trimmed.fastq | samtools view -S -b - | samtools sort - $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted

# Shift reads from tagmentation samples
if [[ $SAMPLE_NAME == *ATAC-seq* ]] || [[ $SAMPLE_NAME == *_CM_* ]]
    then
    echo "Shifting reads for sample:" $SAMPLE_NAME
    samtools view -h $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.bam | python2.7 projects/chipmentation/src/scripts/shift_reads.py | samtools view -S -b - | samtools sort - $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted
fi

# Mark duplicates
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

# HOMER
# make tag dirs
mkdir -p $PROJECTDIR/homer/${SAMPLE_NAME}_homer
if [[ $SAMPLE_NAME == *ATAC-seq* ]] || [[ $SAMPLE_NAME == *_CM_* ]]
    then
    echo "Making tag directory for sample: " $SAMPLE_NAME
    $HOMERDIR/makeTagDirectory $PROJECTDIR/homer/${SAMPLE_NAME}_homer $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.shifted.dup.bam
    # make UCSC tracks
    echo "Making tracks for sample: " $SAMPLE_NAME
    $HOMERDIR/makeUCSCfile $PROJECTDIR/homer/${SAMPLE_NAME}_homer -name $SAMPLE_NAME -o auto
else
    echo "Making tag directory for sample: " $SAMPLE_NAME
    $HOMERDIR/makeTagDirectory $PROJECTDIR/homer/${SAMPLE_NAME}_homer $PROJECTDIR/mapped/$SAMPLE_NAME.trimmed.bowtie2.sorted.dup.bam
    # make UCSC tracks
    echo "Making tracks for sample: " $SAMPLE_NAME
    $HOMERDIR/makeUCSCfile $PROJECTDIR/homer/${SAMPLE_NAME}_homer -name $SAMPLE_NAME -o auto
fi

date