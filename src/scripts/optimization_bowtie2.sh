#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

# Optional parameters
#SBATCH --cpus-per-task=16
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

### Start work on samples 
# Fastqc on raw data
mkdir -p $PROJECTDIR/raw/fastqc
fastqc -o $PROJECTDIR/raw/fastqc $SAMPLE_FILE

# Convert to fastq
mkdir -p $PROJECTDIR/raw/fastq
bamtools convert -in $SAMPLE_FILE -format fastq > $PROJECTDIR/raw/fastq/$SAMPLE_NAME.fastq

# Trimming 
# prepare folder for trimmed data
mkdir -p $PROJECTDIR/raw/trimmed

# Remove adapters, trim reads when needed
java -jar `which trimmomatic-0.32.jar` SE \
$PROJECTDIR/raw/fastq/$SAMPLE_NAME.fastq ${PROJECTDIR}/raw/trimmed/${SAMPLE_NAME}.trimmed.q10.fastq \
ILLUMINACLIP:/home/arendeiro/adapters/TruSeq3-SE.fa:1:40:15 \
LEADING:3 TRAILING:3 \
SLIDINGWINDOW:4:10 \
MINLEN:36

# Fastqc on trimmed data
mkdir -p $PROJECTDIR/raw/trimmed/fastqc
fastqc -o $PROJECTDIR/raw/trimmed/fastqc $PROJECTDIR/raw/trimmed/$SAMPLE_NAME.trimmed.q10.fastq

# Map with bowtie2 end-to-end sensitive settings
mkdir -p $PROJECTDIR/trimmed/mapped/

bowtie2 --very-fast-local -p 16 -x $GENOMEREF $PROJECTDIR/raw/trimmed/$SAMPLE_NAME.trimmed.q10.fastq | samtools view -S -b - | samtools sort - $PROJECTDIR/trimmed/mapped/$SAMPLE_NAME.trimmed.q10.bowtie2.vfl.sorted.bam
bowtie2 --fast-local -p 16 -x $GENOMEREF $PROJECTDIR/raw/trimmed/$SAMPLE_NAME.trimmed.q10.fastq | samtools view -S -b - | samtools sort - $PROJECTDIR/trimmed/mapped/$SAMPLE_NAME.trimmed.q10.bowtie2.fl.sorted.bam
bowtie2 --sensitive-local -p 16 -x $GENOMEREF $PROJECTDIR/raw/trimmed/$SAMPLE_NAME.trimmed.q10.fastq | samtools view -S -b - | samtools sort - $PROJECTDIR/trimmed/mapped/$SAMPLE_NAME.trimmed.q10.bowtie2.sl.sorted.bam
bowtie2 --very-sensitive-local -p 16 -x $GENOMEREF $PROJECTDIR/raw/trimmed/$SAMPLE_NAME.trimmed.q10.fastq | samtools view -S -b - | samtools sort - $PROJECTDIR/trimmed/mapped/$SAMPLE_NAME.trimmed.q10.bowtie2.vsl.sorted.bam

bowtie2 --very-fast -p 16 -x $GENOMEREF $PROJECTDIR/raw/trimmed/$SAMPLE_NAME.trimmed.q10.fastq | samtools view -S -b - | samtools sort - $PROJECTDIR/trimmed/mapped/$SAMPLE_NAME.trimmed.q10.bowtie2.vf.sorted.bam
bowtie2 --fast -p 16 -x $GENOMEREF $PROJECTDIR/raw/trimmed/$SAMPLE_NAME.trimmed.q10.fastq | samtools view -S -b - | samtools sort - $PROJECTDIR/trimmed/mapped/$SAMPLE_NAME.trimmed.q10.bowtie2.f.sorted.bam
bowtie2 --sensitive -p 16 -x $GENOMEREF $PROJECTDIR/raw/trimmed/$SAMPLE_NAME.trimmed.q10.fastq | samtools view -S -b - | samtools sort - $PROJECTDIR/trimmed/mapped/$SAMPLE_NAME.trimmed.q10.bowtie2.s.sorted.bam
bowtie2 --very-sensitive -p 16 -x $GENOMEREF $PROJECTDIR/raw/trimmed/$SAMPLE_NAME.trimmed.q10.fastq | samtools view -S -b - | samtools sort - $PROJECTDIR/trimmed/mapped/$SAMPLE_NAME.trimmed.q10.bowtie2.vs.sorted.bam

date