#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2000
#SBATCH --nodes=1

#SBATCH --job-name=mapping_%j
#SBATCH --output=/home/arendeiro/logs/mapping_%j.out
#SBATCH --error=/home/arendeiro/logs/mapping_%j.err
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
SAMPLE_FILE_1=$2
SAMPLE_FILE_2=$3
### Specify paths
PROJECTDIR=/home/arendeiro/data/human/chipmentation
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19

# Map with bowtie2 end-to-end sensitive settings
mkdir -p $PROJECTDIR/mapped/
bowtie2 --very-fast-local -p 16 -x $GENOMEREF -1 $SAMPLE_FILE_1 -2 $SAMPLE_FILE_2 | samtools view -S -b - | samtools sort - $PROJECTDIR/mapped/$SAMPLE_NAME.bowtie2.vfl.sorted
bowtie2 --fast-local -p 16 -x $GENOMEREF -1 $SAMPLE_FILE_1 -2 $SAMPLE_FILE_2 | samtools view -S -b - | samtools sort - $PROJECTDIR/mapped/$SAMPLE_NAME.bowtie2.fl.sorted
bowtie2 --sensitive-local -p 16 -x $GENOMEREF -1 $SAMPLE_FILE_1 -2 $SAMPLE_FILE_2 | samtools view -S -b - | samtools sort - $PROJECTDIR/mapped/$SAMPLE_NAME.bowtie2.sl.sorted
bowtie2 --very-sensitive-local -p 16 -x $GENOMEREF -1 $SAMPLE_FILE_1 -2 $SAMPLE_FILE_2 | samtools view -S -b - | samtools sort - $PROJECTDIR/mapped/$SAMPLE_NAME.bowtie2.vsl.sorted

bowtie2 --very-fast -p 16 -x $GENOMEREF -1 $SAMPLE_FILE_1 -2 $SAMPLE_FILE_2 | samtools view -S -b - | samtools sort - $PROJECTDIR/mapped/$SAMPLE_NAME.bowtie2.vf.sorted
bowtie2 --fast -p 16 -x $GENOMEREF -1 $SAMPLE_FILE_1 -2 $SAMPLE_FILE_2 | samtools view -S -b - | samtools sort - $PROJECTDIR/mapped/$SAMPLE_NAME.bowtie2.f.sorted
bowtie2 --sensitive -p 16 -x $GENOMEREF -1 $SAMPLE_FILE_1 -2 $SAMPLE_FILE_2 | samtools view -S -b - | samtools sort - $PROJECTDIR/mapped/$SAMPLE_NAME.bowtie2.s.sorted
bowtie2 --very-sensitive -p 16 -x $GENOMEREF -1 $SAMPLE_FILE_1 -2 $SAMPLE_FILE_2 | samtools view -S -b - | samtools sort - $PROJECTDIR/mapped/$SAMPLE_NAME.bowtie2.vs.sorted

date