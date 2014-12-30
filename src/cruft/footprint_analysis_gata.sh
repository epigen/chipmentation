#!/bin/bash
#SBATCH --partition=longq
#SBATCH --ntasks=1
#SBATCH --time=72:00:00

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32000
#SBATCH --nodes=1

#SBATCH --job-name=footprintingPeaks
#SBATCH --output=/home/arendeiro/logs/footprintingPeaks_%j.out

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
CONTROL_NAME=$2

### Specify paths
SAMPLES_FILE=/home/arendeiro/projects/chipmentation/samples_peaks.txt
PROJECTDIR=/home/arendeiro/data/human/chipmentation
GENOMESIZE=/fhgfs/prod/ngs_resources/genomes/hg19/hg19_chromLengths_sorted.txt
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19.fa
CONSERVATION=/home/arendeiro/reference/Homo_sapiens/phyloP/placentalMammals
CAGE=/home/arendeiro/reference/Homo_sapiens/cage_K562_cell_tss_f10.bed

CHIP_SAMPLE=${SAMPLE_NAME/CM/ChIP}
DNase_SAMPLE=DNase_UWashington_K562_mergedReplicates

### Start work on samples 

# center on motifs
annotatePeaks.pl $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks/${SAMPLE_NAME}_peaks.narrowPeak hg19 -size 4000 \
-center $PROJECTDIR/motifs/$SAMPLE_NAME/homerResults/motif2.motif | \
awk -v OFS='\t' '{print $2, $3, $4, $1}' | \
sortBed > $PROJECTDIR/bed/${SAMPLE_NAME}.GATAmotif.bed

# get read count of PU1 in windows
echo "Getting read counts for sample: " $SAMPLE_NAME
bedtools coverage -d \
-a $PROJECTDIR/mapped/merged/${SAMPLE_NAME}.5prime.bed \
-b $PROJECTDIR/bed/${SAMPLE_NAME}.GATAmotif.bed \
> $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage.GATAmotif.bed

date
