#!/bin/bash
#SBATCH --partition=develop
#SBATCH --ntasks=1
#SBATCH --time=12:00:00

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=64000
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
BEDTOOLSDIR=/home/arendeiro/.local/software/bedtools2/bin
HOMERDIR=/home/arendeiro/.local/software/homer-4.6/bin
CONSERVATION=/home/arendeiro/reference/Homo_sapiens/phyloP/placentalMammals
CAGE=/home/arendeiro/reference/Homo_sapiens/cage_K562_cell_tss_f10.bed

CHIP_SAMPLE=PU1_K562_10mio_ChIP
DNase_SAMPLE=DNase_UWashington_K562_mergedReplicates

### Start work on samples 

# center on motifs
annotatePeaks.pl $PROJECTDIR/peaks/${SAMPLE_NAME}.narrowPeak hg19 -size 4000 -center $PROJECTDIR/motifs/$SAMPLE_NAME/homerResults/motif1.motif | \
awk -v OFS='\t' '{print $2, $3, $4, $1}' | \
sortBed > $PROJECTDIR/bed/${SAMPLE_NAME}.motif.bed

# get summits
#awk -v OFS='\t' '{print $1, $2 + $10, $2 + $10 + 1, $4}' $PROJECTDIR/peaks/${SAMPLE_NAME}.narrowPeak > $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.bed
# get 4kb window around summits
#bedtools slopBed -b 2000 -i $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.bed -g $GENOMESIZE > $PROJECTDIR/peaks/${SAMPLE_NAME}.summits.2kb.bed

# get only 5' position of reads
echo "Getting 5' read positions for sample: " $SAMPLE_NAME
bedtools bamtobed -i $PROJECTDIR/mapped/merged/${SAMPLE_NAME}.bam | \
python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
> $PROJECTDIR/mapped/${SAMPLE_NAME}.5prime.bed

echo "Getting 5' read positions for sample: " $CONTROL_NAME
bedtools bamtobed -i $PROJECTDIR/mapped/merged/${CONTROL_NAME}.bam | \
python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
> $PROJECTDIR/mapped/${CONTROL_NAME}.5prime.bed

echo "Getting 5' read positions for sample: " $CHIP_SAMPLE
bedtools bamtobed -i $PROJECTDIR/mapped/merged/${CHIP_SAMPLE}.bam | \
python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
> $PROJECTDIR/mapped/${CHIP_SAMPLE}.5prime.bed

echo "Getting 5' read positions for DNase sample."
bedtools bamtobed -i $PROJECTDIR/mapped/$DNase_SAMPLE.bam | \
python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
> $PROJECTDIR/mapped/$DNase_SAMPLE.5prime.bed

# get read count of PU1 CM in windows
echo "Getting read counts for sample: " $SAMPLE_NAME
bedtools coverage -d \
-a $PROJECTDIR/mapped/${SAMPLE_NAME}.5prime.bed \
-b $PROJECTDIR/bed/${SAMPLE_NAME}.motif.bed \
> $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage.bed

# get read count of ChIP in windows
echo "Getting read counts for sample: " $CHIP_SAMPLE
bedtools coverage -d \
-a $PROJECTDIR/mapped/${CHIP_SAMPLE}.5prime.bed \
-b $PROJECTDIR/bed/${SAMPLE_NAME}.motif.bed \
> $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_ChIP.bed

# get read count of IgG in windows
echo "Getting read counts for sample: " $CONTROL_NAME
bedtools coverage -d \
-a $PROJECTDIR/mapped/${CONTROL_NAME}.5prime.bed \
-b $PROJECTDIR/bed/${SAMPLE_NAME}.motif.bed \
> $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_IgG.bed

# get read count of DNase in windows
echo "Getting read counts for sample: " $DNase_SAMPLE
bedtools coverage -d \
-a $PROJECTDIR/mapped/${DNase_SAMPLE}.5prime.bed \
-b $PROJECTDIR/bed/${SAMPLE_NAME}.motif.bed \
> $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_DNase.bed

# get nucleotide composition in windows
## center at single nucleotide resolution
## get 12bp windows around each nucleotide position in the 2kb windows
## see nucleotide composition in those 12bp windows
echo "Getting nucleotide composition for sample: " $SAMPLE_NAME
awk -v OFS='\t' '{print $1, $2 + $5, $2 + $5 + 1, $4}' $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage.bed | \
bedtools slop -b 6 -i stdin -g $GENOMESIZE | \
bedtools nuc -fi $GENOMEREF -bed stdin > $PROJECTDIR/bed/${SAMPLE_NAME}_peak_nucleotide_compos.bed

# get conservation score
# independently by chromossome
echo "Getting single basepairs for sample:" $SAMPLE_NAME
cat $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage.bed | python /home/arendeiro/projects/chipmentation/src/lib/get_single_bp.py > $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_sbp.bed

echo "Getting Conservation scores for sample: " $SAMPLE_NAME
for chr in `ls $CONSERVATION/*.bed`
do 
    CHR=`basename $chr`
    echo "... chromossome " $CHR
    bedtools intersect -wa -wb -a $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_sbp.bed -b $chr \
    > $TMPDIR/${SAMPLE_NAME}_peak_conservation.${CHR}
done
# concatenate and sort
echo "Concatenation conservation scores for sample: " $SAMPLE_NAME
cat $TMPDIR/${SAMPLE_NAME}_peak_conservation.* > $PROJECTDIR/bed/${SAMPLE_NAME}_peak_conservation.bed

# For each position in each peak put all toghether (paste relevant columns from each file, remove headers when existing)
cut -f 6 $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_IgG.bed | paste $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage.bed - > $TMPDIR/tmp
cut -f 5,6,7,8,9,10 $PROJECTDIR/bed/${SAMPLE_NAME}_peak_nucleotide_compos.bed | tail -n +2 | paste $TMPDIR/tmp - > $TMPDIR/tmp2
cut -f 11 $PROJECTDIR/bed/${SAMPLE_NAME}_peak_conservation.bed | paste $TMPDIR/tmp2 - > $PROJECTDIR/bed/${SAMPLE_NAME}_peak.2kb_coverage.bed



