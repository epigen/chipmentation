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
# SAMPLE_NAME=PU1_K562_10mio_CM
# CONTROL_NAME=IgG_K562_10mio_CM
SAMPLE_NAME=$1
CONTROL_NAME=$2
if [[ $SAMPLE_NAME == *CTCF* ]]; then
    ENCODE_SAMPLE=wgEncodeHaibTfbsK562CtcfcPcr1xAln
elif [[ $SAMPLE_NAME == *PU1* ]]; then
    ENCODE_SAMPLE=wgEncodeHaibTfbsK562Pu1Pcr1xAln
fi

### Specify paths
SAMPLES_FILE=/home/arendeiro/projects/chipmentation/samples_peaks.txt
PROJECTDIR=/home/arendeiro/data/human/chipmentation
GENOMESIZE=/fhgfs/prod/ngs_resources/genomes/hg19/hg19_chromLengths_sorted.txt
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19.fa
CONSERVATION=/home/arendeiro/reference/Homo_sapiens/phyloP/placentalMammals

CHIP_SAMPLE=${SAMPLE_NAME/CM/ChIP}
DNase_SAMPLE=DNase_UWashington_K562_mergedReplicates

### Start work on samples 

# center on motifs
annotatePeaks.pl $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks/${SAMPLE_NAME}_peaks.narrowPeak hg19 -size 4000 \
-center $PROJECTDIR/motifs/$SAMPLE_NAME/homerResults/motif1.motif | \
awk -v OFS='\t' '{print $2, $3, $4, $1, $6, $5}' | 
python /home/arendeiro/projects/chipmentation/src/lib/fix_bedfile_genome_boundaries.py | \
sortBed > $PROJECTDIR/bed/${SAMPLE_NAME}.motifStrand.bed

# get only 5' position of reads
echo "Getting 5' read positions for sample: " $SAMPLE_NAME
sbatch /home/arendeiro/jobScripts/bamTo5primeBed.job.sh \
$PROJECTDIR/mapped/merged/${SAMPLE_NAME}.bam \
$PROJECTDIR/mapped/merged/${SAMPLE_NAME}.5prime.bed
echo "Getting 5' read positions for sample: " $CONTROL_NAME
sbatch /home/arendeiro/jobScripts/bamTo5primeBed.job.sh \
$PROJECTDIR/mapped/merged/${CONTROL_NAME}.bam \
$PROJECTDIR/mapped/merged/${CONTROL_NAME}.5prime.bed
echo "Getting 5' read positions for sample: " $CHIP_SAMPLE
sbatch /home/arendeiro/jobScripts/bamTo5primeBed.job.sh \
$PROJECTDIR/mapped/merged/${CHIP_SAMPLE}.bam \
$PROJECTDIR/mapped/merged/${CHIP_SAMPLE}.5prime.bed
echo "Getting 5' read positions for DNase sample."
sbatch /home/arendeiro/jobScripts/bamTo5primeBed.job.sh \
$PROJECTDIR/mapped/$DNase_SAMPLE.bam \
$PROJECTDIR/mapped/merged/$DNase_SAMPLE.5prime.bed
echo "Getting 5' read positions for Encode sample."
sbatch /home/arendeiro/jobScripts/bamTo5primeBed.job.sh \
/home/arendeiro/data/human/encode/chip-seq/$ENCODE_SAMPLE.bam \
$PROJECTDIR/mapped/merged/$ENCODE_SAMPLE.5prime.bed


# get read count of PU1 in windows
echo "Getting read counts for sample: " $SAMPLE_NAME
sbatch /home/arendeiro/jobScripts/bedToolsCoverage.job.sh \
$PROJECTDIR/mapped/merged/${SAMPLE_NAME}.5prime.bed \
$PROJECTDIR/bed/${SAMPLE_NAME}.motifStrand.bed \
$PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_CM.bed

# get read count of IgG in windows
echo "Getting read counts for sample: " $CONTROL_NAME
sbatch /home/arendeiro/jobScripts/bedToolsCoverage.job.sh \
$PROJECTDIR/mapped/merged/${CONTROL_NAME}.5prime.bed \
$PROJECTDIR/bed/${SAMPLE_NAME}.motifStrand.bed \
$PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_IgG.bed
# get read count of ChIP in windows
echo "Getting read counts for sample: " $CHIP_SAMPLE
sbatch /home/arendeiro/jobScripts/bedToolsCoverage.job.sh \
$PROJECTDIR/mapped/merged/${CHIP_SAMPLE}.5prime.bed \
$PROJECTDIR/bed/${SAMPLE_NAME}.motifStrand.bed \
$PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_ChIP.bed
# get read count of DNase in windows
echo "Getting read counts for sample: " $DNase_SAMPLE
sbatch /home/arendeiro/jobScripts/bedToolsCoverage.job.sh \
$PROJECTDIR/mapped/merged/${DNase_SAMPLE}.5prime.bed \
$PROJECTDIR/bed/${SAMPLE_NAME}.motifStrand.bed \
$PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_DNase.bed
# get read count of Encode sample in windows
echo "Getting read counts for sample: " $ENCODE_SAMPLE
sbatch /home/arendeiro/jobScripts/bedToolsCoverage.job.sh \
$PROJECTDIR/mapped/merged/$ENCODE_SAMPLE.5prime.bed \
$PROJECTDIR/bed/${SAMPLE_NAME}.motifStrand.bed \
$PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_Encode.bed

# get nucleotide composition in windows
## center at single nucleotide resolution
## get 12bp windows around each nucleotide position in the 2kb windows
## see nucleotide composition in those 12bp windows
# echo "Getting nucleotide composition for sample: " $SAMPLE_NAME
# awk -v OFS='\t' '{print $1, $2 + $5, $2 + $5 + 1, $4}' $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage.bed | \
# bedtools slop -b 6 -i stdin -g $GENOMESIZE | \
# bedtools nuc -fi $GENOMEREF -bed stdin > $PROJECTDIR/bed/${SAMPLE_NAME}_peak_nucleotide_compos.bed

# # get conservation score
# # independently by chromossome
# echo "Getting single basepairs for sample:" $SAMPLE_NAME
# cat $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage.bed | python /home/arendeiro/projects/chipmentation/src/lib/get_single_bp.py > $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_sbp.bed

# echo "Getting Conservation scores for sample: " $SAMPLE_NAME
# for chr in `ls $CONSERVATION/*.bed`
# do 
#     CHR=`basename $chr`
#     echo "... chromossome " $CHR
#     bedtools intersect -wa -wb -a $PROJECTDIR/bed/${SAMPLE_NAME}.motif.bed -b $chr \
#     > $TMPDIR/${SAMPLE_NAME}_wholepeak_conservation.${CHR}
# done

# # concatenate
# echo "Concatenation conservation scores for sample: " $SAMPLE_NAME
# cat $TMPDIR/${SAMPLE_NAME}_wholepeak_conservation.* > $PROJECTDIR/bed/${SAMPLE_NAME}_wholepeak_conservation.bed


# Parse coverage into a csv 
for TECH in CM ChIP IgG DNase Encode
do
    sbatch /home/arendeiro/jobScripts/parseCoverage.job.sh \
    $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage_${TECH}.bed
done

#### Co-binding prediction (GATA on PU.1)
if [[ $SAMPLE_NAME == *PU1* ]]; then
    # center on motifs
    echo "Centering on second most abundant motif in sample: " $SAMPLE_NAME
    annotatePeaks.pl $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks/${SAMPLE_NAME}_peaks.narrowPeak hg19 -size 4000 \
    -center $PROJECTDIR/motifs/$SAMPLE_NAME/homerResults/motif2.motif | \
    awk -v OFS='\t' '{print $2, $3, $4, $1, $6, $5}' | 
    sortBed \
    > $PROJECTDIR/bed/${SAMPLE_NAME}.GATA.motifStrand.bed

    # get read count of PU1 CM around GATA motif
    echo "Getting read counts for sample: " $SAMPLE_NAME
    bedtools coverage -d \
    -a $PROJECTDIR/mapped/merged/${SAMPLE_NAME}.5prime.bed \
    -b $PROJECTDIR/bed/${SAMPLE_NAME}.GATA.motifStrand.bed \
    > $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage.GATAmotif.bed

    # parse to csv
    sbatch /home/arendeiro/jobScripts/parseCoverage.job.sh \
    $PROJECTDIR/bed/${SAMPLE_NAME}_peak_coverage.GATAmotif.bed
fi


date
