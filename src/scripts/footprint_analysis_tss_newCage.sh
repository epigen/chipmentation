#!/bin/bash
#SBATCH --partition=develop
#SBATCH --ntasks=1
#SBATCH --time=12:00:00

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64000
#SBATCH --nodes=1

#SBATCH --job-name=footprintingTSSs
#SBATCH --output=/home/arendeiro/logs/footprintingTSSs_%j.out

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
CAGE=/home/arendeiro/reference/Homo_sapiens/DPIcluster_hg19_20120116.permissive_set.TATA_CpG_annotated.onlyTATA.K562.bed
CHIP_SAMPLE=H3K4me3_K562_500k_ChIP
DNase_SAMPLE=DNase_UWashington_K562_mergedReplicates
TATA=/home/arendeiro/reference/Homo_sapiens/tata-promoters.bed


### Start work


# Make cage annotation
#grep -v '#' /fhgfs/groups/lab_bock/shared/data/posfiles/DPIcluster_hg19_20120116.permissive_set.TATA_CpG_annotated.osc | \
#tail -n +2 | \
#grep -v 'TATA-less' \
#> /home/arendeiro/reference/Homo_sapiens/DPIcluster_hg19_20120116.permissive_set.TATA_CpG_annotated.onlyTATA.bed

#bedtools intersect -wa \
#-a /home/arendeiro/reference/Homo_sapiens/DPIcluster_hg19_20120116.permissive_set.TATA_CpG_annotated.onlyTATA.bed \
#-b /fhgfs/groups/lab_bock/shared/data/posfiles/wgEncodeRikenCageK562CellPapTssHmm.bedRnaElements \
#> /home/arendeiro/reference/Homo_sapiens/DPIcluster_hg19_20120116.permissive_set.TATA_CpG_annotated.onlyTATA.K562.bed


##### DOnt do now
# Filter CAGE peaks by score
#awk -v OFS='\t' '$8 >= 3' /home/arendeiro/reference/Homo_sapiens/cage.txt | \
# Get central position
#awk -v OFS='\t' '{a=int(($3+$4)/2+0.5); $3=a; $4=a+1;print}' | \
#cut -f 2,3,4,5,6,7,8,9 > /home/arendeiro/reference/Homo_sapiens/cage_K562_cell_tss_filtered3.bed

# Get cage peaks in H3K4me3 peaks
bedtools intersect -wa \
-a $CAGE \
-b $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks/${SAMPLE_NAME}_peaks.narrowPeak \
> $PROJECTDIR/cageTSSs/${SAMPLE_NAME}_cageTSSs_in_peaks.bed

# Get 4kb window around cage peaks
bedtools slop -b 2000 \
-i $PROJECTDIR/cageTSSs/${SAMPLE_NAME}_cageTSSs_in_peaks.bed \
-g $GENOMESIZE \
> $PROJECTDIR/cageTSSs/${SAMPLE_NAME}_cageTSSs_in_peaks.2kbSlop.bed

# Filter TATA-containing promoters
#bedtools intersect -wa \
#-a $PROJECTDIR/cageTSSs/${SAMPLE_NAME}_cageTSSs_in_peaks.2kbSlop.bed \
#-b $TATA \
#> $PROJECTDIR/cageTSSs/${SAMPLE_NAME}_cageTSSs_in_peaks.2kbSlop.TATA.bed


# get only 5' position of reads
echo "Getting 5' read positions for sample: " $SAMPLE_NAME
bedtools bamtobed -i $PROJECTDIR/mapped/merged/${SAMPLE_NAME}.bam | \
python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
> $PROJECTDIR/bed/${SAMPLE_NAME}.5prime.bed
echo "Getting 5' read positions for sample: " $CHIP_SAMPLE
bedtools bamtobed -i $PROJECTDIR/mapped/merged/${CHIP_SAMPLE}.bam | \
python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
> $PROJECTDIR/bed/${CHIP_SAMPLE}.5prime.bed
echo "Getting 5' read positions for sample: " $CONTROL_NAME
bedtools bamtobed -i $PROJECTDIR/mapped/merged/${CONTROL_NAME}.bam | \
python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
> $PROJECTDIR/bed/${CONTROL_NAME}.5prime.bed
echo "Getting 5' read positions for DNase sample."
bedtools bamtobed -i $PROJECTDIR/mapped/$DNase_SAMPLE.bam | \
python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
> $PROJECTDIR/mapped/$DNase_SAMPLE.5prime.bed

# get read count of H3K4me3 in windows at single basepair resolution
## Chipmentation signal
echo "Getting read counts for sample: " $SAMPLE_NAME
bedtools coverage -d \
-a $PROJECTDIR/bed/${SAMPLE_NAME}.5prime.bed \
-b $PROJECTDIR/cageTSSs/${SAMPLE_NAME}_cageTSSs_in_peaks.2kbSlop.bed \
> $PROJECTDIR/bed/${SAMPLE_NAME}_tss_coverage.bed
## ChIP signal
echo "Getting read counts for sample: " $CHIP_SAMPLE
bedtools coverage -d \
-a $PROJECTDIR/bed/${CHIP_SAMPLE}.5prime.bed \
-b $PROJECTDIR/cageTSSs/${SAMPLE_NAME}_cageTSSs_in_peaks.2kbSlop.bed \
> $PROJECTDIR/bed/${SAMPLE_NAME}_tss_coverage_ChIP.bed
## IgG
echo "Getting read counts for sample: " $CONTROL_NAME
bedtools coverage -d \
-a $PROJECTDIR/bed/${CONTROL_NAME}.5prime.bed \
-b $PROJECTDIR/cageTSSs/${SAMPLE_NAME}_cageTSSs_in_peaks.2kbSlop.bed \
> $PROJECTDIR/bed/${SAMPLE_NAME}_tss_coverage_IgG.bed
## DNase
echo "Getting read counts for sample: " $DNase_SAMPLE
bedtools coverage -d \
-a $PROJECTDIR/mapped/${DNase_SAMPLE}.5prime.bed \
-b $PROJECTDIR/cageTSSs/${SAMPLE_NAME}_cageTSSs_in_peaks.2kbSlop.bed \
> $PROJECTDIR/bed/${SAMPLE_NAME}_tss_coverage_DNase.bed


# Get single basepair position within TSSs for nucl compos. coverage and conservation
echo "Getting single basepairs for sample:" $SAMPLE_NAME
cat $PROJECTDIR/bed/${SAMPLE_NAME}_tss_coverage.bed | \
python /home/arendeiro/projects/chipmentation/src/lib/get_single_bp_for_tss.py \
> $PROJECTDIR/bed/${SAMPLE_NAME}_tss_coverage_sbp.bed

# get nucleotide composition in windows
## center at single nucleotide resolution
## get 12bp windows around each nucleotide position in the 4kb
## see nucleotide composition in those 12bp windows
echo "Getting nucleotide composition for sample: " $SAMPLE_NAME
awk -v OFS='\t' '{print $1, $2 + $5, $2 + $5 + 1, $4}' $PROJECTDIR/bed/${SAMPLE_NAME}_tss_coverage_sbp.bed | \
bedtools slop -b 6 -i stdin -g $GENOMESIZE | \
bedtools nuc -fi $GENOMEREF -bed stdin \
> $PROJECTDIR/bed/${SAMPLE_NAME}_tss_nucleotide_compos.bed

# get conservation score
# independently by chromossome

echo "Getting Conservation scores for sample: " $SAMPLE_NAME
for chr in `ls $CONSERVATION/*.bed`
do 
    CHR=`basename $chr`
    echo "... chromossome " $CHR
    bedtools intersect -wa -wb -a $PROJECTDIR/bed/${SAMPLE_NAME}_tss_coverage_sbp.bed -b $chr \
    > $TMPDIR/${SAMPLE_NAME}_tss_conservation.${CHR}
done
# concatenate and sort
echo "Concatenation conservation scores for sample: " $SAMPLE_NAME
cat $TMPDIR/${SAMPLE_NAME}_tss_conservation.* > $PROJECTDIR/bed/${SAMPLE_NAME}_tss_conservation.bed

# For each position in each peak put all toghether (paste relevant columns from each file, remove headers when existing)
#cut -f 5,6,7,8,9,10 $PROJECTDIR/bed/${SAMPLE_NAME}_tss_nucleotide_compos.bed | \
#tail -n +2 | paste $PROJECTDIR/bed/${SAMPLE_NAME}_tss_coverage.bed - > $TMPDIR/tmp
#cut -f 11 $PROJECTDIR/bed/${SAMPLE_NAME}_tss_conservation.bed | paste $TMPDIR/tmp - \
#> $PROJECTDIR/bed/${SAMPLE_NAME}_tss.2kb_coverage.bed

# get TSSs with CpG islands

date