

# Get nucleossomes
ROOT=/media/afr/cemm-backup/chipmentation
SAMPLE=K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19
PEAKS=${ROOT}/data/peaks/${SAMPLE}/${SAMPLE}_peaks.narrowPeak
BAM=${ROOT}/data/mapped/${SAMPLE}.trimmed.bowtie2.shifted.dups.bam
GENOME=/media/afr/cemm-backup/reference/hg19/forBowtie2/hg19.fa
OUTPUT=${ROOT}/data/nucleoATAC/$SAMPLE
mkdir -p $OUTPUT

samtools index $BAM
nucleoatac run --bed $PEAKS --bam $BAM --fasta $GENOME --out $OUTPUT --cores 16


# 
ROOT=/media/afr/cemm-backup/chipmentation
SAMPLE=K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19
PEAKS=${ROOT}/data/peaks/${SAMPLE}/${SAMPLE}_peaks.narrowPeak
BAM=${ROOT}/data/mapped/${SAMPLE}.trimmed.bowtie2.shifted.dups.bam
GENOME=/media/afr/cemm-backup/reference/hg19/forBowtie2/hg19.fa
OUTPUT=${ROOT}/data/nucleoATAC/${SAMPLE}-occ
mkdir -p $OUTPUT

samtools index $BAM
nucleoatac occ --lower 50 --bed $PEAKS --bam $BAM --out $OUTPUT --write_peaks --cores 4

# VPLOTS
# pu1

ROOT=/media/afr/cemm-backup/chipmentation
SAMPLE=K562_10M_CM_PU1_nan_PE_1_1_hg19
PEAKS=${ROOT}/data/peaks/K562_10M_CM_PU1_nan_nan_0_0_hg19/K562_10M_CM_PU1_nan_nan_0_0_hg19_peaks.narrowPeak
BAM=${ROOT}/data/mapped/${SAMPLE}.trimmed.bowtie2.shifted.dups.bam
GENOME=/media/afr/cemm-backup/reference/hg19/forBowtie2/hg19.fa
OUTPUT=${ROOT}/data/nucleoATAC/$SAMPLE
mkdir -p $OUTPUT

samtools index $BAM
pyatac vplot --cores 4 --lower 0 --upper 250 --flank 200 --plot_extra --bed $PEAKS --bam $BAM --out $OUTPUT

# pu1 motif centered
ROOT=/media/afr/cemm-backup/chipmentation
SAMPLE=K562_10M_CM_PU1_nan_PE_1_1_hg19
PEAKS=${ROOT}/data/peaks/K562_10M_CM_PU1_nan_nan_0_0_hg19/K562_10M_CM_PU1_nan_nan_0_0_hg19_peaks.motifCentered.bed
BAM=${ROOT}/data/mapped/${SAMPLE}.trimmed.bowtie2.shifted.dups.bam
GENOME=/media/afr/cemm-backup/reference/hg19/forBowtie2/hg19.fa
OUTPUT=${ROOT}/data/nucleoATAC/$SAMPLE.motifCentered
mkdir -p $OUTPUT

pyatac vplot --cores 4 --lower 0 --upper 250 --flank 200 --plot_extra --bed $PEAKS --bam $BAM --out $OUTPUT


# vplot 1k4 with top peaks
ROOT=/media/afr/cemm-backup/chipmentation
SAMPLE=K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19
PEAKS=${ROOT}/data/peaks/K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19/K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19_peaks.narrowPeak
PEAKS_TOP=${ROOT}/data/peaks/K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19/K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19_peaks.top.narrowPeak
BAM=${ROOT}/data/mapped/${SAMPLE}.trimmed.bowtie2.shifted.dups.bam
GENOME=/media/afr/cemm-backup/reference/hg19/forBowtie2/hg19.fa
OUTPUT=${ROOT}/data/nucleoATAC/${SAMPLE}.vplot
mkdir -p $OUTPUT

sort -k9nr $PEAKS | head -n 10000 > $PEAKS_TOP

pyatac vplot --cores 4 --lower 0 --upper 250 --flank 200 --plot_extra --bed $PEAKS_TOP --bam $BAM --out $OUTPUT



# vplot 1k4 on predicted nucleosomes
ROOT=/media/afr/cemm-backup/chipmentation
SAMPLE=K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19
PEAKS=${ROOT}/data/nucleoATAC/K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19.nfrpos.UNFINISHED.bed
BAM=${ROOT}/data/mapped/${SAMPLE}.trimmed.bowtie2.shifted.dups.bam
GENOME=/media/afr/cemm-backup/reference/hg19/forBowtie2/hg19.fa
OUTPUT=${ROOT}/data/nucleoATAC/${SAMPLE}.nucleosomes.UNFINISHED.vplot
mkdir -p $OUTPUT

pyatac vplot --cores 4 --lower 0 --upper 250 --flank 200 --plot_extra --bed $PEAKS --bam $BAM --out $OUTPUT


pyatac bias --fasta $GENOME


